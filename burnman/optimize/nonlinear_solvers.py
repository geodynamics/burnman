from __future__ import annotations
import numpy as np
import numpy.typing as npt
from typing import Any
from scipy.linalg import lu_factor, lu_solve
from types import SimpleNamespace
from collections import namedtuple


class DampedNewtonSolver:
    """
    Solver for the multivariate nonlinear system F(x) = 0 with Jacobian J(x),
    using the damped affine invariant modification to Newton's method
    (Deuflhard, 1974; 1975; 2004).

    This implementation follows the algorithm described in
    Nowak and Weimann (1991): Technical Report TR-91-10, Algorithm B,
    modified to accept linear inequality constraints.

    Linear inequality constraints are provided by a tuple argument named
    linear_constraints. The constraints are satisfied if A*x + b <= 0.
    If any constraints are violated by the current candidate solution using
    step size lambda, lambda is reduced to satisfy all constraints.

    If the current iterate lies on one or more constraints and the Newton step
    violates one or more of those constraints, the next step is calculated via
    the method of Lagrangian multipliers, minimizing the L2-norm of F(x_{i+1})
    subject to the violated constraints.

    Successful termination of the solver requires meeting all of the following
    criteria:
      - all(np.abs(dx (simplified Newton step)) < tol)
      - all(np.abs(dx (full Newton step)) < sqrt(10*tol))
      (to avoid pathological cases)
      - lambda = lambda_bounds(dx, x)[1] (lambda = 1 for a full Newton step)

    If these criteria are not satisfied, iterations continue until one of the
    following occurs:
      - lambda is reduced to its minimum value (indicating a very nonlinear problem)
      - successive iterations have descent vectors which violate the constraints
      - the maximum number of iterations is reached
    """

    def __init__(
        self,
        F,
        J,
        guess,
        tol: float = 1.0e-6,
        max_iterations: int = 100,
        lambda_bounds=lambda dx, x: (1.0e-8, 1.0),
        linear_constraints=(0.0, np.array([-1.0])),
        store_iterates: bool = False,
        regularization: float = 0.0,
        cond_lu_thresh: float = 1e12,
        cond_lstsq_thresh: float = 1e15,
    ):
        """
        Initialize the solver instance.

        :param F: Function that evaluates the system of nonlinear equations F(x) = 0
            at a given x.
        :type F: callable[[np.ndarray], np.ndarray]

        :param J: Function that evaluates the Jacobian matrix J(x) of F(x)
            at a given x.
        :type J: callable[[np.ndarray], np.ndarray]

        :param guess: Initial guess for the solution vector x.
        :type guess: np.ndarray

        :param tol: Convergence tolerance for Newton iterations, defaults to 1.0e-6.
        :type tol: float, optional

        :param max_iterations: Maximum number of Newton iterations to perform,
            defaults to 100.
        :type max_iterations: int, optional

        :param lambda_bounds: Function returning (min_lambda, max_lambda)
            for the damping factor given the current Newton step dx and point x,
            defaults to `lambda dx, x: (1.0e-8, 1.0)`.
        :type lambda_bounds: callable[[np.ndarray, np.ndarray], tuple[float, float]], optional

        :param linear_constraints: Tuple (A, b) representing linear inequality
            constraints A.x + b <= 0, defaults to (0.0, np.array([-1.0])).
        :type linear_constraints: tuple[np.ndarray, np.ndarray], optional

        :param store_iterates: If True, stores all intermediate iterates, function
            evaluations, and lambda values in the solution object, defaults to False.
        :type store_iterates: bool, optional

        :param regularization: Regularization parameter for the KKT system
            in Lagrangian solves, defaults to 0.0.
        :type regularization: float, optional

        :param cond_lu_thresh: Condition number threshold below which LU decomposition
            is considered stable, defaults to 1e12.
        :type cond_lu_thresh: float, optional

        :param cond_lstsq_thresh: Condition number threshold below which
            least-squares fallback is considered stable, defaults to 1e15.
        :type cond_lstsq_thresh: float, optional
        """
        self.F = F
        self.J = J
        self.guess = guess
        self.tol = tol
        self.max_iterations = max_iterations
        self.lambda_bounds = lambda_bounds
        self.linear_constraints = linear_constraints
        self.store_iterates = store_iterates
        self.regularization = regularization
        self.cond_lu_thresh = cond_lu_thresh
        self.cond_lstsq_thresh = cond_lstsq_thresh
        self.eps = 2.0 * np.finfo(float).eps

        assert np.all(
            self._constraints(self.guess) < self.eps
        ), "The starting guess is outside the supplied constraints."

        if not isinstance(self.tol, float):
            assert len(self.tol) < len(
                self.guess
            ), "tol must either be a float or an array like guess."

    def _constraints(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return np.dot(self.linear_constraints[0], x) + self.linear_constraints[1]

    def _update_lmda(
        self,
        x: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        h: npt.NDArray[np.float64],
        lmda_bounds: tuple[float, float],
    ) -> float:
        assert lmda_bounds[1] < 1.0 + self.eps
        assert lmda_bounds[0] > 1.0e-8 - self.eps

        lmda_j = min(1.0 / (h + self.eps), lmda_bounds[1])
        return max(lmda_j, lmda_bounds[0])

    def _solve_subject_to_constraints(
        self,
        x: npt.NDArray[np.float64],
        jac_x: npt.NDArray[np.float64],
        c_x: npt.NDArray[np.float64],
        c_prime: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """
        Solve a constrained Newton correction step using the method of
        Lagrange multipliers (KKT system).

        This method computes a step ``dx`` that minimizes the linearized
        residual ||J(x)·dx|| subject to linear equality constraints derived
        from the currently active inequality constraints.

        The system is solved using the KKT (Karush-Kuhn-Tucker) formulation:

        .. math::

            \\begin{bmatrix}
                J^T J + \\alpha I & A^T \\\\
                A & 0
            \\end{bmatrix}
            \\begin{bmatrix}
                dx \\\\
                \\lambda
            \\end{bmatrix}
            =
            - \\begin{bmatrix}
                0 \\\\
                c(x)
            \\end{bmatrix}

        where:

        * ``J`` is the Jacobian at ``x``
        * ``A`` is the constraint Jacobian (``c_prime``)
        * ``c(x)`` is the constraint evaluation
        * ``\\lambda`` are the Lagrange multipliers
        * ``\\alpha`` = ``self.regularization`` is an optional regularization parameter

        The KKT system is solved using one of three strategies depending on
        the estimated condition number of the matrix:

        1. **LU factorization** if ``cond < cond_lu_thresh``
        2. **Least-squares solve** if ``cond < cond_lstsq_thresh``
        3. **SVD-based pseudo-inverse** for ill-conditioned cases

        :param x: Current solution vector.
        :type x: np.ndarray
        :param jac_x: Current Jacobian matrix J(x).
        :type jac_x: np.ndarray
        :param c_x: Values of the active constraints at x.
        :type c_x: np.ndarray
        :param c_prime: Jacobian of the active constraints (A in Ax + b = 0).
        :type c_prime: np.ndarray

        :return: A 3-tuple containing:

            * **x_new** (np.ndarray) -- Updated solution ``x + dx``.
            * **lambdas** (np.ndarray) -- Lagrange multipliers for active constraints.
            * **condition_number** (float) -- Estimated condition number of the KKT matrix.

        :rtype: tuple[np.ndarray, np.ndarray, float]
        """
        n_x = x.shape[0]
        n_c = c_x.shape[0]
        JTJ_reg = jac_x.T @ jac_x + self.regularization * np.eye(n_x)
        norm = n_x * n_x / np.linalg.norm(JTJ_reg)
        KKT = np.block([[JTJ_reg * norm, c_prime.T], [c_prime, np.zeros((n_c, n_c))]])
        rhs = -np.concatenate([np.zeros(n_x), c_x])

        condition_number = np.linalg.cond(KKT)
        if condition_number < self.cond_lu_thresh:
            dx_lambda = lu_solve(lu_factor(KKT), rhs)
        elif condition_number < self.cond_lstsq_thresh:
            dx_lambda, *_ = np.linalg.lstsq(KKT, rhs, rcond=None)
        else:
            U, s, Vt = np.linalg.svd(KKT, full_matrices=False)
            s_inv = np.where(s > 1e-12, 1.0 / s, 0.0)
            dx_lambda = Vt.T @ (s_inv * (U.T @ rhs))

        dx = dx_lambda[:n_x]
        return x + dx, dx_lambda[n_x:], condition_number

    def _constrain_step_to_feasible_region(
        self,
        x: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        n_constraints: int,
        lmda: float,
        x_j: npt.NDArray[np.float64],
    ) -> tuple[npt.NDArray[np.float64], float]:
        """
        Project a trial Newton step back into the feasible region defined
        by linear inequality constraints A.x + b <= 0.

        This method checks whether the trial point x_j = x + lambda.dx violates
        any constraints. If so, it computes the maximum allowable step scaling
        factor to remain feasible, reduces lambda accordingly, and updates the
        trial iterate.

        The scaling factor is computed per violated constraint as:

        .. math::

            \\lambda_i = \\frac{c_x[i]}{c_x[i] - c_{x_j}[i]}

        where c_x and c_{x_j} are the constraint function values at x and x_j.
        The smallest lambda_i is used to rescale the step to just touch the first
        violated constraint.

        :param x: Current solution vector.
        :type x: np.ndarray
        :param dx: Full Newton step direction.
        :type dx: np.ndarray
        :param n_constraints: Total number of linear inequality constraints.
        :type n_constraints: int
        :param lmda: Current damping factor lambda for the trial step.
        :type lmda: float
        :param x_j: Current trial iterate x + lambda.dx.
        :type x_j: np.ndarray

        :return: A 3-tuple containing:

            * **lmda** (float)
              -- Updated damping factor lambda that ensures feasibility.
            * **x_j** (np.ndarray)
              -- Adjusted trial point within the feasible region.
            * **violated_constraints** (list[tuple[int, float]])
              -- List of (index, lambda_i) for each violated constraint,
              sorted by lambda_i.

        :rtype: tuple[float, np.ndarray, list[tuple[int, float]]]
        """
        c_x_j = self._constraints(x_j)
        c_x = self._constraints(x)
        violated_constraints = sorted(
            [
                (i, c_x[i] / (c_x[i] - c_x_j[i]))
                for i in range(n_constraints)
                if c_x_j[i] >= self.eps
            ],
            key=lambda x: x[1],
        )
        lmda *= violated_constraints[0][1]
        x_j = x + lmda * dx
        return lmda, x_j, violated_constraints

    def _lagrangian_walk_along_constraints(
        self,
        sol: Any,
        dx: npt.NDArray[np.float64],
        luJ: Any,
        dx_norm: float,
        violated_constraints: list[int],
    ) -> tuple[npt.NDArray[np.float64], float]:
        """
        Attempt to find a constrained Newton step when a step along the
        standard Newton direction would immediately violate active linear
        inequality constraints (A.x + b <= 0).
        Uses the method of Lagrange multipliers, attemping to "walk along"
        the active constraints to remain in the feasible region while
        decreasing the residual norm ||F(x)||.

        :param sol: Current solver state with fields x and F.
        :type sol: SimpleNamespace
        :param dx: Current Newton step direction.
        :type dx: np.ndarray
        :param luJ: LU factorization of the current Jacobian, as returned by
            ``scipy.linalg.lu_factor``.
        :type luJ: tuple
        :param dx_norm: L2 norm of the current Newton step dx.
        :type dx_norm: float
        :param lmda_bounds: Tuple (min_lambda, max_lambda) for the damping factor.
        :type lmda_bounds: tuple[float, float]
        :param violated_constraints: List of (index, fraction) for constraints
            that would be violated by the current Newton step.
        :type violated_constraints: list[tuple[int, float]]

        :return: Updated damping factor, updated values, full Newton step,
            and flag indicating whether the solver encountered a persistent
            constraint violation or reached the minimum lambda.
        :rtype: tuple[float, np.ndarray, np.ndarray, bool]
        """
        active_constraint_indices = [
            i for i, vc in violated_constraints if vc < self.eps
        ]
        inactive_constraint_indices = [
            i for i, vc in violated_constraints if vc >= self.eps
        ]
        c_newton = self._constraints(sol.x + dx)[active_constraint_indices]
        c_A = self.linear_constraints[0][active_constraint_indices]
        x_n = sol.x + dx
        persistent_bound_violation = False

        if len(c_A) > 0 and np.linalg.matrix_rank(c_A) == len(dx):
            n_act = len(active_constraint_indices)
            for i_rm in range(n_act):
                potential_active_indices = [
                    active_constraint_indices[i] for i in range(n_act) if i != i_rm
                ]
                c_newton = self._constraints(sol.x + dx)[potential_active_indices]
                c_A = self.linear_constraints[0][potential_active_indices]
                x_m = self._solve_subject_to_constraints(x_n, sol.J, c_newton, c_A)[0]
                if self._constraints(x_m)[active_constraint_indices[i_rm]] < 0.0:
                    break
        else:
            x_m = self._solve_subject_to_constraints(x_n, sol.J, c_newton, c_A)[0]

        dx = x_m - sol.x
        lmda_bounds_new = self.lambda_bounds(dx, sol.x)
        lmda = lmda_bounds_new[1]
        x_j = sol.x + lmda * dx

        # Check feasibility
        x_j_min = sol.x + lmda_bounds_new[0] * dx
        F_j_min = self.F(x_j_min)
        dxbar_j_min = lu_solve(luJ, -F_j_min)
        dxbar_j_min_norm = np.linalg.norm(dxbar_j_min, ord=2)

        if dxbar_j_min_norm > dx_norm or np.linalg.norm(dx, ord=2) < self.eps:
            persistent_bound_violation = True

        # Check newly violated inactive constraints
        n_inactive = len(inactive_constraint_indices)
        c_x_j = self._constraints(x_j)[inactive_constraint_indices]
        if not np.all(c_x_j < self.eps):
            c_x = self._constraints(sol.x)[inactive_constraint_indices]
            violated_constraints = sorted(
                [
                    (i, c_x[i] / (c_x[i] - c_x_j[i]))
                    for i in range(n_inactive)
                    if c_x_j[i] >= self.eps
                ],
                key=lambda x: x[1],
            )
            lmda *= violated_constraints[0][1]
            x_j = sol.x + lmda * dx

        return lmda, x_j, dx, persistent_bound_violation

    def _check_convergence(
        self,
        dxbar_j: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        lmda: float,
        lmda_bounds: tuple[float, float],
    ) -> bool:
        if (
            all(np.abs(dxbar_j) < self.tol)
            and all(np.abs(dx) < np.sqrt(10.0 * self.tol))
            and np.abs(lmda - lmda_bounds[1]) < self.eps
        ):
            return True
        return False

    def _posteriori_loop(
        self,
        x: npt.NDArray[np.float64],
        F: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        dx_norm: float,
        dxbar_j: npt.NDArray[np.float64],
        dxbar_j_norm: float,
        x_j: npt.NDArray[np.float64],
        luJ: Any,
        lmda: float,
        lmda_bounds: tuple[float, float],
        converged: bool,
        minimum_lmda: bool,
        persistent_bound_violation: bool,
        require_posteriori_loop: bool,
    ) -> tuple[npt.NDArray[np.float64], float, bool, bool, bool]:
        """
        Perform the a posteriori step-size control loop of Deuflhard's
        damped Newton method.

        After computing a trial step x_j = x + lambda.dx, this loop checks
        whether the simplified Newton step (dx̄_j) decreases the residual norm
        ||F(x)|| monotonically. If not, the damping factor lambda is reduced
        and the trial step is recomputed until either monotonicity is restored
        or the minimum lambda bound is reached. This procedure prevents
        divergence and stabilizes the Newton iteration in highly nonlinear
        regions.

        :param x: Current solution vector before update.
        :type x: np.ndarray
        :param F: Current function evaluation F(x).
        :type F: np.ndarray
        :param dx: Newton step vector.
        :type dx: np.ndarray
        :param dx_norm: Euclidean norm of dx.
        :type dx_norm: float
        :param dxbar_j: Simplified Newton step computed at the candidate next step.
        :type dxbar_j: np.ndarray
        :param dxbar_j_norm: Euclidean norm of dxbar_j.
        :type dxbar_j_norm: float
        :param x_j: Candidate next solution vector after applying step dx scaled by lambda.
        :type x_j: np.ndarray
        :param luJ: LU factorization of the Jacobian matrix at current x.
        :type luJ: tuple (lu, piv)
        :param lmda: Current damping factor (step length scaling).
        :type lmda: float
        :param lmda_bounds: Tuple (min_lambda, max_lambda) specifying allowed range of damping factors.
        :type lmda_bounds: tuple(float, float)
        :param converged: Flag indicating whether convergence criteria have been met.
        :type converged: bool
        :param minimum_lmda: Flag indicating whether the minimum lambda bound has been reached.
        :type minimum_lmda: bool
        :param persistent_bound_violation: Flag indicating if step violates constraints persistently.
        :type persistent_bound_violation: bool
        :param require_posteriori_loop: Flag indicating if the a posteriori loop should run.
        :type require_posteriori_loop: bool

        :return: Tuple containing updated values:
            - x (np.ndarray): Updated solution vector.
            - F (np.ndarray): Updated function evaluation.
            - dxbar (np.ndarray): Updated simplified Newton step.
            - dxprev (np.ndarray): Previous Newton step.
            - lmda (float): Updated damping factor.
            - minimum_lmda (bool): Updated minimum lambda flag.
            - converged (bool): Updated convergence flag.

        :rtype: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, bool, bool]
        """
        dxbar, dxprev = [1.0], [1.0]  # will be updated
        while (
            require_posteriori_loop
            and not minimum_lmda
            and not persistent_bound_violation
        ):
            if dxbar_j_norm <= dx_norm:
                if dxbar_j_norm < self.eps:
                    converged = True
                dxbar, x, F = dxbar_j, x_j, self.F(x_j)
                require_posteriori_loop = False
                dxprev = dx
            else:
                if np.abs(lmda - lmda_bounds[0]) < self.eps:
                    minimum_lmda = True
                h_j = (
                    (2.0 / lmda)
                    * np.linalg.norm((dxbar_j - (1.0 - lmda) * dx), ord=2)
                    / dx_norm
                )
                lmda_j = min(lmda_bounds[1], 1.0 / h_j)
                lmda = min(lmda_j, lmda / 2.0)
                lmda = max(lmda, lmda_bounds[0])

                x_j = x + lmda * dx
                F_j = self.F(x_j)
                dxbar_j = lu_solve(luJ, -F_j)
                dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

        return x, F, dxbar, dxprev, lmda, minimum_lmda, converged

    def _termination_info(
        self,
        converged: bool,
        minimum_lmda: bool,
        persistent_bound_violation: bool,
        lmda_bounds: tuple[float, float],
        n_it: int,
        max_iterations: int,
        violated_constraints: list[int],
    ) -> tuple[int, str, bool]:
        if converged:
            return (
                True,
                0,
                f"The solver successfully found a root after {n_it} iterations",
            )
        if minimum_lmda:
            return (
                False,
                1,
                f"The function is too non-linear for lower lambda bound ({lmda_bounds[0]})",
            )
        if persistent_bound_violation:
            return (
                False,
                2,
                (
                    "The descent vector crosses the constraints with the following indices: "
                    f"{[i for i, lmda in violated_constraints]}"
                ),
            )
        if n_it == max_iterations:
            return False, 3, f"The solver reached max_iterations ({max_iterations})"
        raise Exception("Unknown termination of solver")

    def solve(self) -> SimpleNamespace:
        """
        Execute the damped Newton solver to find a root of F(x) = 0,
        optionally subject to linear inequality constraints A·x + b <= 0.

        If the solver fails to converge, it terminates with an informative
        status code and description.

        :return: Solution object with the following attributes:

            * **x** (np.ndarray) -- Final solution vector.
            * **F** (np.ndarray) -- Function evaluation at the solution, F(x).
            * **F_norm** (float) -- Euclidean (L2) norm of F(x).
            * **J** (np.ndarray) -- Jacobian matrix evaluated at the solution.
            * **success** (bool) -- True if the solver converged.
            * **code** (int) -- Termination code:

                * 0 -- Successful convergence
                * 1 -- Failure: lambda reached its minimum bound
                * 2 -- Failure: descent direction violates constraints
                * 3 -- Failure: maximum iterations reached

            * **text** (str) -- Human-readable termination description.
            * **n_it** (int) -- Number of Newton iterations performed.
            * **iterates** (namedtuple, optional) -- Present if
              ``store_iterates=True``. Contains iteration history:

                * **x** (list[np.ndarray]) -- Solution for each iteration.
                * **F** (list[np.ndarray]) -- Function evaluation for each iteration.
                * **lmda** (list[float]) -- Step length scaling factor for each iteration.

        :rtype: SimpleNamespace
        """
        sol = namedtuple(
            "Solution", ["x", "n_it", "F", "F_norm", "J", "code", "text", "success"]
        )
        sol.x = self.guess
        sol.F = self.F(sol.x)

        if self.store_iterates:
            sol.iterates = namedtuple("iterates", ["x", "F", "lmda"])
            sol.iterates.x, sol.iterates.F, sol.iterates.lmda = [sol.x], [sol.F], [0.0]

        lmda, dxprev, dxbar = 0.0, [1.0], [1.0]
        sol.n_it = 0
        n_constraints = len(self._constraints(sol.x))
        minimum_lmda = converged = persistent_bound_violation = False

        while (
            sol.n_it < self.max_iterations
            and not minimum_lmda
            and not persistent_bound_violation
            and not converged
        ):
            sol.J = self.J(sol.x)
            luJ = lu_factor(sol.J)
            dx = lu_solve(luJ, -sol.F)
            dx_norm = np.linalg.norm(dx, ord=2)

            lmda_bounds = self.lambda_bounds(dx, sol.x)
            h = (
                lmda
                * np.linalg.norm((dxbar - dx), ord=2)
                * dx_norm
                / (np.linalg.norm(dxprev, ord=2) * np.linalg.norm(dxbar, ord=2))
            )
            lmda = self._update_lmda(sol.x, dx, h, lmda_bounds)

            x_j = sol.x + lmda * dx
            c_x_j = self._constraints(x_j)

            if not np.all(c_x_j < self.eps):
                lmda, x_j, violated_constraints = (
                    self._constrain_step_to_feasible_region(
                        sol.x, dx, n_constraints, lmda, x_j
                    )
                )

            if lmda < self.eps:
                lmda, x_j, dx, persistent_bound_violation = (
                    self._lagrangian_walk_along_constraints(
                        sol, dx, luJ, dx_norm, violated_constraints
                    )
                )

            # Evaluate simplified Newton step
            F_j = self.F(x_j)
            dxbar_j = lu_solve(luJ, -F_j)
            dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

            converged = self._check_convergence(dxbar_j, dx, lmda, lmda_bounds)
            require_posteriori_loop = not converged

            loop_vars = self._posteriori_loop(
                sol.x,
                sol.F,
                dx,
                dx_norm,
                dxbar_j,
                dxbar_j_norm,
                x_j,
                luJ,
                lmda,
                lmda_bounds,
                converged,
                minimum_lmda,
                persistent_bound_violation,
                require_posteriori_loop,
            )

            sol.x, sol.F, dxbar, dxprev, lmda, minimum_lmda, converged = loop_vars

            sol.n_it += 1
            if self.store_iterates:
                sol.iterates.x.append(sol.x)
                sol.iterates.F.append(sol.F)
                sol.iterates.lmda.append(lmda)

        # Final adjustment for constraints
        if converged and not persistent_bound_violation:
            sol.x = x_j + dxbar_j
            if not np.all(self._constraints(sol.x) <= 0.0):
                sol.x -= dxbar_j

        sol.F = self.F(sol.x)
        sol.F_norm = np.linalg.norm(sol.F, ord=2)
        sol.J = self.J(sol.x)

        if self.store_iterates:
            sol.iterates.x = np.array(sol.iterates.x)
            sol.iterates.F = np.array(sol.iterates.F)

        sol.success, sol.code, sol.text = self._termination_info(
            converged,
            minimum_lmda,
            persistent_bound_violation,
            lmda_bounds,
            sol.n_it,
            self.max_iterations,
            locals().get("violated_constraints", []),
        )
        return sol


def damped_newton_solve(
    F,
    J,
    guess,
    tol: float = 1.0e-6,
    max_iterations: int = 100,
    lambda_bounds=lambda dx, x: (1.0e-8, 1.0),
    linear_constraints=(0.0, np.array([-1.0])),
    store_iterates: bool = False,
) -> SimpleNamespace:
    """
    Helper function to run the DampedNewtonSolver.

    :param F: Function that evaluates the system of nonlinear equations F(x) = 0
        at a given x.
    :type F: callable[[np.ndarray], np.ndarray]

    :param J: Function that evaluates the Jacobian matrix J(x) of F(x)
        at a given x.
    :type J: callable[[np.ndarray], np.ndarray]

    :param guess: Initial guess for the solution vector x.
    :type guess: np.ndarray

    :param tol: Convergence tolerance for Newton iterations, defaults to 1.0e-6.
    :type tol: float, optional

    :param max_iterations: Maximum number of Newton iterations to perform,
        defaults to 100.
    :type max_iterations: int, optional

    :param lambda_bounds: Function returning (min_lambda, max_lambda)
        for the damping factor given the current Newton step dx and point x,
        defaults to `lambda dx, x: (1.0e-8, 1.0)`.
    :type lambda_bounds: callable[[np.ndarray, np.ndarray], tuple[float, float]], optional

    :param linear_constraints: Tuple (A, b) representing linear inequality
        constraints A.x + b <= 0, defaults to (0.0, np.array([-1.0])).
    :type linear_constraints: tuple[np.ndarray, np.ndarray], optional

    :param store_iterates: If True, stores all intermediate iterates, function
        evaluations, and lambda values in the solution object, defaults to False.
    :type store_iterates: bool, optional


    :rtype: SimpleNamespace
    """
    solver = DampedNewtonSolver(
        F,
        J,
        guess,
        tol,
        max_iterations,
        lambda_bounds,
        linear_constraints,
        store_iterates,
    )
    return solver.solve()
