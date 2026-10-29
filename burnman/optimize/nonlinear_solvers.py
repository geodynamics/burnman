# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.
import numpy as np
import numpy.typing as npt
from typing import Any, Optional
from scipy.linalg import lu_factor, lu_solve
from enum import Enum
from dataclasses import dataclass, field


@dataclass
class Iterates:
    """
    Container for storing iteration history of the solver.
    Attributes are lists of the iterates x, function evaluations F,
    and step lengths relative to the full Newton step lmda.
    """

    x: list = field(default_factory=list)
    F: list = field(default_factory=list)
    lmda: list = field(default_factory=list)

    def append(self, x, F, lmda):
        """
        Append an iteration's data to the history.
        :param x: Current iterate.
        :param F: Current function evaluation.
        :param lmda: Current step length.
        """
        self.x.append(x)
        self.F.append(F)
        self.lmda.append(lmda)

    def close(self):
        """
        Finalize the iteration history by converting lists to arrays.
        """
        self.x = np.array(self.x)
        self.F = np.array(self.F)
        self.lmda = np.array(self.lmda)

    def __repr__(self):
        """
        Human-readable string representation of the Iterates instance.
        """
        s = "Iterates:\n"
        for i in range(len(self.x)):
            s += f"Iterate {i}: x={self.x[i]}, F={self.F[i]}, lmda={self.lmda[i]}\n"
        return s


class TerminationCode(Enum):
    """
    Enumeration of termination codes for the nonlinear solver.
    These are used to indicate the reason for solver termination.
    The codes are:

        - SUCCESS (0): The solver successfully found a root.
        - TOO_NONLINEAR (1): The function is too non-linear for lower lambda bound.
        - CONSTRAINT_VIOLATION (2): The descent vector crosses the constraints.
        - MAX_ITERATIONS (3): The solver reached the maximum number of iterations.
        - SINGULAR_FAIL (4): The system was effectively singular and failed to converge.
        - SINGULAR_SUCCESS (5): The solver found a root but the system is singular.

    """

    SUCCESS = 0
    TOO_NONLINEAR = 1
    CONSTRAINT_VIOLATION = 2
    MAX_ITERATIONS = 3
    SINGULAR_FAIL = 4
    SINGULAR_SUCCESS = 5


@dataclass
class Solution:
    """
    Container for the solution of a nonlinear system of equations.

    :param x: Final solution vector.
    :type x: np.ndarray, optional
    :param n_it: Number of iterations performed.
    :type n_it: int, optional, default is 0
    :param F: Function evaluation at the solution.
    :type F: np.ndarray, optional
    :param F_norm: Euclidean norm of F at the solution.
    :type F_norm: float, optional
    :param J: Jacobian matrix at the solution.
    :type J: np.ndarray, optional
    :param code: Termination code. Valid codes and their values
        are given by the TerminationCode enum.
    :type code: TerminationCode, optional
    :param text: Human-readable termination description.
    :type text: str, optional
    :param success: True if the solver converged.
    :type success: bool, optional, default is False
    :param store_iterates: If True, iteration history is stored in an
        `iterates` attribute. The Iterates instance contains lists
        of all iterates x, function evaluations F, and step lengths
        relative to the full Newton step lmda.
    :type store_iterates: bool, optional, default is False
    """

    x: Optional[np.ndarray] = None
    n_it: int = 0
    F: Optional[np.ndarray] = None
    F_norm: Optional[float] = None
    J: Optional[np.ndarray] = None
    code: Optional[TerminationCode] = None
    text: Optional[str] = None
    success: bool = False
    store_iterates: bool = False

    def __post_init__(self):
        if self.store_iterates:
            self.iterates = Iterates()

    def __repr__(self):
        """
        Human-readable string representation of the Solution instance.
        """
        s = f"{self.message}\n"
        s += f"x = {self.x}\n"
        s += f"F = {self.F}\n"
        s += f"F_norm = {self.F_norm}\n"
        s += f"J = {self.J}\n"
        try:
            s += f"{self.iterates}"
        except AttributeError:
            pass
        return s

    def __str__(self):
        """
        String representation of the Solution instance.
        """
        return self.__repr__()

    def terminate(
        self,
        converged: bool,
        minimum_lmda: bool,
        persistent_bound_violation: bool,
        lmda_bounds: tuple[float, float],
        n_it: int,
        max_iterations: int,
        violated_constraints: list[int],
        singular_jacobian: bool,
        condition_number: float,
    ) -> None:
        """
        Set the termination status of the solver.
        """
        if self.store_iterates:
            self.iterates.close()

        self.success = True if converged else False

        if converged and not singular_jacobian:
            self.code = TerminationCode.SUCCESS
            self.message = (
                f"The solver successfully found a root after {n_it} iterations"
            )
        elif minimum_lmda:
            self.code = TerminationCode.TOO_NONLINEAR
            self.message = f"The function is too non-linear for lower lambda bound ({lmda_bounds[0]})"
        elif persistent_bound_violation and not singular_jacobian:
            self.code = TerminationCode.CONSTRAINT_VIOLATION
            self.message = (
                "The descent vector crosses the constraints with the following indices: "
                f"{[i for i, lmda in violated_constraints]}"
            )
        elif n_it == max_iterations:
            self.code = TerminationCode.MAX_ITERATIONS
            self.message = f"The solver reached max_iterations ({max_iterations})"
        elif singular_jacobian and not converged:
            self.code = TerminationCode.SINGULAR_FAIL
            self.message = (
                f"The system was effectively singular (cond={condition_number:.2e}) "
                f"and tried to leave the feasible region at iteration {n_it} by "
                "crossing the constraints with the following indices: "
                f"{[i for i, _ in violated_constraints]}.\n"
                "You may want to try a different initial guess "
                "or check the Jacobian implementation. "
                "It may be that the problem is ill-posed."
            )
        elif singular_jacobian and converged:
            self.code = TerminationCode.SINGULAR_SUCCESS
            self.message = (
                f"The solver successfully found a root after {n_it} iterations,"
                f"but the system is effectively singular (cond={condition_number:.2e}). "
                "You should check the problem formulation before trusting this result."
            )
        else:
            raise Exception("Unrecognised termination condition for nonlinear solver")

        return None


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
    Note that the sign of b is opposite to that used in
    scipy.optimize.LinearConstraint. If any constraints are violated by the
    current candidate solution using step size lambda, lambda is reduced to
    satisfy all constraints.

    If the current iterate lies on one or more constraints and the Newton step
    violates one or more of those constraints, the next step is calculated via
    the method of Lagrangian multipliers, minimizing the L2-norm of F(x_{i+1})
    subject to the violated constraints.

    Successful termination of the solver requires meeting all of the following
    criteria:

        - all(np.abs(dx (simplified Newton step)) < tol)
        - all(np.abs(dx (full Newton step)) < np.sqrt(10*tol))
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
        tol: float | npt.NDArray = 1.0e-6,
        F_tol: float | npt.NDArray = 1.0e-8,
        max_iterations: int = 100,
        lambda_bounds=lambda dx, x: (1.0e-8, 1.0),
        linear_constraints=(0.0, np.array([-1.0])),
        store_iterates: bool = False,
        regularization: float = np.finfo(float).eps,
        cond_lu_thresh: float = 1e12,
        constraint_thresh: float = 2 * np.finfo(float).eps,
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

        :param tol: Convergence tolerance for each component of the
            Newton step dx, defaults to 1.0e-6.
        :type tol: float, or numpy.ndarray, optional

        :param F_tol: Convergence tolerance for each component of F(x),
            defaults to 1.0e-8.
        :type F_tol: float, or numpy.ndarray, optional

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

        :param regularization: Regularization parameter for the Jacobian.
            This parameter scales a diagonal matrix, which is added to the Jacobian
            when its condition number exceeds cond_lu_thresh,
            or by default for solves subject to one or more linear constraints.
            This helps stabilize the solution of the linear system in ill-conditioned cases.
            Must be a small positive value; defaults to numpy float epsilon.
        :type regularization: float, optional

        :param cond_lu_thresh: Condition number threshold below which LU decomposition
            is considered stable, defaults to 1e12.
        :type cond_lu_thresh: float, optional

        :param cond_lstsq_thresh: Condition number threshold below which
            least-squares fallback is considered stable, defaults to 1e15.
        :type cond_lstsq_thresh: float, optional

        :param constraint_thresh: Threshold for considering a constraint
            "active" when determining step feasibility, defaults to 2*eps.
        :type constraint_thresh: float, optional
        """

        self.sol = Solution(x=guess, F=F(guess), store_iterates=store_iterates)
        self.F = F
        self.J = J
        self.tol = tol
        self.F_tol = F_tol
        self.constraint_thresh = constraint_thresh
        self.max_iterations = max_iterations
        self.lambda_bounds = lambda_bounds
        self.linear_constraints = linear_constraints
        self.regularization = regularization
        self.cond_lu_thresh = cond_lu_thresh
        self.eps = 2.0 * np.finfo(float).eps
        self.max_condition_number = 1.0 / np.finfo(float).eps

        assert np.all(
            self._constraints(guess) < self.eps
        ), "The starting guess is outside the supplied constraints."

        if not isinstance(self.tol, float):
            self.tol = np.asarray(self.tol, dtype=np.float64)
            assert self.tol.ndim == 1, "tol must be a float or a 1D array."
            assert len(self.tol) == len(
                guess
            ), f"tol must either be a float or an array like guess. Got len(tol)={len(self.tol)} and len(guess)={len(guess)}."

    def _constraints(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return np.dot(self.linear_constraints[0], x) + self.linear_constraints[1]

    def _update_lmda(
        self,
        x: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        h: npt.NDArray[np.float64],
        lmda_bounds: tuple[float, float],
    ) -> float:
        """
        Update the damping factor lambda based on the current Newton step.
        """
        assert (
            lmda_bounds[1] < 1.0 + self.eps
        ), f"max_lambda must be <= 1.0, got {lmda_bounds[1]}"
        assert (
            lmda_bounds[0] > 1.0e-8 - self.eps
        ), f"min_lambda must be > 1.0e-8, got {lmda_bounds[0]}"

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
        Compute a constrained Newton correction step using the
        Karush-Kuhn-Tucker (KKT) formulation.

        This method solves for the update ``dx`` that minimizes the linearized
        residual ``||J(x).dx||`` subject to the active linear equality constraints
        ``A.dx + c(x) = 0``:

            [ J^TJ + alpha * I   c'^T ] [dx]   =  -[ 0 ]
            [    c'        0  ] [lambda]         [c(x)]

        where:
            * ``J`` is the Jacobian at ``x`` (``jac_x``)
            * ``c_prime`` is the active constraint Jacobian
            * ``c(x)`` are the constraint values
            * ``lambda`` are the Lagrange multipliers
            * ``alpha = self.regularization`` is a regularization parameter

        The KKT system is solved adaptively based on its condition number:
            1. LU factorization for well-conditioned systems
            2. SVD-based pseudo-inverse for ill-conditioned systems

        :param x: Current solution vector.
        :type x: numpy.ndarray
        :param jac_x: Jacobian of residuals at ``x``.
        :type jac_x: numpy.ndarray
        :param c_x: Values of the active constraints at ``x``.
        :type c_x: numpy.ndarray
        :param c_prime: Jacobian of the active constraints.
        :type c_prime: numpy.ndarray

        :return: Tuple ``(x_new, lambdas, condition_number)`` where:
          - **x_new** (*numpy.ndarray*) - Updated solution vector ``x + dx``.
          - **lambdas** (*numpy.ndarray*) - Lagrange multipliers for active constraints.
          - **condition_number** (*float*) - Estimated condition number of the KKT matrix.
        :rtype: tuple[numpy.ndarray, numpy.ndarray, float]
        """

        n_x = x.shape[0]
        n_c = c_x.shape[0]
        JTJ_reg = jac_x.T @ jac_x + self.regularization * np.eye(n_x)
        scale = np.linalg.norm(JTJ_reg)
        KKT = np.block([[JTJ_reg / scale, c_prime.T], [c_prime, np.zeros((n_c, n_c))]])
        rhs = -np.concatenate([np.zeros(n_x), c_x])

        condition_number = np.linalg.cond(KKT)
        if condition_number < self.cond_lu_thresh:
            dx_lambda = lu_solve(lu_factor(KKT), rhs)
        else:
            U, s, Vt = np.linalg.svd(KKT, full_matrices=False)
            tol = self.eps * max(KKT.shape) * np.max(s)
            s_inv = np.where(s > tol, 1.0 / s, 0.0)
            dx_lambda = Vt.T @ (s_inv * (U.T @ rhs))

        dx = dx_lambda[:n_x]
        lambdas = dx_lambda[n_x:]

        return x + dx, lambdas, condition_number

    def _constrain_step_to_feasible_region(
        self,
        x: npt.NDArray[np.float64],
        dx: npt.NDArray[np.float64],
        n_constraints: int,
        lmda: float,
        x_j: npt.NDArray[np.float64],
    ) -> tuple[npt.NDArray[np.float64], float]:
        """
        Project a trial step back into the feasible region of linear inequality constraints.

        Given a trial point x_j = x + lambda*dx, this method checks for constraint
        violations and rescales the step to remain feasible. The scaling factor is
        computed per violated constraint, and the smallest factor is applied to
        lambda to ensure the trial point stays within the feasible region.

        :param x: Current solution vector.
        :type x: numpy.ndarray
        :param dx: Newton step direction.
        :type dx: numpy.ndarray
        :param n_constraints: Number of linear inequality constraints.
        :type n_constraints: int
        :param lmda: Current step scaling factor.
        :type lmda: float
        :param x_j: Trial point x + lambda*dx.
        :type x_j: numpy.ndarray

        :return: Tuple ``(lmda, x_j, violated_constraints)`` where:
          - **lmda** (*float*) - Updated scaling factor ensuring feasibility.
          - **x_j** (*numpy.ndarray*) - Adjusted trial point within feasible region.
          - **violated_constraints** (*list[tuple[int, float]]*) - List of
            (constraint index, scaling factor) for violated constraints, sorted by factor.
        :rtype: tuple[float, numpy.ndarray, list[tuple[int, float]]]
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
        dx: npt.NDArray[np.float64],
        luJ: Any,
        dx_norm: float,
        violated_constraints: list[tuple[int, float]],
    ) -> tuple[float, npt.NDArray[np.float64], npt.NDArray[np.float64], bool]:
        """
        Attempt a constrained Newton step along active linear constraints
        to remain feasible while decreasing the residual norm.

        :param dx: Newton step direction.
        :param luJ: LU factorization of current Jacobian (from `lu_factor`).
        :param dx_norm: L2 norm of the Newton step.
        :param violated_constraints: List of (index, fraction) for constraints
            that would be violated by the current step.

        :return: Tuple of (lambda, adjusted trial point x_j, full Newton step dx,
                persistent_bound_violation flag).
        """
        sol = self.sol

        # Split constraints into active and inactive based on proximity to boundary
        active_idx = [
            i for i, vc in violated_constraints if vc < self.constraint_thresh
        ]
        inactive_idx = [
            i for i, vc in violated_constraints if vc >= self.constraint_thresh
        ]

        # Evaluate active constraints and corresponding Jacobian
        c_active = self._constraints(sol.x + dx)[active_idx]
        A_active = self.linear_constraints[0][active_idx]
        x_n = sol.x + dx
        persistent_violation = False

        # Solve KKT system along active constraints if well-posed
        if len(A_active) > 0 and np.linalg.matrix_rank(A_active) == len(dx):
            n_act = len(active_idx)
            # Attempt to remove one active constraint at a time if necessary
            for i_rm in range(n_act):
                keep_idx = [active_idx[j] for j in range(n_act) if j != i_rm]
                c_subset = self._constraints(sol.x + dx)[keep_idx]
                A_subset = self.linear_constraints[0][keep_idx]
                x_m = self._solve_subject_to_constraints(
                    x_n, sol.J, c_subset, A_subset
                )[0]
                if self._constraints(x_m)[active_idx[i_rm]] < 0:
                    break
        else:
            x_m = self._solve_subject_to_constraints(x_n, sol.J, c_active, A_active)[0]

        # Update step and damping factor
        dx = x_m - sol.x
        lmda_min, lmda_max = self.lambda_bounds(dx, sol.x)
        lmda = lmda_max
        x_j = sol.x + lmda * dx

        # Check feasibility at minimum lambda
        try:
            x_j_min = sol.x + lmda_min * dx
            F_j_min = self.F(x_j_min)
            dxbar_j_min = lu_solve(luJ, -F_j_min)
            if np.linalg.norm(dxbar_j_min) > dx_norm or np.linalg.norm(dx) < self.eps:
                persistent_violation = True
        except Exception:
            persistent_violation = True

        # Check that inactive constraints are not now violated
        # If they are, rescale lambda
        if inactive_idx:
            c_inactive_new = self._constraints(x_j)[inactive_idx]
            if not np.all(c_inactive_new < self.eps):
                c_inactive_old = self._constraints(sol.x)[inactive_idx]
                violated_new = sorted(
                    [
                        (i, c_inactive_old[i] / (c_inactive_old[i] - c_inactive_new[i]))
                        for i in range(len(inactive_idx))
                        if c_inactive_new[i] >= self.eps
                    ],
                    key=lambda t: t[1],
                )
                # Rescale lambda to maintain feasibility
                lmda *= violated_new[0][1]
                x_j = sol.x + lmda * dx

        return lmda, x_j, dx, persistent_violation

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

    def solve(self) -> Solution:
        """
        Execute the damped Newton solver to find a root of F(x) = 0,
        optionally subject to linear inequality constraints A.x + b <= 0.

        If the solver fails to converge, it terminates with an informative
        status code and description.

        :return: Solution object. See `Solution` class for details of its attributes.
        :rtype: Solution instance
        """
        sol = self.sol

        if sol.store_iterates:
            sol.iterates.append(sol.x, sol.F, 0.0)

        # Initialize variables for the optimization loop
        lmda, dxprev, dxbar = 0.0, [1.0], [1.0]
        n_constraints = len(self._constraints(sol.x))
        minimum_lmda = converged = persistent_bound_violation = False

        while (
            sol.n_it < self.max_iterations
            and not minimum_lmda
            and not persistent_bound_violation
            and not converged
        ):
            sol.J = self.J(sol.x)
            condition_number = np.linalg.cond(sol.J)

            # Regularize ill-conditioned Jacobian
            if condition_number > self.cond_lu_thresh:
                sol.J = sol.J + np.eye(sol.J.shape[0]) * self.regularization

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
                        dx, luJ, dx_norm, violated_constraints
                    )
                )

            # Evaluate simplified Newton step
            F_j = self.F(x_j)
            dxbar_j = lu_solve(luJ, -F_j)
            dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

            # Additional convergence check on F(x)
            converged = self._check_convergence(dxbar_j, dx, lmda, lmda_bounds)
            if converged and not all(np.abs(F_j) < self.F_tol):
                converged = False

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
            if sol.store_iterates:
                sol.iterates.append(sol.x, sol.F, lmda)

        # Final adjustment for constraints
        # and recompute F, J, condition number without regularization
        if condition_number < self.max_condition_number:
            if converged and not persistent_bound_violation:
                sol.x = x_j + dxbar_j
                if not np.all(self._constraints(sol.x) <= 0.0):
                    sol.x -= dxbar_j

            sol.F = self.F(sol.x)
            sol.F_norm = np.linalg.norm(sol.F, ord=2)
            sol.J = self.J(sol.x)
            condition_number = np.linalg.cond(sol.J)

        sol.terminate(
            converged=converged,
            minimum_lmda=minimum_lmda,
            persistent_bound_violation=persistent_bound_violation,
            lmda_bounds=lmda_bounds,
            n_it=sol.n_it,
            max_iterations=self.max_iterations,
            violated_constraints=locals().get("violated_constraints", []),
            singular_jacobian=condition_number >= self.max_condition_number,
            condition_number=condition_number,
        )

        return sol


def damped_newton_solve(
    F,
    J,
    guess,
    tol: float | npt.NDArray = 1.0e-6,
    F_tol: float | npt.NDArray = 1.0e-8,
    max_iterations: int = 100,
    lambda_bounds=lambda dx, x: (1.0e-8, 1.0),
    linear_constraints=(0.0, np.array([-1.0])),
    store_iterates: bool = False,
) -> Solution:
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

    :return: Solution object. See `Solution` class for details of its attributes.
    :rtype: Solution instance
    """
    solver = DampedNewtonSolver(
        F,
        J,
        guess,
        tol,
        F_tol,
        max_iterations,
        lambda_bounds,
        linear_constraints,
        store_iterates,
    )
    return solver.solve()
