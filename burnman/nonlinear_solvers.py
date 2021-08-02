import numpy as np
from scipy.linalg import lu_factor, lu_solve
from collections import namedtuple


def solve_constraint_lagrangian(x, jac_x, c_x, c_prime):
    """
    Function which solves the problem
    minimize || J.dot(x_mod - x) ||
    subject to C(x_mod) = 0
    via the method of Lagrange multipliers.

    Parameters
    ----------
    x : 1D numpy array
        Parameter values at x
    jac_x : 2D numpy array.
        The (estimated, approximate or exact)
        value of the Jacobian J(x)
    c_x : 1D numpy array
        Values of the constraints at x
    c_prime : 2D array of floats
        The Jacobian of the constraints
        (A, where A.x + b = 0)

    Returns
    -------
    x_mod : 1D numpy array
        The parameter values which minimizes the L2-norm
        of any function which has the Jacobian jac_x.
    lagrange_multipliers : 1D numpy array
        The multipliers for each of the equality
        constraints
    """
    n_x = len(x)
    n = n_x + len(c_x)
    A = np.zeros((n, n))
    b = np.zeros(n)

    JTJ = jac_x.T.dot(jac_x)
    A[:n_x, :n_x] = JTJ/np.linalg.norm(JTJ)*n*n  # includes scaling
    A[:n_x, n_x:] = c_prime.T
    A[n_x:, :n_x] = c_prime
    b[n_x:] = c_x

    luA = lu_factor(A)
    dx_m = lu_solve(luA, -b)  # lu_solve computes the solution of ax = b

    x_mod = x + dx_m[:n_x]
    lagrange_multipliers = dx_m[n_x:]
    return (x_mod, lagrange_multipliers)


def damped_newton_solve(F, J, guess, tol=1.e-6,
                        max_iterations=100,
                        lambda_bounds=lambda dx, x: (1.e-8, 1.),
                        linear_constraints=(0., np.array([-1.])),
                        store_iterates=False):
    """
    Solver for the multivariate nonlinear system F(x)=0
    with Jacobian J(x), using the damped affine invariant modification
    to Newton's method (Deuflhard, 1974;1975;2004).
    Here we follow the algorithm as described in Nowak and Weimann (1991):
    [Technical Report TR-91-10, Algorithm B], modified to accept
    linear inequality constraints.

    Linear inequality constraints are provided by the arrays constraints_A and
    constraints_b. The constraints are satisfied if A*x + b <= 0.
    If any constraints are not satisfied by the current
    value of lambda, lambda is reduced to satisfy all the constraints.

    If a current iterate starting point (x_i) lies on one or more constraints
    and the Newton step violates one or more of those constraints, then
    the next step is calculated via the method of Lagrangian multipliers,
    minimizing the L2-norm of F(x_i+1) subject to the violated constraints.

    Successful termination of the solver is based on three criteria:
    - all(np.abs(dx (simplified newton step) < tol))
    - all(np.abs(dx (full Newton step) < sqrt(10*tol))) [avoiding pathological behaviour] and
    - lambda = lambda_bounds(dx, x)[1] (lambda = 1 for a full Newton step).

    If these criteria are not satisfied, iterations continue until one of the following
    occurs:
    - the value of lmda is reduced to its minimum value
      (this happens if the problem is very nonlinear)
    - successive iterations have descent vectors which violate the constraints
    - the maximum number of iterations (given by max_iterations) is reached.

    Information on the root (or lack of root) obtained by the solver is provided
    in the returned namedtuple.


    Parameters
    ----------
    F : function of x
        Returns the system function F(x)
        as a 1D numpy array.
    J : function of x
        Returns the Jacobian function J(x)
        as a 2D numpy array.
    guess : 1D numpy array
        Starting guess for the solver.
    tol : float or array of floats [1.e-6]
        Tolerance(s) for termination.
    max_iterations : integer [100]
        Maximum number of iterations for solver.
    lambda_bounds: function of dx and x
        Returns a tuple of floats (1.e-8, 1.) corresponding
        to the minimum and maximum allowed fractions of the
        full newton step (dx).
    linear_constraints : tuple of a 2D numpy array (A) and 1D numpy array (b)
        Constraints are satisfied if A.x + b < eps


    Returns
    -------
    sol : namedtuple
        Includes the following attributes:
        x : 1D numpy array of floats
            The solution vector.
        F : 1D numpy array of floats
            The evaluated function F(x).
        F_norm : float
            Euclidean norm of F(x).
        J : 2D numpy array of floats
            The evaluated Jacobian J(x).
        n_it : integer
            Number of iterations.
        code : integer
            Numerical description of the solver termination.
                0 -> Successful convergence
                1 -> Failure due to solver hitting lower lambda bound
                2 -> Failure due to descent vector crossing constraints
                3 -> Failure due to solver reaching maximum number of iterations
        text : string
            Description of the solver termination.
        success : bool
            Solution convergence boolean.
        iterates : namedtuple
            Only present if store_iterates=True
            Includes the following attributes:
            x : list of 1D numpy arrays of floats
                The parameters for each iteration
            F : list of 2D numpy arrays of floats
                The function for each iteration
            lmda : list of floats
                The value of the damping parameter for each iteration

    This function is available as ``burnman.damped_newton_solve``.
    """

    # Make sure damping factor is within bounds, and that the bounds are reasonable
    # Problem classes in Nowak and Weimann (1991); [lmda_min, lmda_max]:
    # linear: [0.1, 1.]
    # mildly nonlinear: [1.e-4, 1.]
    # highly nonlinear: [1.e-2, 1.e-4]
    # extremely nonlinear: [1.e-4, 1.e-8]
    eps = 2.*np.finfo(float).eps

    def update_lmda(x, dx, h, lmda_bounds):
        assert lmda_bounds[1] < 1. + eps, 'The highest upper bound for lambda is 1. (a full Newton step)'
        assert lmda_bounds[0] > 1.e-8 - eps, 'The lowest lower bound for lambda is 1.e-8 (suitable only for extremely nonlinear systems)'

        lmda_j = min(1./(h + eps), lmda_bounds[1])  # this is lmda_j^0
        return max(lmda_j, lmda_bounds[0])

    def constraints(x):
        return np.dot(linear_constraints[0], x) + linear_constraints[1]

    assert np.all(constraints(guess) < eps), 'The starting guess is outside the supplied constraints.'

    if not isinstance(tol, float):
        assert len(tol) < len(guess), 'tol must either be a float or an array like guess.'

    sol = namedtuple('Solution', ['x', 'n_it', 'F', 'F_norm', 'J', 'code', 'text', 'success'])

    # evaluate system
    sol.x = guess
    sol.F = F(sol.x)

    if store_iterates:
        sol.iterates = namedtuple('iterates', ['x', 'F', 'lmda'])
        sol.iterates.x = [sol.x]
        sol.iterates.F = [sol.F]
        sol.iterates.lmda = [0.]

    # Begin Newton loop

    # Some dummy variables for the first h calculation (h = 0)
    lmda = 0.
    dxprev = [1.]
    dxbar = [1.]

    sol.n_it = 0
    n_constraints = len(constraints(sol.x))
    minimum_lmda = False
    converged = False
    persistent_bound_violation = False
    while (sol.n_it < max_iterations
           and not minimum_lmda
           and not persistent_bound_violation
           and not converged):

        sol.J = J(sol.x)  # evaluate Jacobian
        luJ = lu_factor(sol.J)  # storing the factorisation saves time later
        dx = lu_solve(luJ, -sol.F)  # compute ordinary Newton step
        dx_norm = np.linalg.norm(dx, ord=2)
        lmda_bounds = lambda_bounds(dx, sol.x)
        h = (lmda*np.linalg.norm((dxbar - dx), ord=2) * dx_norm
             / (np.linalg.norm(dxprev, ord=2) * np.linalg.norm(dxbar, ord=2)))
        lmda = update_lmda(sol.x, dx, h, lmda_bounds)

        # Create the (k+1)^0 values
        x_j = sol.x + lmda*dx

        # Check that all constraints are satisfied. If not, adjust lambda.
        # This must be done just before every call to F() *if* lambda has been increased:
        c_x_j = constraints(x_j)
        if not np.all(c_x_j < eps):  # x allowed to lie on constraints but not in forbidden area
            c_x = constraints(sol.x)
            violated_constraints = sorted([(i, c_x[i] / (c_x[i] - c_x_j[i])) for i in range(n_constraints) if c_x_j[i] >= eps], key=lambda x: x[1])
            lmda = lmda * violated_constraints[0][1]
            x_j = sol.x + lmda*dx

        # If the same current iterate is on a constraint,
        # and a very small lambda causes the next iterate to leave the
        # feasible region, then a new step direction must be found,
        # along with a new guess for lmda
        # We do this here using Lagrange multipliers
        if lmda < eps:
            active_constraint_indices = [i for i, vc in violated_constraints if vc < eps]
            inactive_constraint_indices = [i for i, vc in violated_constraints if vc >= eps]
            c_newton = constraints(sol.x + dx)[active_constraint_indices]
            c_A = linear_constraints[0][active_constraint_indices]
            x_n = sol.x + dx  # newton iterate
            if np.linalg.matrix_rank(c_A) == len(dx):  # if true, we must leave a constraint here
                n_act = len(active_constraint_indices)
                for i_rm in range(n_act):
                    potential_active_indices = [active_constraint_indices[i]
                                                for i in range(n_act) if i != i_rm]
                    c_newton = constraints(sol.x + dx)[potential_active_indices]
                    c_A = linear_constraints[0][potential_active_indices]
                    x_m = solve_constraint_lagrangian(x_n, sol.J, c_newton, c_A)[0]
                    if constraints(x_m)[active_constraint_indices[i_rm]] < 0.:
                        break
            else:
                x_m = solve_constraint_lagrangian(x_n, sol.J, c_newton, c_A)[0]

            dx = x_m - sol.x
            lmda_bounds = lambda_bounds(dx, sol.x)
            lmda = lmda_bounds[1]  # no a-priori maximum limit
            x_j = sol.x + lmda*dx

            # Check that the solution is still able to converge, i.e.
            # that the constraints aren't stopping our approach to a potential root
            x_j_min = sol.x + lmda_bounds[0]*dx  # because lmda must be getting smaller, no need to check constraints
            F_j_min = F(x_j_min)
            dxbar_j_min = lu_solve(luJ, -F_j_min)
            dxbar_j_min_norm = np.linalg.norm(dxbar_j_min, ord=2)

            # Newton step size must be decreasing and dx must be non-zero
            if dxbar_j_min_norm > dx_norm or np.linalg.norm(dx, ord=2) < eps:
                persistent_bound_violation = True

            # Now we need to check for newly violated constraints
            n_inactive = len(inactive_constraint_indices)
            c_x_j = constraints(x_j)[inactive_constraint_indices]
            if not np.all(c_x_j < eps):  # x allowed to lie on constraints but not in forbidden area
                c_x = constraints(sol.x)[inactive_constraint_indices]
                violated_constraints = sorted([(i, c_x[i] / (c_x[i] - c_x_j[i])) for i in range(n_inactive) if c_x_j[i] >= eps], key=lambda x: x[1])
                lmda = lmda * violated_constraints[0][1]
                x_j = sol.x + lmda*dx

        F_j = F(x_j)
        dxbar_j = lu_solve(luJ, -F_j)  # this is the simplified newton step
        dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

        if ((all(np.abs(dxbar_j) < tol)                  # <- Success requirements
             and all(np.abs(dx) < np.sqrt(10.*tol))      # <- avoids pathological cases
             and np.abs(lmda - lmda_bounds[1]) < eps)):  # <- end on a maximal newton step
            require_posteriori_loop = False              # <- No need for the a posteriori loop
            converged = True                             # <- Successful convergence
        else:
            require_posteriori_loop = True

        # Begin the a posteriori loop
        while (require_posteriori_loop and not minimum_lmda
               and not persistent_bound_violation):
            # Monotonicity check
            # always based on the Newton step, even if on a constraint
            if dxbar_j_norm <= dx_norm:
                if dxbar_j_norm < eps:  # <- occasionally the simplified newton step finds the exact solution
                    converged = True
                dxbar = dxbar_j
                sol.x = x_j
                sol.F = F_j

                require_posteriori_loop = False  # return to Newton step
                sol.n_it += 1  # move to next iteration
                dxprev = dx  # to calculate the next value of h
            else:
                if np.abs(lmda - lmda_bounds[0]) < eps:
                    minimum_lmda = True
                h_j = (2./lmda)*np.linalg.norm((dxbar_j - (1. - lmda)*dx), ord=2)/dx_norm
                lmda_j = min(lmda_bounds[1], 1./h_j)
                lmda = min(lmda_j, lmda/2.)
                lmda = max(lmda, lmda_bounds[0])  # allows a check of monotonicity once at minimum lmda

                x_j = sol.x + lmda*dx  # because lmda must be getting smaller, no need to check constraints
                F_j = F(x_j)
                dxbar_j = lu_solve(luJ, -F_j)
                dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

        if store_iterates:
            sol.iterates.x.append(sol.x)
            sol.iterates.F.append(sol.F)
            sol.iterates.lmda.append(lmda)

    if converged and not persistent_bound_violation:
        sol.x = x_j + dxbar_j
        # Even if the solver succeeds, there may be a small chance that the last simplified Newton step
        # shifts the solution just outside the constraints.
        # If so, shift the solution back to the allowed region
        c_x = constraints(sol.x)
        if not np.all(c_x <= 0.):  # x allowed to lie on constraints but not in forbidden area
            sol.x -= dxbar_j

    sol.F = F(sol.x)
    sol.F_norm = np.linalg.norm(sol.F, ord=2)
    sol.J = J(sol.x)

    if store_iterates:
        sol.iterates.x = np.array(sol.iterates.x)
        sol.iterates.F = np.array(sol.iterates.F)

    sol.success = False
    if converged:
        sol.success = True
        sol.code = 0
        sol.text = 'The solver successfully found a root after {0} iterations'.format(sol.n_it)
    elif minimum_lmda:
        sol.code = 1
        sol.text = 'The function is too non-linear for lower lambda bound ({0})'.format(lmda_bounds[0])
    elif persistent_bound_violation:
        sol.code = 2
        sol.text = 'The descent vector crosses the constraints with the following indices: {0}'.format([i for i, lmda in violated_constraints])
    elif sol.n_it == max_iterations:
        sol.code = 3
        sol.text = 'The solver reached max_iterations ({0})'.format(max_iterations)
    else:
        raise Exception('Unknown termination of solver')
    return sol
