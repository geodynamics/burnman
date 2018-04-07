import numpy as np
from scipy.linalg import lu_factor, lu_solve
from collections import namedtuple

def damped_newton_solve(F, J, guess, tol=1.e-6,
                        max_iterations=100,
                        lambda_bounds=lambda dx, x: (1.e-8, 1.),
                        constraints=lambda x: np.array([-1.])):
    """
    Solver for the multivariate nonlinear system F(x)=0 
    with Jacobian J(x), using the damped affine invariant modification
    to Newton's method (Deuflhard, 1974;1975;2004).
    Here we follow the algorithm as described in Nowak and Weimann (1991):
    [Technical Report TR-91-10, Algorithm B]

    Inequality constraints are provided by the function constraints(x),
    which returns a 1D numpy array. The constraints are satisfied if all the
    returned values are <=0. If any constraints are not satisfied by the current
    value of lambda, lambda is reduced to satisfy all the constraints.

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
    constraints : function of x
        Returns the LHS of the inequality constraints(x)
        as a 1D numpy array. The constraints are satisfied if
        all the elements of the array are less than or equal to
        zero.

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
        
    This function is available as ``burnman.damped_newton_solve``.    
    """

    # Make sure damping factor is within bounds, and that the bounds are reasonable
    # Problem classes in Nowak and Weimann (1991); [lmda_min, lmda_max]:
    # linear: [0.1, 1.]
    # mildly nonlinear: [1.e-4, 1.]
    # highly nonlinear: [1.e-2, 1.e-4]
    # extremely nonlinear: [1.e-4, 1.e-8]
    eps = 2.*np.finfo(float).eps
    assert np.all(constraints(guess) < eps), 'The starting guess is outside the supplied constraints.'

    if not isinstance(tol, float):
        assert len(tol) < len(guess), 'tol must either be a float or an array like guess.'

    sol = namedtuple('Solution', ['x', 'n_it', 'F', 'F_norm', 'J', 'code', 'text', 'success'])

    # evaluate system
    sol.x = guess
    sol.F = F(sol.x) 

    # Begin Newton loop
    sol.n_it = 0
    n_constraints = len(constraints(sol.x))
    minimum_lmda = False
    converged = False
    bound_violation = False
    persistent_bound_violation = False
    while (sol.n_it < max_iterations and
           not minimum_lmda and
           not persistent_bound_violation and
           not converged):

        sol.J = J(sol.x) # evaluate Jacobian
        luJ = lu_factor(sol.J) # storing the factorisation saves time later
        dx = lu_solve(luJ, -sol.F) # compute ordinary Newton step
        dx_norm = np.linalg.norm(dx, ord=2)

        lmda_bounds = lambda_bounds(dx, sol.x)
        assert lmda_bounds[1] < 1. + eps, 'The highest upper bound for lambda is 1. (a full Newton step)'
        assert lmda_bounds[0] > 1.e-8 - eps, 'The lowest lower bound for lambda is 1.e-8 (suitable only for extremely nonlinear systems)'
        
        # Calculate a priori damping factor
        if sol.n_it > 0:
            h = (lmda * np.linalg.norm((dxbar - dx), ord=2) * np.linalg.norm(dx, ord=2) /
                 (np.linalg.norm(dxprev, ord=2) * np.linalg.norm(dxbar, ord=2)))
            lmda_j = min(1./(h+eps), lmda_bounds[1]) # this is lmda_j^0
        else:
            lmda_j = lmda_bounds[1] # this is lmda_0_0

        lmda = max(lmda_j, lmda_bounds[0])

        # Create the (k+1)^0 values 
        x_j = sol.x + lmda*dx        

        # Check that all constraints are satisfied.
        # If not, adjust lambda. This must be done just before every call to F() *if* lambda has been increased:
        c_x_j = constraints(x_j)
        if not np.all(c_x_j < eps): # x allowed to lie on constraints but not in forbidden area
            c_x = constraints(sol.x)
            lmda = lmda * min([c_x[i] / (c_x[i] - c_x_j[i]) for i in range(n_constraints) if c_x_j[i]>=eps])
            x_j = sol.x + lmda*dx
            if bound_violation:
                persistent_bound_violation=True
            bound_violation=True
        else:
            bound_violation=False # reset if a violation does not recur

        F_j = F(x_j)

        dxbar_j = lu_solve(luJ, -F_j)
        dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)

        if (((all(np.abs(dxbar_j) < tol) and                     # <- Success requirements
              all(np.abs(dx) < np.sqrt(10.*tol))) or             # <- avoids pathological cases
             dxbar_j_norm < eps) and                             # <- occasionally the simplified newton step finds the exact solution
            np.abs(lmda - lmda_bounds[1]) < eps) :               # <- end on a maximal newton step
            require_posteriori_loop = False                      # <- No need for the a posteriori loop
            converged = True                                     # <- Successful convergence
        else:
            require_posteriori_loop = True

        # Begin the a posteriori loop
        while require_posteriori_loop and not minimum_lmda:
            # Monotonicity check
            if dxbar_j_norm <= dx_norm: 
                dxbar = dxbar_j
                sol.x = x_j
                sol.F = F_j

                require_posteriori_loop = False # return to Newton step
                sol.n_it += 1 # move to next iteration
                dxprev = dx # to calculate the next value of h
            else:
                if np.abs(lmda - lmda_bounds[0]) < eps:
                    minimum_lmda = True
                h_j = (2./lmda)*np.linalg.norm((dxbar_j - (1. - lmda)*dx), ord=2)/dx_norm
                lmda_j = min(lmda_bounds[1], 1./h_j)
                lmda = min(lmda_j, lmda/2.)
                lmda = max(lmda, lmda_bounds[0]) # allows a check of monotonicity once at minimum lmda

                x_j = sol.x + lmda*dx # because lmda must be getting smaller, no need to check constraints
                F_j = F(x_j)
                dxbar_j = lu_solve(luJ, -F_j)
                dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)


    if not persistent_bound_violation:
        sol.x = x_j + dxbar_j
        # Even if the solver succeeds, there may be a small chance that the solution lies
        # just outside the constraints. If so, print a warning and shift the solution back
        # to the allowed region
        c_x = constraints(sol.x)
        if not np.all(c_x < eps): # x allowed to lie on constraints but not in forbidden area
            sol.x -= dxbar_j
            print('Warning: The solution appears to lie just outside the chosen constraints.')

    sol.F = F(sol.x)
    sol.F_norm = np.linalg.norm(sol.F, ord=2)
    sol.J = J(sol.x)

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
        sol.text = 'The descent vector crosses one or more constraints'
    elif sol.n_it == max_iterations:
        sol.code = 3
        sol.text = 'The solver reached max_iterations ({0})'.format(max_iterations)
    else:
        raise Exception('Unknown termination of solver')
    return sol
