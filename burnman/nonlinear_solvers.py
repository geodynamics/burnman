import numpy as np
from scipy.linalg import lu_factor, lu_solve
from collections import namedtuple

def damped_newton_solve(F, J, x, tol=1.e-6,
                        max_iterations=100,
                        lmda_bounds=[1.e-8, 1.],
                        constraints=lambda x: -1.):
    """
    Solver for the multivariate nonlinear system F(x)=0 
    with Jacobian J(x) and constraints C(x), using the 
    damped affine invariant modification to Newton's method 
    (Deuflhard, 1974;1975;2004), 
    Here we follow the algorithm as described in Nowak and Weimann (1991):
    [Technical Report TR-91-10, Algorithm B]

    The iteration continues until three conditions are satisfied:
    ||dx (simplified newton step)|| < tol
    ||dx (full Newton step)|| < sqrt(10*tol) [avoiding pathological behaviour] and
    lambda = lambda_max (lambda = 1 for a full Newton step).

    Inequality constraints are provided by the function C(x), which returns a 
    1D numpy array. The constraints are satisfied if the values are <=0.
    If any constraints are not satisfied by the current value of lambda, 
    lambda is reduced to satisfy all the constraints.
    """

    # Make sure damping factor is within bounds, and that the bounds are reasonable
    # Problem classes in Nowak and Weimann (1991); lmda_max, lmda_min:
    # linear: 1, 1
    # mildly nonlinear: 1, 1e-4
    # highly nonlinear: 1e-2, 1e-4
    # extremely nonlinear: 1e-4, 1e-8
    eps = 2.*np.finfo(float).eps
    assert lmda_bounds[1] < 1. + eps, 'The highest upper bound for lmda is 1. (a full Newton step)'
    assert lmda_bounds[0] > 1.e-8 - eps, 'The lowest lower bound for lmda is 1.e-8 (suitable only for extremely nonlinear systems)'
    assert np.all(constraints(x) < eps), 'The starting guess is outside the supplied constraints.'
           
    lmda_j = lmda_bounds[1] # this is lmda_0_0

    sol = namedtuple('Solution', ['x', 'n_it', 'F', 'F_norm', 'J', 'code', 'text'])

    # evaluate system
    sol.F = F(x) 
            
    # Begin Newton loop
    sol.n_it = 0
    minimum_lmda = False
    converged = False
    persistent_bound_violation = False
    while (sol.n_it < max_iterations and
           not minimum_lmda and
           not persistent_bound_violation and
           not converged):
           
        sol.J = J(x) # evaluate Jacobian
        luJ = lu_factor(sol.J) # storing the factorisation saves time later
        dx = lu_solve(luJ, -sol.F) # compute ordinary Newton correction

        # Calculate a priori damping factor
        if sol.n_it > 0:
            
            h = (lmda * np.linalg.norm((dxbar - dx), ord=2) * np.linalg.norm(dx, ord=2) /
                 (np.linalg.norm(dxprev, ord=2) * np.linalg.norm(dxbar, ord=2)))
            lmda_j = min(1./h, lmda_bounds[1]) # this is lmda_j^0


        lmda = max(lmda_j, lmda_bounds[0])

        # Create the (k+1)^0 values 
        x_j = x + lmda*dx        

        # Check that all constraints are satisfied.
        # If not, adjust lambda. This must be done just before every call to F() *if* lambda has been increased:
        c_x_j = constraints(x_j)
        if not np.all(c_x_j < eps): # x allowed to lie on constraints but not in forbidden area
            c_x = constraints(x)
            if not np.all(c_x < -eps): # terminate if x tries to slide off constraint into disallowed area
                persistent_bound_violation = True
                x_j = x
            else:
                lmda = lmda * min(c_x / (c_x - c_x_j))
                x_j = x + lmda*dx  
                
        F_j = F(x_j)
        dxbar_j = lu_solve(luJ, -F_j)
                
            
        dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)
        dx_norm = np.linalg.norm(dx, ord=2)
        if (dxbar_j_norm < tol and                             # <- Success requirements
            dx_norm < np.sqrt(10.*tol) and                     # <- avoids pathological cases
            np.abs(lmda - lmda_bounds[1]) < eps):              # <- end on a full newton step
            require_posteriori_loop = False                    # <- No need for the a posteriori loop
            converged = True                                   # <- Successful convergence
        else:
            require_posteriori_loop = True
            

        # Begin the a posteriori loop
        while require_posteriori_loop and not minimum_lmda:
            # Monotonicity check
            if dxbar_j_norm <= dx_norm: 
                dxbar = dxbar_j
                x = x_j
                sol.F = F_j
                sol.n_it += 1
                require_posteriori_loop = False # return to Newton step
            else:
                if np.abs(lmda - lmda_bounds[0]) < eps:
                    minimum_lmda = True
                h_j = (2./lmda)*np.linalg.norm((dxbar_j - (1. - lmda)*dx), ord=2)/dx_norm
                lmda_j = min(lmda_bounds[1], 1./h_j)
                lmda = min(lmda_j, lmda/2.)
                lmda = max(lmda, lmda_bounds[0]) # allows a check of monotonicity once at minimum lmda
            
                x_j = x + lmda*dx # because lmda must be getting smaller, no need to check constraints
                F_j = F(x_j)
                dxbar_j = lu_solve(luJ, -F_j)
                dxbar_j_norm = np.linalg.norm(dxbar_j, ord=2)
            
        dxprev = dx
        sol.F_norm = np.linalg.norm(sol.F, ord=2)

    if not persistent_bound_violation:
        sol.x = x_j + dxbar_j
        # Even if the solver succeeds, there may be a small chance that the solution lies
        # just outside the constraints. If so, print a warning and shift the solution back
        # to the allowed region
        c_x = constraints(x)
        if not np.all(c_x < eps): # x allowed to lie on constraints but not in forbidden area
            sol.x -= dxbar_j
            print('Warning: The solution appears to lie just outside the chosen constraints.')
    else:
        sol.x = x
            
    sol.F = F(sol.x)
    sol.res = sol.F
    
    sol.success = False
    if converged:
        sol.success = True
        sol.code = 0 
        sol.text = 'solution converged after {0} iterations'.format(sol.n_it)
    elif minimum_lmda:
        sol.code = 1 
        sol.text = 'solution is too non-linear for lower lambda bound ({0})'.format(lmda_bounds[0])
    elif persistent_bound_violation:
        sol.code = 2 
        sol.text = 'if solution exists, it is outside the region allowed by the constraints'
    elif sol.n_it == max_iterations:
        sol.code = 3
        sol.text = 'solution reached max_iterations ({0})'.format(max_iterations)
    else:
        raise Exception('Unknown termination of solver')
    return sol
