import unittest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from util import BurnManTest
from burnman.nonlinear_solvers import damped_newton_solve
import numpy as np

class test_solvers(BurnManTest):
    # The following tests are taken from More et al, 1981; doi:10.1145/355934.355936
    # (Testing Unconstrained Optimization Software)
    def test_dns_rosenbrock_1(self):
        F = lambda x: np.array([10.*(x[1] - x[0]*x[0]), 1. - x[0]])
        J = lambda x: np.array([[-20.*x[0], 10.],
                                [-1., 0.]])

        sol = damped_newton_solve(F, J, guess=np.array([-1.2, 1.]))
        self.assertArraysAlmostEqual(sol.x, [1., 1.])
        assert(sol.success)

    def test_dns_rosenbrock_1_w_constraint_violation(self):
        F = lambda x: np.array([10.*(x[1] - x[0]*x[0]), 1. - x[0]])
        J = lambda x: np.array([[-20.*x[0], 10.],
                                [-1., 0.]])
        C = lambda x: np.array([x[0]])

        sol = damped_newton_solve(F, J, guess=np.array([-1.2, 1.]), constraints=C)

        self.assertFloatEqual(sol.x[0], 0.)
        assert(sol.code == 2)

    def test_dns_rosenbrock_1_w_temporary_constraint_violation(self):
        F = lambda x: np.array([10.*(x[1] - x[0]*x[0]), 1. - x[0]])
        J = lambda x: np.array([[-20.*x[0], 10.],
                                [-1., 0.]])
        C = lambda x: np.array([x[0] - 1., -3. - x[1]])

        # Here the solver takes two iterations, the first of which
        # has a full Newton step which violates the constraints
        # [1., -3.84]
        sol = damped_newton_solve(F, J, guess=np.array([-1.2, 1.]), constraints=C)

        # The solution lies on one of the constraints
        self.assertArraysAlmostEqual(sol.x, [1., 1.])

    def test_dns_freudenstein_roth_2(self):
        F = lambda x: np.array([-13. + x[0] + ((5. - x[1])*x[1] - 2.)*x[1],
                                -29. + x[0] + ((x[1] + 1.)*x[1] - 14.)*x[1]])
        J = lambda x: np.array([[1., (10. - 3.*x[1])*x[1] - 2.],
                                [1., (3.*x[1] + 2.)*x[1] - 14.]])
        
        # sol = damped_newton_solve(F, J, guess=np.array([0.5, -2.])) # does not converge
        sol = damped_newton_solve(F, J, guess=np.array([10., 10.]))
        self.assertArraysAlmostEqual(sol.x, [5., 4.])
        assert(sol.success)

    def test_dns_powell_badly_scaled_3(self):
        F = lambda x: np.array([1.e4*x[0]*x[1] - 1., np.exp(-x[0]) + np.exp(-x[1]) - 1.0001])
        J = lambda x: np.array([[1.e4*x[1], 1.e4*x[0]],
                                [-np.exp(-x[0]), -np.exp(-x[1])]])
        
        sol = damped_newton_solve(F, J, guess=np.array([0., 1.]))
        self.assertArraysAlmostEqual(sol.x, [1.09815933e-05, 9.10614674])
        assert(sol.success)

    def test_dns_powell_singular_13(self):
        F = lambda x: np.array([x[0] + 10.*x[1], np.sqrt(5.)*(x[2] - x[3]), np.power(x[1] - 2.*x[2], 2.), np.sqrt(10.)*np.power(x[0] - x[3], 2.)])
        J = lambda x: np.array([[1., 10., 0., 0.],
                                [0., 0., np.sqrt(5.), -np.sqrt(5.)],
                                [0., 2.*x[1] - 4.*x[2], 8.*x[2] - 4.*x[1], 0.],
                                [np.sqrt(10.)*(2.*x[0] - 2.*x[3]), 0., 0., np.sqrt(10.)*(2.*x[3] - 2.*x[0])]])
        
        sol = damped_newton_solve(F, J, guess=np.array([3., -1., 0., 1.]))
        self.assertArraysAlmostEqual(sol.x + 1., [1., 1., 1., 1.])
        assert(sol.success)


    def test_dns_broyden_tridiagonal_30(self):
        def F(x):
            xpad = np.concatenate(([0.], x, [0.]))
            f = np.zeros((len(xpad)))
            for i in range(1, len(xpad)-1):
                f[i] = (3. - 2.*xpad[i])*xpad[i] - xpad[i-1] - 2.*xpad[i+1] + 1.
            return f[1:-1]

        def J(x):
            xpad = np.concatenate(([0.], x, [0.]))
            j = np.zeros((len(xpad), len(xpad)))
            for i in range(1, len(xpad)-1):
                j[i,i-1] = -1.
                j[i,i] = 3. - 4.*xpad[i]
                j[i,i+1] = -2.
                
            return j[1:-1, 1:-1]

        # check convergence for several test cases
        expected_solutions = [[-0.28077641],
                              [-0.45328926, -0.38540505],
                              [-0.52677285, -0.56764891, -0.41031222],
                              [-0.55457673, -0.63942044, -0.59070079, -0.41526838],
                              [-0.5648284, -0.66627372, -0.66091704, -0.59505005, -0.41620111]]
        
        for n in range(1, 6):
            guess = -1.*np.ones((n))
            sol = damped_newton_solve(F, J, guess=guess)
            self.assertArraysAlmostEqual(sol.x, expected_solutions[n-1])
            assert(sol.success)

if __name__ == '__main__':
    unittest.main()
