import unittest
from util import BurnManTest
import numpy as np

from burnman.optimize.nonlinear_solvers import damped_newton_solve


class test_solvers(BurnManTest):
    # The following tests are taken from More et al, 1981; doi:10.1145/355934.355936
    # (Testing Unconstrained Optimization Software)
    def test_dns_rosenbrock_1(self):
        F = lambda x: np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]])
        J = lambda x: np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])

        sol = damped_newton_solve(F, J, guess=np.array([-1.2, 1.0]))
        self.assertArraysAlmostEqual(sol.x, [1.0, 1.0])
        assert sol.success

    def test_dns_rosenbrock_1_w_constraint_violation(self):
        F = lambda x: np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]])
        J = lambda x: np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])
        C = (np.array([[1.0, 0.0]]), np.array([0.0]))

        sol = damped_newton_solve(
            F, J, guess=np.array([-1.2, 1.0]), linear_constraints=C
        )
        self.assertFloatEqual(sol.x[0], 0.0)
        assert sol.code == 2

    def test_dns_rosenbrock_1_w_temporary_constraint_violation(self):
        F = lambda x: np.array([10.0 * (x[1] - x[0] * x[0]), 1.0 - x[0]])
        J = lambda x: np.array([[-20.0 * x[0], 10.0], [-1.0, 0.0]])
        C = (np.array([[1.0, 0.0], [0.0, -1.0]]), np.array([-1.0, -3.0]))
        # Here the solver takes two iterations, the first of which
        # has a full Newton step which violates one of the constraints
        # [1., -3.84]
        sol = damped_newton_solve(
            F, J, guess=np.array([-1.2, 1.0]), linear_constraints=C
        )
        # The solution lies on one of the constraints
        self.assertArraysAlmostEqual(sol.x, [1.0, 1.0])

    def test_dns_freudenstein_roth_2(self):
        F = lambda x: np.array(
            [
                -13.0 + x[0] + ((5.0 - x[1]) * x[1] - 2.0) * x[1],
                -29.0 + x[0] + ((x[1] + 1.0) * x[1] - 14.0) * x[1],
            ]
        )
        J = lambda x: np.array(
            [
                [1.0, (10.0 - 3.0 * x[1]) * x[1] - 2.0],
                [1.0, (3.0 * x[1] + 2.0) * x[1] - 14.0],
            ]
        )

        # sol = damped_newton_solve(F, J, guess=np.array([0.5, -2.])) # does not converge
        sol = damped_newton_solve(F, J, guess=np.array([10.0, 10.0]))
        self.assertArraysAlmostEqual(sol.x, [5.0, 4.0])
        assert sol.success

    def test_dns_powell_badly_scaled_3(self):
        F = lambda x: np.array(
            [1.0e4 * x[0] * x[1] - 1.0, np.exp(-x[0]) + np.exp(-x[1]) - 1.0001]
        )
        J = lambda x: np.array(
            [[1.0e4 * x[1], 1.0e4 * x[0]], [-np.exp(-x[0]), -np.exp(-x[1])]]
        )

        sol = damped_newton_solve(F, J, guess=np.array([0.0, 1.0]))
        self.assertArraysAlmostEqual(sol.x, [1.09815933e-05, 9.10614674])
        assert sol.success

    def test_dns_powell_singular_13(self):
        F = lambda x: np.array(
            [
                x[0] + 10.0 * x[1],
                np.sqrt(5.0) * (x[2] - x[3]),
                np.power(x[1] - 2.0 * x[2], 2.0),
                np.sqrt(10.0) * np.power(x[0] - x[3], 2.0),
            ]
        )
        J = lambda x: np.array(
            [
                [1.0, 10.0, 0.0, 0.0],
                [0.0, 0.0, np.sqrt(5.0), -np.sqrt(5.0)],
                [0.0, 2.0 * x[1] - 4.0 * x[2], 8.0 * x[2] - 4.0 * x[1], 0.0],
                [
                    np.sqrt(10.0) * (2.0 * x[0] - 2.0 * x[3]),
                    0.0,
                    0.0,
                    np.sqrt(10.0) * (2.0 * x[3] - 2.0 * x[0]),
                ],
            ]
        )

        sol = damped_newton_solve(F, J, guess=np.array([3.0, -1.0, 0.0, 1.0]))
        self.assertArraysAlmostEqual(sol.x + 1.0, [1.0, 1.0, 1.0, 1.0])
        assert sol.success

    def test_dns_broyden_tridiagonal_30(self):
        def F(x):
            xpad = np.concatenate(([0.0], x, [0.0]))
            f = np.zeros((len(xpad)))
            for i in range(1, len(xpad) - 1):
                f[i] = (
                    (3.0 - 2.0 * xpad[i]) * xpad[i]
                    - xpad[i - 1]
                    - 2.0 * xpad[i + 1]
                    + 1.0
                )
            return f[1:-1]

        def J(x):
            xpad = np.concatenate(([0.0], x, [0.0]))
            j = np.zeros((len(xpad), len(xpad)))
            for i in range(1, len(xpad) - 1):
                j[i, i - 1] = -1.0
                j[i, i] = 3.0 - 4.0 * xpad[i]
                j[i, i + 1] = -2.0

            return j[1:-1, 1:-1]

        # check convergence for several test cases
        expected_solutions = [
            [-0.28077641],
            [-0.45328926, -0.38540505],
            [-0.52677285, -0.56764891, -0.41031222],
            [-0.55457673, -0.63942044, -0.59070079, -0.41526838],
            [-0.5648284, -0.66627372, -0.66091704, -0.59505005, -0.41620111],
        ]

        for n in range(1, 6):
            guess = -1.0 * np.ones((n))
            sol = damped_newton_solve(F, J, guess=guess)
            self.assertArraysAlmostEqual(sol.x, expected_solutions[n - 1])
            assert sol.success

    # The following tests use the generalised Rosenbrock function:
    # f(x, y) = (a - x)^2 + b(y - x^2)^2
    # solving for f'=0
    def test_dns_rosenbrock_generalised(self):
        a = 1.0
        b = 15.0  # takes >100 iterations for b>16, ~570 iterations for b=100
        # this is a tricky root-finding problem when x0 < 0, because the gradient in the descent valley is very small compared to the gradient at the valley edges. Only very small fractions of a full Newton step can be taken if global convergence is to be guaranteed.

        F = lambda x: np.array(
            [
                2.0 * (x[0] - a) + 4.0 * b * x[0] * (x[0] * x[0] - x[1]),
                2.0 * b * (x[1] - x[0] * x[0]),
            ]
        )

        J = lambda x: np.array(
            [
                [12.0 * b * x[0] * x[0] - 4.0 * b * x[1] + 2.0, -4.0 * b * x[0]],
                [-4.0 * b * x[0], 2.0 * b],
            ]
        )
        guess = np.array([-1, 1.0])
        sol = damped_newton_solve(
            F, J, guess=guess, max_iterations=600, tol=1.0e-15, store_iterates=True
        )

        self.assertArraysAlmostEqual(sol.x, [a, a * a])
        assert sol.success

    def test_dns_rosenbrock_generalised_w_descent_path_constraint(self):
        # Sometimes there are inequality constraints on feasible regions.
        # The damped newton solver uses Lagrangian multipliers to
        # maximize the approach to F=0 (in the L1-norm sense)
        # whilst satisfying all the constraints.
        a = 1.0
        b = 15.0
        F = lambda x: np.array(
            [
                2.0 * (x[0] - a) + 4.0 * b * x[0] * (x[0] * x[0] - x[1]),
                2.0 * b * (x[1] - x[0] * x[0]),
            ]
        )

        J = lambda x: np.array(
            [
                [12.0 * b * x[0] * x[0] - 4.0 * b * x[1] + 2.0, -4.0 * b * x[0]],
                [-4.0 * b * x[0], 2.0 * b],
            ]
        )

        C = (np.array([[0.0, 1.0]]), np.array([-1.1]))

        guess = np.array([-4, -3.0])
        sol = damped_newton_solve(
            F, J, guess=guess, linear_constraints=C, max_iterations=100, tol=1.0e-15
        )

        self.assertArraysAlmostEqual(sol.x, [a, a * a])
        assert sol.success

    def test_dns_rosenbrock_generalised_w_two_descent_path_constraints(self):
        # This test is like the one above, except that there are
        # two successive constraints along which the solver must move
        a = 1.0
        b = 15.0
        F = lambda x: np.array(
            [
                2.0 * (x[0] - a) + 4.0 * b * x[0] * (x[0] * x[0] - x[1]),
                2.0 * b * (x[1] - x[0] * x[0]),
            ]
        )

        J = lambda x: np.array(
            [
                [12.0 * b * x[0] * x[0] - 4.0 * b * x[1] + 2.0, -4.0 * b * x[0]],
                [-4.0 * b * x[0], 2.0 * b],
            ]
        )

        C = (np.array([[0.0, 1.0], [-1.0, 1.0]]), np.array([-1.1, -3.5]))

        guess = np.array([-4, -3.0])
        sol = damped_newton_solve(
            F, J, guess=guess, linear_constraints=C, max_iterations=200, tol=1.0e-15
        )

        self.assertArraysAlmostEqual(sol.x, [a, a * a])
        assert sol.success

    def test_dns_rosenbrock_generalised_w_constraint_violation(self):
        # This test is like the one above, except that there are
        # two successive constraints along which the solver must move
        a = 1.0
        b = 15.0
        F = lambda x: np.array(
            [
                2.0 * (x[0] - a) + 4.0 * b * x[0] * (x[0] * x[0] - x[1]),
                2.0 * b * (x[1] - x[0] * x[0]),
            ]
        )

        J = lambda x: np.array(
            [
                [12.0 * b * x[0] * x[0] - 4.0 * b * x[1] + 2.0, -4.0 * b * x[0]],
                [-4.0 * b * x[0], 2.0 * b],
            ]
        )

        y_max = 0.6
        C = (np.array([[0.0, 1.0]]), np.array([-y_max]))

        guess = np.array([-4, -3.0])
        sol = damped_newton_solve(
            F,
            J,
            guess=guess,
            linear_constraints=C,
            max_iterations=200,
            tol=1.0e-15,
            store_iterates=True,
        )

        self.assertFloatEqual(sol.x[1], y_max)
        assert sol.code == 2


if __name__ == "__main__":
    unittest.main()
