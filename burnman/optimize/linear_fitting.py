# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import importlib
import numpy as np
from scipy.linalg import inv, sqrtm
import warnings

try:
    cp = importlib.import_module("cvxpy")
except ImportError as err:
    print(
        f"Warning: {err}. " "For full functionality of BurnMan, please install cvxpy."
    )


def weighted_constrained_least_squares(
    A, b, Cov_b=None, equality_constraints=None, inequality_constraints=None
):
    """
    Solves a weighted, constrained least squares problem using cvxpy.
    The objective function is to minimize the following:
    sum_squares(Cov_b^(-1/2).A.x - Cov_b^(-1/2).b))
    subject to
    C.x == c
    D.x <= d

    :param A: An array defining matrix A in the objective function above.
    :type A: 2D numpy array

    :param b: An array defining vector b in the objective function above.
    :type b: numpy array

    :param Cov_b: A covariance matrix associated with b
    :type Cov_b: 2D numpy array

    :param equality_constraints: A list containing the matrices C and c
        in the objective function above.
    :type equality_constraints: list containing a 2D array and 1D array

    :param inequality_constraints: A list containing the matrices D and d
        in the objective function above.
    :type inequality_constraints: list containing a 2D array and 1D array

    :returns: Tuple containing the optimized phase amounts (1D numpy.array),
        a covariance matrix corresponding to the optimized phase amounts
        (2D numpy.array), and the weighted residual of the fitting procedure
        (a float).
    :rtype: tuple
    """

    if Cov_b is None:
        Cov_b = np.eye(len(b))

    # Create the standard weighted least squares objective function
    # (https://stats.stackexchange.com/a/333551)
    n_vars = A.shape[1]
    m = inv(sqrtm(Cov_b))
    mA = m @ A
    mb = m @ b
    x = cp.Variable(n_vars)
    objective = cp.Minimize(cp.sum_squares(mA @ x - mb))

    constraints = []
    if equality_constraints is not None:
        n_eq_csts = len(equality_constraints[0])
        constraints = [
            equality_constraints[0][i] @ x == equality_constraints[1][i]
            for i in range(n_eq_csts)
        ]

    if inequality_constraints is not None:
        n_ineq_csts = len(inequality_constraints[0])
        constraints.extend(
            [
                inequality_constraints[0][i] @ x <= inequality_constraints[1][i]
                for i in range(n_ineq_csts)
            ]
        )

    # Set up the problem and solve it
    warns = []
    if len(constraints) > 1:
        prob = cp.Problem(objective, constraints)
    else:
        prob = cp.Problem(objective)

    try:
        with warnings.catch_warnings(record=True) as w:
            res = prob.solve()
            popt = np.array([x.value[i] for i in range(len(A.T))])
            warns.extend(w)
    except Exception:
        try:
            with warnings.catch_warnings(record=True) as w:
                res = prob.solve(solver=cp.ECOS)
                popt = np.array([x.value[i] for i in range(len(A.T))])
                warns.extend(w)
        except Exception as e:
            raise Exception(e)

    # Calculate the covariance matrix
    # (also from https://stats.stackexchange.com/a/333551)
    inv_Cov_b = np.linalg.inv(Cov_b)
    pcov = np.linalg.inv(A.T.dot(inv_Cov_b.dot(A)))

    return (popt, pcov, res)
