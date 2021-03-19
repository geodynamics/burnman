# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
import cvxpy as cp
from scipy.linalg import inv, sqrtm
from collections import Counter
import warnings
from .tools import merge_two_dicts


def weighted_constrained_least_squares(A, b, Cov_b,
                                       equality_constraints=None,
                                       inequality_constraints=None):
    """

    """

    # Create the standard weighted least squares objective function
    # (https://stats.stackexchange.com/a/333551)
    n_vars = A.shape[1]
    m = inv(sqrtm(Cov_b))
    mA = m@A
    mb = m@b
    x = cp.Variable(n_vars)
    objective = cp.Minimize(cp.sum_squares(mA@x - mb))

    constraints = []
    if equality_constraints is not None:
        n_eq_csts = len(equality_constraints[0])
        constraints = [equality_constraints[0][i]@x
                       == equality_constraints[1][i]
                       for i in range(n_eq_csts)]

    if inequality_constraints is not None:
        n_ineq_csts = len(inequality_constraints[0])
        constraints.extend([inequality_constraints[0][i]@x
                            <= inequality_constraints[1][i]
                            for i in range(n_ineq_csts)])

    # Set up the problem and solve it
    warns = []
    prob = cp.Problem(objective, constraints)
    try:
        with warnings.catch_warnings(record=True) as w:
            res = prob.solve(solver=cp.ECOS)
            popt = np.array([x.value[i] for i in range(len(A.T))])
            warns.extend(w)
    except Exception:
        print('ECOS Solver failed. Trying default solver.')
        try:
            with warnings.catch_warnings(record=True) as w:
                res = prob.solve()
                popt = np.array([x.value[i] for i in range(len(A.T))])
                warns.extend(w)
        except Exception as e:
            raise Exception(e)

    # Calculate the covariance matrix
    # (also from https://stats.stackexchange.com/a/333551)
    inv_Cov_b = np.linalg.inv(Cov_b)
    pcov = np.linalg.inv(A.T.dot(inv_Cov_b.dot(A)))

    return (popt, pcov, res)


def fit_composition_to_solution(fitted_variables,
                                variable_values, variable_covariances,
                                solution):
    """
    It is assumed that any elements not in composition were not measured (but may exist in unknown quantities).
    If distinct oxidation states or site occupancies were measured
    (by Moessbauer, for example), then formulae should be modified

    fitted elements should be a list of strings
    composition and compositional uncertainties should be given as arrays
    formulae should either be given as arrays (with columns corresponding to elements) or as dictionaries
    If the compositional uncertainties can either be sigmas or covariance matrices

    The composition and associated uncertainties in endmember *amounts*, not proportions. If normalize=True, then the endmember amounts are normalized to a total of one.
    """

    if type(formulae[0]) is dict or type(formulae[0]) is Counter:
        stoichiometric_matrix = np.array([[f[e] if e in f else 0. for e in fitted_elements] for f in formulae])
    else:
        stoichiometric_matrix = formulae

    b = composition

    if len(compositional_uncertainties.shape) == 1:
        b_uncertainties = np.diag(compositional_uncertainties *
                                  compositional_uncertainties)
    else:
        b_uncertainties = compositional_uncertainties

    if np.linalg.det(b_uncertainties) < 1.e-30: # ensure uncertainty matrix is not singular
        #warnings.warn('The compositional covariance matrix for the {0} solution is nearly singular or not positive-definite (determinant = {1}). '
        #              'This is likely to be because your fitting parameters are not independent. '
        #              'For now, we increase all diagonal components by 1%. '
        #              'However, you may wish to redefine your problem.'.format(name, np.linalg.det(b_uncertainties)))

        b_uncertainties += np.diag(np.diag(b_uncertainties))*0.01

    endmember_constraints = lambda site_occ: [{'type': 'ineq', 'fun': lambda x, eq=eq: eq.dot(x)}
                                              for eq in site_occ]
    cons = endmember_constraints(endmember_site_occupancies.T)    

    fn = lambda A, *proportions: A.dot(proportions)
    popt, pcov = curve_fit(fn, A, b,
                           p0=np.array([0. for i in range(len(A.T))]),
                           sigma=b_uncertainties, absolute_sigma=True)

    res = np.sqrt((A.dot(popt) - b).dot(np.linalg.solve(b_uncertainties, A.dot(popt) - b)))

    # Check constraints
    if any([c['fun'](popt)<0. for c in cons]):
        warnings.warn('Warning: Simple least squares predicts an unfeasible solution composition for {0} solution. '
                      'Recalculating with site constraints. The covariance matrix must be treated with caution.'.format(name))
        fn = lambda x, A, b, b_uncertainties: np.sqrt((A.dot(popt) - b).dot(np.linalg.solve(b_uncertainties, A.dot(popt) - b)))
        sol = minimize(fn, popt, args=(A, b, b_uncertainties), method='SLSQP',constraints=cons)
        popt = sol.x
        res = sol.fun

    if normalize:
        sump = sum(popt)
        popt /= sump
        pcov /= sump*sump
        res /= sump
    return (popt, pcov, res)
