# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from .linear_fitting import weighted_constrained_least_squares


def fit_composition_to_solution(solution,
                                fitted_variables,
                                variable_values, variable_covariances,
                                variable_conversions=None,
                                normalize=True):
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
    n_vars = len(fitted_variables)
    n_mbrs = len(solution.endmembers)

    solution_variables = solution.elements
    solution_variables.extend(solution.solution_model.site_names)

    solution_matrix = np.hstack((solution.stoichiometric_matrix,
                                 solution.solution_model.endmember_noccupancies))

    n_sol_vars = solution_matrix.shape[1]

    if variable_conversions is not None:
        print(variable_conversions)
        solution_matrix = np.hstack((solution_matrix,
                                     np.zeros((solution_matrix.shape[0],
                                               len(variable_conversions)))))

        for i, (new_var, conversion_dict) in enumerate(variable_conversions.items()):
            print(new_var, conversion_dict)
            assert (new_var not in solution_variables)
            solution_variables.append(new_var)

            for var in conversion_dict.keys():
                solution_matrix[:, n_sol_vars+i] += solution_matrix[:, solution_variables.index(var)]

    # Now, construct A using the fitted variables
    A = np.zeros((n_vars, solution_matrix.shape[0]))
    for i, var in enumerate(fitted_variables):
        A[i, :] = solution_matrix[:, solution_variables.index(var)]

    b = variable_values
    Cov_b = variable_covariances

    # Define the constraints
    # Ensure that element abundances / site occupancies
    # are exactly equal to zero if the user specifies that
    # they are equal to zero.
    S, S_index = np.unique(A, axis=0, return_index=True)
    S = np.array([s for i, s in enumerate(S)
                  if np.abs(b[S_index[i]]) < 1.e-10
                  and any(np.abs(s) > 1.e-10)])
    equality_constraints = [S, np.zeros(len(S))]

    # Ensure all site occupancies are non-negative
    T = np.array([-t for t in np.unique(solution.solution_model.endmember_occupancies.T, axis=0)
                  if any(np.abs(t) > 1.e-10)])
    inequality_constraints = [T, np.zeros(len(T))]

    popt, pcov, res = weighted_constrained_least_squares(A, b, Cov_b,
                                                         equality_constraints,
                                                         inequality_constraints)

    if normalize:
        sump = sum(popt)
        popt /= sump
        pcov /= sump * sump
        res /= sump

    # Convert the variance-covariance matrix from endmember amounts to
    # endmember proportions
    dpdx = (np.eye(n_mbrs) - popt).T  # = (1. - p[i] if i == j else -p[i])
    pcov = dpdx.dot(pcov).dot(dpdx.T)
    return (popt, pcov, res)


def fit_phase_proportions_to_bulk_composition(phase_compositions,
                                              bulk_composition):

    n_phases = len(phase_compositions[0])
    inequality_constraints = [-np.eye(n_phases), np.zeros(n_phases)]
    popt, pcov, res = weighted_constrained_least_squares(phase_compositions,
                                                         bulk_composition,
                                                         None,
                                                         None,
                                                         inequality_constraints)
    return (popt, pcov, res)
