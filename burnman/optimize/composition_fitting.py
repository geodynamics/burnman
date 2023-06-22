# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from copy import deepcopy
from collections import namedtuple
from .linear_fitting import weighted_constrained_least_squares
from ..utils.chemistry import dictionarize_formula
from ..utils.chemistry import process_solution_chemistry
from ..classes.solution import Solution


class DummyCompositionSolution(Solution):
    """
    This is a dummy base class for a solution object to
    facilitate composition fitting when no solution model has
    been prepared. The model is initialized with appropriate
    chemical formulae for each endmember, and can do all basic
    compositional processing that doesn't involve any material
    properties.

    :param endmember_element_formulae: Formulae for each of the independent endmembers.
        e.g. ['Mg2SiO4', 'Fe2SiO4'].
    :type endmember_element_formulae: list of str
    :param endmember_site_formulae: Site formulae for each of the independent endmembers,
        in the same order as endmember_element_formulae. e.g. ['[Mg]2SiO4', '[Fe]2SiO4'].
    :type endmember_site_formulae: list of str
    """

    def __init__(self, endmember_element_formulae, endmember_site_formulae):
        self.endmember_formulae = [
            dictionarize_formula(f) for f in endmember_element_formulae
        ]
        self.solution_model = type(
            "Dimension", (object,), {"formulas": endmember_site_formulae}
        )()
        self.solution_model.endmembers = [None for f in endmember_site_formulae]
        process_solution_chemistry(self.solution_model)


def fit_composition_to_solution(
    solution,
    fitted_variables,
    variable_values,
    variable_covariances,
    variable_conversions=None,
    normalize=True,
):
    """
    Takes a Solution object and a set of variable names and
    associates values and covariances and finds the molar fractions of the
    solution which provide the best fit (in a least-squares sense)
    to the variable values.

    The fitting applies appropriate non-negativity constraints
    (i.e. no species can have a negative occupancy on a site).

    :param solution: The solution to use in the fitting procedure.
    :type solution: :class:`burnman.Solution`

    :param fitted_variables: A list of the variables used to find the best-fit
        molar fractions of the solution. These should either be elements
        such as "Fe", site_species such as "Fef_B" which would correspond to a
        species labelled Fef on the second site,
        or user-defined variables which are arithmetic sums of
        elements and/or site_species defined in "variable_conversions".
    :type fitted_variables: list of str

    :param variable_values: Numerical values of the fitted variables.
        These should be given as amounts; they do not need to be normalized.
    :type variable_values: numpy.array

    :param variable_covariances: Covariance matrix of the variables.
    :type variable_covariances: 2D numpy.array

    :param variable_conversions: A dictionary converting any user-defined variables
        into an arithmetic sum of element and site-species amounts. For example,
        {'Mg_equal': {'Mg_A': 1., 'Mg_B': -1.}}, coupled with Mg_equal = 0
        would impose a constraint that the amount of Mg would be equal on
        the first and second site in the solution.
    :type variable_conversions: dict of dict, or None

    :param normalize: If True, normalizes the optimized molar fractions to sum to unity.
    :type normalize: bool

    :returns: Optimized molar fractions, corresponding covariance matrix and the
        weighted residual.
    :rtype: tuple of 1D numpy.array, 2D numpy.array and float
    """

    solution_variables = deepcopy(solution.elements)
    solution_variables.extend(solution.solution_model.site_names)

    solution_matrix = np.hstack(
        (solution.stoichiometric_matrix, solution.solution_model.endmember_noccupancies)
    )

    n_vars = len(fitted_variables)
    n_mbrs, n_sol_vars = solution_matrix.shape

    if variable_conversions is not None:
        solution_matrix = np.hstack(
            (
                solution_matrix,
                np.zeros((solution_matrix.shape[0], len(variable_conversions))),
            )
        )

        for i, (new_var, conversion_dict) in enumerate(variable_conversions.items()):
            assert new_var not in solution_variables
            solution_variables.append(new_var)

            for var in conversion_dict.keys():
                solution_matrix[:, n_sol_vars + i] += solution_matrix[
                    :, solution_variables.index(var)
                ]

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
    S = np.array(
        [
            s
            for i, s in enumerate(S)
            if np.abs(b[S_index[i]]) < 1.0e-10 and any(np.abs(s) > 1.0e-10)
        ]
    )
    equality_constraints = [S, np.zeros(len(S))]

    # Ensure all site occupancies are non-negative
    T = np.array(
        [
            -t
            for t in np.unique(solution.solution_model.endmember_occupancies.T, axis=0)
            if any(np.abs(t) > 1.0e-10)
        ]
    )
    inequality_constraints = [T, np.zeros(len(T))]

    popt, pcov, res = weighted_constrained_least_squares(
        A, b, Cov_b, equality_constraints, inequality_constraints
    )

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


def fit_phase_proportions_to_bulk_composition(phase_compositions, bulk_composition):
    """
    Performs weighted constrained least squares on a set of phase compositions
    to find the amount of those phases that best-fits a given bulk composition.

    The fitting applies appropriate non-negativity constraints
    (i.e. no phase can have a negative abundance in the bulk).

    :param phase_compositions: The composition of each phase.
        Can be in weight or mole amounts.
    :type phase_compositions: 2D numpy.array

    :param bulk_composition: The bulk composition of the composite.
        Must be in the same units as the phase compositions.
    :type bulk_composition: numpy.array

    :returns: Optimized molar fractions, corresponding covariance matrix and the
        weighted residual.
    :rtype: tuple of 1D numpy.array, 2D numpy.array and float
    """

    n_phases = len(phase_compositions[0])
    inequality_constraints = [-np.eye(n_phases), np.zeros(n_phases)]
    popt, pcov, res = weighted_constrained_least_squares(
        phase_compositions, bulk_composition, None, None, inequality_constraints
    )
    return (popt, pcov, res)
