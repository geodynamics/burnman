# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from sympy import Matrix
from scipy.linalg import block_diag

import logging
import importlib
from ..classes.polytope import MaterialPolytope, independent_row_indices
from ..classes.solution import Solution
from ..classes.composite import Composite
from .solution import transform_solution_to_new_basis

try:
    cp = importlib.import_module('cvxpy')
except ImportError as err:
    print(f'Warning: {err}. '
          'For full functionality of BurnMan, please install cvxpy.')


def solution_polytope_from_charge_balance(charges, charge_total,
                                          return_fractions=False):
    """
    Creates a polytope object from a list of the charges for each species on
    each site and the total charge for all site-species.

    Parameters
    ----------
    charges : 2D list of floats
        2D list containing the total charge for species j on site i,
        including the site multiplicity. So, for example,
        a solution with the site formula [Mg,Fe]3[Mg,Al,Si]2Si3O12 would
        have the following list: [[6., 6.], [4., 6., 8.]].

    charge_total : float
        The total charge for all site-species per formula unit.
        The example given above would have charge_total = 12.

    return_fractions : boolean
        Determines whether the created polytope object returns its
        attributes (such as endmember occupancies) as fractions or as floats.
        Default is False.

    Returns
    -------
    polytope : :class:`burnman.polytope.MaterialPolytope` object
        A polytope object corresponding to the parameters provided.
    """
    n_sites = len(charges)
    all_charges = np.concatenate(charges)
    n_site_elements = len(all_charges)
    equalities = np.empty((n_sites+1, n_site_elements+1))
    equalities[:-1, 0] = -1
    i = 0
    for i_site, site_charges in enumerate(charges):
        equalities[i_site, 1:] = [1 if (j >= i and j < i+len(site_charges))
                                  else 0 for j in range(n_site_elements)]
        i += len(site_charges)

    equalities[-1, 0] = -charge_total
    equalities[-1, 1:] = all_charges

    pos_constraints = np.concatenate((np.zeros((len(equalities[0])-1, 1)),
                                      np.identity(len(equalities[0]) - 1)),
                                     axis=1)
    return MaterialPolytope(equalities, pos_constraints,
                            return_fractions=return_fractions)


def solution_polytope_from_endmember_occupancies(endmember_occupancies,
                                                 return_fractions=False):
    """
    Creates a polytope object from a list of independent endmember occupancies.

    Parameters
    ----------
    endmember_occupancies : 2D numpy array
        2D list containing the site-species occupancies j for endmember i.
        So, for example, a solution with independent endmembers
        [Mg]3[Al]2Si3O12, [Mg]3[Mg0.5Si0.5]2Si3O12, [Fe]3[Al]2Si3O12
        might have the following array:
        [[1., 0., 1., 0., 0.],
        [1., 0., 0., 0.5, 0.5],
        [0., 1., 1., 0., 0.]],
        where the order of site-species is Mg_A, Fe_A, Al_B, Mg_B, Si_B.

    return_fractions : boolean
        Determines whether the created polytope object returns its
        attributes (such as endmember occupancies) as fractions or as floats.
        Default is False.

    Returns
    -------
    polytope : :class:`burnman.polytope.MaterialPolytope` object
        A polytope object corresponding to the parameters provided.
    """
    n_sites = sum(endmember_occupancies[0])
    n_occs = endmember_occupancies.shape[1]

    nullspace = np.array(Matrix(endmember_occupancies).nullspace(),
                         dtype=float)

    equalities = np.zeros((len(nullspace)+1, n_occs+1))
    equalities[0, 0] = -n_sites
    equalities[0, 1:] = 1

    if len(nullspace) > 0:
        try:
            equalities[1:, 1:] = nullspace
        except ValueError:
            equalities[1:, 1:] = nullspace[:, :, 0]

    pos_constraints = np.concatenate((np.zeros((len(equalities[0])-1, 1)),
                                      np.identity(len(equalities[0]) - 1)),
                                     axis=1)

    return MaterialPolytope(equalities, pos_constraints,
                            return_fractions=return_fractions,
                            independent_endmember_occupancies=endmember_occupancies)


def composite_polytope_at_constrained_composition(composite, composition,
                                                  return_fractions=False):
    """
    Creates a polytope object from a Composite object and a composition.
    This polytope describes the complete set of valid composite
    endmember amounts that satisfy the compositional constraints.

    Parameters
    ----------
    composite : :class:`burnman.Composite` object
        A composite containing one or more Solution and Mineral objects.

    composition : dictionary
        A dictionary containing the amounts of each element.

    return_fractions : boolean
        Determines whether the created polytope object returns its
        attributes (such as endmember occupancies) as fractions or as floats.
        Default is False.

    Returns
    -------
    polytope : :class:`burnman.polytope.MaterialPolytope` object
        A polytope object corresponding to the parameters provided.
    """
    c_array = np.empty((composite.n_elements, 1))
    c_array[:, 0] = [-composition[e] if e in composition else 0.
                     for e in composite.elements]

    equalities = np.concatenate((c_array, composite.stoichiometric_array.T),
                                axis=1)

    eoccs = []
    for i, ph in enumerate(composite.phases):
        if isinstance(ph, Solution):
            eoccs.append(ph.solution_model.endmember_occupancies.T)
        else:
            eoccs.append(np.ones((1, 1)))

    eoccs = block_diag(*eoccs)
    inequalities = np.concatenate((np.zeros((len(eoccs), 1)), eoccs),
                                  axis=1)

    return MaterialPolytope(equalities, inequalities,
                            number_type='float',
                            return_fractions=return_fractions)


def simplify_composite_with_composition(composite, composition):
    """
    Takes a composite and a composition, and returns the simplest composite
    object that spans the solution space at the given composition.

    For example, if the composition is given as {'Mg': 2., 'Si': 1.5, 'O': 5.},
    and the composite is given as a mix of Mg,Fe olivine and pyroxene
    solutions, this function will return a composite that only contains
    the Mg-bearing endmembers.

    Parameters
    ----------
    composite : :class:`burnman.Composite` object
        The initial Composite object

    composition : dictionary
        A dictionary containing the amounts of each element

    Returns
    -------
    simple_composite : :class:`burnman.Composite` object
        The simplified Composite object
    """
    polytope = composite_polytope_at_constrained_composition(composite,
                                                             composition,
                                                             return_fractions=True)

    composite_changed = False
    new_phases = []
    mbr_amounts = polytope.endmember_occupancies
    i = 0
    for i_ph, n_mbrs in enumerate(composite.endmembers_per_phase):
        ph = composite.phases[i_ph]

        amounts = mbr_amounts[:, i:i+n_mbrs].astype(float)
        i += n_mbrs

        rank = np.linalg.matrix_rank(amounts, tol=1.e-8)

        if rank < n_mbrs:

            if isinstance(ph, Solution) and rank > 0:

                if len(amounts) > 1:
                    c_mean = np.mean(amounts, axis=0)
                else:
                    c_mean = amounts[0]

                poly = solution_polytope_from_endmember_occupancies(
                    ph.solution_model.endmember_occupancies)
                dmbrs = poly.endmembers_as_independent_endmember_amounts

                x = cp.Variable(dmbrs.shape[0])
                objective = cp.Minimize(cp.sum_squares(x))
                constraints = [dmbrs.T@x == c_mean, x >= 0]

                prob = cp.Problem(objective, constraints)
                prob.solve()

                mbr_indices = np.argsort(x.value)[::-1]
                ind_indices = [i for i in mbr_indices
                               if x.value[i] > 1.e-6]
                new_basis = dmbrs[ind_indices]

                # And now reduce the new basis if necessary
                new_basis = new_basis[independent_row_indices(new_basis)]

                if len(new_basis) < ph.n_endmembers:
                    logging.info(f'Phase {i_ph} ({ph.name}) is '
                                 'rank-deficient ({rank} < {n_mbrs}). '
                                 'The transformed solution is described '
                                 f'using {len(new_basis)} endmembers.')

                    composite_changed = True
                    soln = transform_solution_to_new_basis(ph, new_basis)
                    new_phases.append(soln)
                else:
                    logging.info('This solution is rank-deficient '
                                 f'({rank} < {n_mbrs}), '
                                 'but its composition requires all '
                                 'independent endmembers.')
            else:
                composite_changed = True
                logging.info(f'Phase {i_ph} ({ph.name}) removed from '
                             'composite (rank = 0).')
        else:
            new_phases.append(ph)

    if composite_changed:
        return Composite(new_phases)
    else:
        return composite
