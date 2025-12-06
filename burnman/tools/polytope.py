# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
from sympy import Matrix
from scipy.linalg import block_diag

import warnings
import logging
import importlib
from ..classes.polytope import MaterialPolytope, independent_row_indices
from ..classes.solution import Solution
from ..classes.composite import Composite
from .solution import transform_solution_to_new_basis

logging.captureWarnings(True)
warnings.filterwarnings(
    "ignore",
    message="The default value of raise_error will change to True in the future.",
    category=FutureWarning,
)

try:
    cp = importlib.import_module("cvxpy")
except ImportError as err:
    print(
        f"Warning: {err}. " "For full functionality of BurnMan, please install cvxpy."
    )


def solution_polytope_from_charge_balance(
    charges, charge_total, return_fractions=False
):
    """
    Creates a polytope object from a list of the charges for each species on
    each site and the total charge for all site-species.

    :param charges: 2D list containing the total charge for species j on site i,
        including the site multiplicity. So, for example,
        a solution with the site formula [Mg,Fe]3[Mg,Al,Si]2Si3O12 would
        have the following list: [[6., 6.], [4., 6., 8.]].
    :type charges: 2D list of floats
    :param charge_total: The total charge for all site-species per formula unit.
        The example given above would have charge_total = 12.
    :type charge_total: float
    :param return_fractions: Determines whether the created polytope object returns its
        attributes (such as endmember occupancies) as fractions or as floats.
        Default is False.
    :type return_fractions: bool

    :returns: A polytope object corresponding to the parameters provided.
    :rtype: :class:`burnman.polytope.MaterialPolytope` object
    """
    n_sites = len(charges)
    all_charges = np.concatenate(charges)
    n_site_elements = len(all_charges)
    equalities = np.empty((n_sites + 1, n_site_elements + 1))
    equalities[:-1, 0] = -1
    i = 0
    for i_site, site_charges in enumerate(charges):
        equalities[i_site, 1:] = [
            1 if (j >= i and j < i + len(site_charges)) else 0
            for j in range(n_site_elements)
        ]
        i += len(site_charges)

    equalities[-1, 0] = -charge_total
    equalities[-1, 1:] = all_charges

    pos_constraints = np.concatenate(
        (np.zeros((len(equalities[0]) - 1, 1)), np.identity(len(equalities[0]) - 1)),
        axis=1,
    )
    return MaterialPolytope(
        equalities, pos_constraints, return_fractions=return_fractions
    )


def solution_polytope_from_endmember_occupancies(
    endmember_occupancies, return_fractions=False
):
    """
    Creates a polytope object from a list of independent endmember occupancies.

    :param endmember_occupancies: 2D list containing the
        site-species occupancies j for endmember i.
        So, for example, a solution with independent endmembers
        [Mg]3[Al]2Si3O12, [Mg]3[Mg0.5Si0.5]2Si3O12, [Fe]3[Al]2Si3O12
        might have the following array:
        [[1., 0., 1., 0., 0.],
        [1., 0., 0., 0.5, 0.5],
        [0., 1., 1., 0., 0.]],
        where the order of site-species is Mg_A, Fe_A, Al_B, Mg_B, Si_B.
    :type endmember_occupancies: 2D numpy array

    :param return_fractions: Determines whether the created polytope object
        returns its attributes (such as endmember occupancies) as fractions
        or as floats.
    :type return_fractions: bool

    :returns: A polytope object corresponding to the parameters provided.
    :rtype: :class:`burnman.polytope.MaterialPolytope` object
    """
    n_sites = sum(endmember_occupancies[0])
    n_occs = endmember_occupancies.shape[1]

    nullspace = np.array(Matrix(endmember_occupancies).nullspace(), dtype=float)

    equalities = np.zeros((len(nullspace) + 1, n_occs + 1))
    equalities[0, 0] = -n_sites
    equalities[0, 1:] = 1

    if len(nullspace) > 0:
        try:
            equalities[1:, 1:] = nullspace
        except ValueError:
            equalities[1:, 1:] = nullspace[:, :, 0]

    pos_constraints = np.concatenate(
        (np.zeros((len(equalities[0]) - 1, 1)), np.identity(len(equalities[0]) - 1)),
        axis=1,
    )

    return MaterialPolytope(
        equalities,
        pos_constraints,
        return_fractions=return_fractions,
        independent_endmember_occupancies=endmember_occupancies,
    )


def composite_polytope_at_constrained_composition(
    composite, composition, return_fractions=False
):
    """
    Creates a polytope object from a Composite object and a composition.
    This polytope describes the complete set of valid composite
    endmember amounts that satisfy the compositional constraints.

    :param composite: A composite containing one or more Solution and Mineral objects.
    :type composite: :class:`burnman.Composite` object

    :param composition: A dictionary containing the amounts of each element.
    :type composition: dict

    :param return_fractions: Determines whether the created polytope object returns its
        attributes (such as endmember occupancies) as fractions or as floats.
    :type return_fractions: bool


    :returns: A polytope object corresponding to the parameters provided.
    :rtype: :class:`burnman.polytope.MaterialPolytope` object
    """
    c_array = np.empty((composite.n_elements, 1))
    c_array[:, 0] = [
        -composition[e] if e in composition else 0.0 for e in composite.elements
    ]

    equalities = np.concatenate((c_array, composite.stoichiometric_array.T), axis=1)

    eoccs = []
    for i, ph in enumerate(composite.phases):
        if isinstance(ph, Solution):
            eoccs.append(ph.solution_model.endmember_occupancies.T)
        else:
            eoccs.append(np.ones((1, 1)))

    eoccs = block_diag(*eoccs)
    inequalities = np.concatenate((np.zeros((len(eoccs), 1)), eoccs), axis=1)

    return MaterialPolytope(equalities, inequalities, return_fractions=return_fractions)


def reorder_to_echelon(A):
    """
    Reorders the rows of a matrix A so that it is in echelon form.
    :param A: The input matrix.
    :type A: 2D numpy array
    :returns: The reordered matrix.
    :rtype: 2D numpy array
    """
    A = np.array(A, dtype=float)
    pivot_cols = [np.argmax(r != 0) if np.any(r != 0) else np.inf for r in A]
    return A[np.argsort(pivot_cols)]


def simplify_composite_with_composition(composite, composition):
    """
    Takes a composite and a composition, and returns the simplest composite
    object that spans the solution space at the given composition.

    For example, if the composition is given as {'Mg': 2., 'Si': 1.5, 'O': 5.},
    and the composite is given as a mix of Mg,Fe olivine and pyroxene
    solutions, this function will return a composite that only contains
    the Mg-bearing endmembers.

    If the solutions have endmember proportions that are consistent with the
    bulk composition, the site occupancies of the new solution models are set to
    the values in the original models.

    :param composite: The initial Composite object.
    :type composite: :class:`burnman.Composite` object

    :param composition: A dictionary containing the amounts of each element.
    :type composition: dict

    :returns: The simplified Composite object
    :rtype: :class:`burnman.Composite` object
    """
    polytope = composite_polytope_at_constrained_composition(
        composite, composition, return_fractions=True
    )

    composite_changed = False
    new_phases = []
    mbr_amounts = polytope.endmember_occupancies
    i = 0
    for i_ph, n_mbrs in enumerate(composite.endmembers_per_phase):
        ph = composite.phases[i_ph]

        amounts = mbr_amounts[:, i : i + n_mbrs].astype(float)
        i += n_mbrs

        rank = np.linalg.matrix_rank(amounts, tol=1.0e-8)

        if rank < n_mbrs:
            if isinstance(ph, Solution) and rank > 0:
                if len(amounts) > 1:
                    c_mean = np.mean(amounts, axis=0)
                else:
                    c_mean = amounts[0]

                poly = solution_polytope_from_endmember_occupancies(
                    ph.solution_model.endmember_occupancies
                )
                dmbrs = poly.endmembers_as_independent_endmember_amounts

                x = cp.Variable(dmbrs.shape[0])
                objective = cp.Minimize(cp.sum_squares(x))
                constraints = [dmbrs.T @ x == c_mean, x >= 0]

                prob = cp.Problem(objective, constraints)
                prob.solve()

                mbr_indices = np.argsort(x.value)[::-1]
                ind_indices = [i for i in mbr_indices if x.value[i] > 1.0e-6]
                new_basis = dmbrs[ind_indices]

                # And now reduce the new basis if necessary
                new_basis = new_basis[independent_row_indices(new_basis)]

                if len(new_basis) < ph.n_endmembers:
                    logging.info(
                        f"Phase {i_ph} ({ph.name}) is "
                        "rank-deficient ({rank} < {n_mbrs}). "
                        "The transformed solution is described "
                        f"using {len(new_basis)} endmembers."
                    )

                    composite_changed = True

                    # Try to preserve the order of endmembers
                    # in the original solution model
                    # by reordering rows of the new_basis matrix
                    new_basis = reorder_to_echelon(new_basis)

                    # If the composition is set and
                    # consistent with the new basis,
                    # make a new solution with the composition
                    # already set.
                    try:
                        sol = np.linalg.lstsq(
                            new_basis.T, ph.molar_fractions, rcond=None
                        )
                        if sol[1][0] > 1.0e-10:
                            f = None
                        else:
                            f = sol[0]
                    except AttributeError:
                        f = None
                    new_name = f"{ph.name} (transformed)"
                    soln = transform_solution_to_new_basis(
                        ph, new_basis, solution_name=new_name, molar_fractions=f
                    )
                    for h, name in enumerate(soln.endmember_names):
                        if name == "User-created endmember":
                            new_name = (
                                f"Derived member (occupancies: {soln.endmembers[h][1]})"
                            )
                            soln.endmembers[h][0].name = new_name
                            soln.endmember_names[h] = new_name
                    new_phases.append(soln)
                else:
                    logging.info(
                        "This solution is rank-deficient "
                        f"({rank} < {n_mbrs}), "
                        "but its composition requires all "
                        "independent endmembers."
                    )
            else:
                composite_changed = True
                logging.info(
                    f"Phase {i_ph} ({ph.name}) removed from " "composite (rank = 0)."
                )
        else:
            new_phases.append(ph)

    if composite_changed:
        return Composite(new_phases)
    else:
        return composite


def greedy_independent_endmember_selection(
    endmember_site_occupancies, site_occupancies, small_fraction_tol=0.0, norm_tol=1e-12
):
    """
    Greedy algorithm to select independent endmembers from a solution to approximate given
    site occupancies through a non-negative linear combination of endmember site occupancies.

    This function starts with the full site occupancies and then iteratively selects endmembers
    to approximate those site occupancies. It loops through all possible endmembers in a number of steps.
    At each step the algorithm selects the endmember that can be subtracted in the largest amount from the
    current residual site occupancies without making any site occupancy negative.
    The process continues until either no endmember can be subtracted in an amount greater than fraction_tol,
    or the norm of the residual site occupancies is less than tol.

    :param endmember_site_occupancies: A 2D array of shape (m, n), where m is the number of endmembers
        and n is the number of sites. Each row corresponds to the site occupancies of an endmember.
    :type endmember_site_occupancies: np.ndarray

    :param site_occupancies: A 1D array of length n, representing the target site occupancies to approximate.
    :type site_occupancies: np.ndarray

    :param small_fraction_tol: Algorithm stops if no endmember can be added with a molar fraction larger than this value.
    :type small_fraction_tol: float, optional, default 0.0

    :param norm_tol: Algorithm stops if the norm of the residual site occupancies is less than this value.
    :type norm_tol: float, optional, default 1e-12

    :return: indices of selected endmembers, their fractions, and the final residual site occupancies.
    :rtype: tuple(list[int], list[float], np.ndarray)
    """

    site_occupancy_residuals = site_occupancies.copy()
    indices = []
    fractions = []

    for _ in range(endmember_site_occupancies.shape[0]):
        # compute fraction for each candidate endmember
        # fraction_i = min_j r_j / s_ij over s_ij > 0; if no s_ij>0 -> fraction_i = 0
        with np.errstate(divide="ignore", invalid="ignore"):
            ratios = np.where(
                endmember_site_occupancies > 0,
                site_occupancy_residuals[np.newaxis, :] / endmember_site_occupancies,
                np.inf,
            )  # shape (m,n)

        # For each row, take minimum ratio over columns where S>0; if all entries inf -> set 0
        fractions_all = np.min(ratios, axis=1)
        fractions_all[np.isinf(fractions_all)] = 0.0

        # pick largest alpha
        i_max = int(np.argmax(fractions_all))
        fraction_max = float(fractions_all[i_max])
        if fraction_max <= small_fraction_tol:
            break  # no further progress possible or desirable

        # update residual
        site_occupancy_residuals = (
            site_occupancy_residuals - fraction_max * endmember_site_occupancies[i_max]
        )

        # ensure non-negativity in elements of the residual
        site_occupancy_residuals[site_occupancy_residuals < 0] = 0.0
        indices.append(i_max)
        fractions.append(fraction_max)

        # check if done
        if np.linalg.norm(site_occupancy_residuals) <= norm_tol:
            break

    return indices, fractions, site_occupancy_residuals
