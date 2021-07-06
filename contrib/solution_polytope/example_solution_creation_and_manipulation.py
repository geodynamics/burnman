# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2019 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_solution_creation_and_manipulation
----------------

This example script demonstrates how the site-occupancy space of
a solid solution is mathematically equivalent to a polytope, and how
this insight can be used to semiautomate solution model creation
and manipulation.

*Uses:*

* :class:`burnman.solutionpolytope.SolutionPolytope`
* :func:`burnman.solutionpolytope.polytope_from_charge_balance`
* :func:`burnman.solutionpolytope.polytope_from_endmember_occupancies`
* :func:`burnman.solutionpolytope.transform_solution_to_new_basis`
* :func:`burnman.solutionpolytope.site_occupancies_to_strings`
* :func:`burnman.solutionpolytope.generate_complete_basis`
* :doc:`mineral_database`


*Demonstrates:*

* Generation of several different polytopes from charge balance constraints
* Interrogation of solution polytopes for endmember site distributions
* Generation of independent endmember sets from polytopes
* Completion of an independent endmember set from a partial set
* Transformation of endmember basis sets
"""

# Initial imports
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.solutionpolytope import polytope_from_charge_balance
from burnman.solutionpolytope import polytope_from_endmember_occupancies
from burnman.solutionpolytope import transform_solution_to_new_basis
from burnman.solutionpolytope import site_occupancies_to_strings
from burnman.solutionpolytope import generate_complete_basis
from burnman.minerals import JH_2015


if __name__ == "__main__":
    """
    First, we show how to generate a polytope from a set of species charges
    and total charge. The function we use here is polytope_from_charge_balance.
    One attribute of the returned SolutionPolytope instance is
    endmember_occupancies, which is a 2D numpy array of the
    site-species occupancies of each endmember.

    The function site_occupancies_to_strings converts the array into a
    list of strings corresponding to the endmember formulae.

    A first example is bridgmanite in the FMASO system, with formula:
    [Mg,Fe,Fe3+,Al3+][Fe3+,Al3+,Si]O3
    """

    # bridgmanite in FMASO system: Mg Fe Fe3+ Al3+ | Fe3+ Al3+ Si
    print('Endmember occupancies for bridgmanite in the FMASO system')
    print('(i.e. [Mg,Fe,Fe3+,Al3+][Fe3+,Al3+,Si]O3)')
    print('Each of the following lines represents a distinct endmember:')
    bridgmanite_polytope = polytope_from_charge_balance([[2, 2, 3, 3],
                                                         [3, 3, 4]], 6,
                                                        return_fractions=False)
    endmember_occupancies = bridgmanite_polytope.endmember_occupancies
    print(endmember_occupancies)

    print('Endmember formulae corresponding to these occupancies:')
    print(site_occupancies_to_strings([['Mg', 'Fe', 'Fef', 'Al'],
                                       ['Fef', 'Al', 'Si']],
                                      [1, 1],
                                      endmember_occupancies))

    """
    The complete set of endmembers does not depend on the nature of
    each site species. Here we calculate the endmember occupancies of
    a two-site pyrope-majorite solution (i.e. Mg3[Mg,Al,Si][Mg,Al,Si]Si3O12).
    """

    # Two site pyrope-majorite
    print('\nEndmember occupancies for two-site pyrope-majorite')
    print('(i.e. Mg3[Mg,Al,Si][Mg,Al,Si]Si3O12)')
    print('Each of the following lines represents a distinct endmember:')
    py_maj_2s_polytope = polytope_from_charge_balance([[2, 3, 4],
                                                       [2, 3, 4]], 6,
                                                      return_fractions=False)
    print(py_maj_2s_polytope.endmember_occupancies)

    """
    We can also calculate a basis set of endmembers (independent endmembers)
    from the polytope
    """
    print('A potential set of independent endmember site occupancies:')
    ind_occupancies = py_maj_2s_polytope.independent_endmember_occupancies
    print(ind_occupancies)

    print('Formulae corresponding to these independent endmember '
          'occupancies:')
    print(site_occupancies_to_strings([['Mg', 'Al', 'Si'],
                                       ['Mg', 'Al', 'Si']],
                                      [1, 1],
                                      ind_occupancies))

    print('The complete set of endmembers expressed as proportions of the '
          'independent endmember set:')
    print(py_maj_2s_polytope.dependent_endmembers_as_independent_endmember_proportions)

    """
    This majorite can be split into 4 sites, enabling ordering at the
    50% pyrope: 50% majorite composition.
    """
    # Four site majorite, independent endmembers
    print('\nIndependent endmember set for four-site pyrope-majorite')
    print('(i.e. '
          'Mg3[Mg,Al,Si]0.5[Mg,Al,Si]0.5[Mg,Al,Si]0.5[Mg,Al,Si]0.5Si3O12):')
    py_maj_4s_polytope = polytope_from_charge_balance([[1, 1.5, 2]]*4, 6,
                                                      return_fractions=False)

    occupancies = py_maj_4s_polytope.endmember_occupancies
    ind_occupancies = py_maj_4s_polytope.independent_endmember_occupancies
    print(site_occupancies_to_strings([['Mg', 'Al', 'Si']]*4, [1]*4,
                                      ind_occupancies))
    print('There are {0} endmembers in total, '
          '{1} of which are independent.'.format(len(occupancies),
                                                 len(ind_occupancies)))

    """
    Let's now look at a published model for NCFMAS majorite,
    taken from Holland et al. (2013).
    """
    print('\nIndependent endmember set for NCFMAS majorite '
          'from Holland et al., 2013')
    print('([Mg,Fe,Ca,Na]3[Mg,Fe2+,Al,Si]2Si3O12):')
    # majorite from Holland et al. (2013)  Mg Fe Ca Na | Mg Fe Al Si
    majorite_polytope = polytope_from_charge_balance([[6, 6, 6, 3],
                                                      [4, 4, 6, 8]], 12,
                                                     return_fractions=False)
    occupancies = majorite_polytope.endmember_occupancies
    ind_occupancies = majorite_polytope.independent_endmember_occupancies
    print(site_occupancies_to_strings([['Mg', 'Fe', 'Ca', 'Na'],
                                       ['Mg', 'Fe', 'Al', 'Si']],
                                      [3, 2],
                                      ind_occupancies))

    print('There are {0} endmembers in total, '
          '{1} of which are independent.'.format(len(occupancies),
                                                 len(ind_occupancies)))

    """
    We now turn to a much more complicated example; the clinoamphibole from
    Green et al. (2016) and Holland et al. (2018).
    """
    print('\n A much more complicated example...')
    print('Clinoamphibole from Green et al., 2016; Holland et al., 2018')
    print('Site species:')
    print('A*1      M13*3   M2*2               M4*2          T1*4    V*2')
    print('v Na K | Mg Fe | Mg Fe Al Fe3 Ti  | Ca Mg Fe Na | Si Al | OH O')

    site_species = [['v', 'Na', 'K'],
                    ['Mg', 'Fe'],
                    ['Mg', 'Fe', 'Al', 'Fef', 'Ti'],
                    ['Ca', 'Mg', 'Fe', 'Na'],
                    ['Si', 'Al'],
                    ['OH', 'O']]
    multiplicities = [1, 3, 2, 2, 4, 2]

    camph_polytope = polytope_from_charge_balance([[0, 1, 1],
                                                   [6, 6],
                                                   [4, 4, 6, 6, 8],
                                                   [4, 4, 4, 2],
                                                   [16, 12],
                                                   [-2, -4]], 28.,
                                                  return_fractions=True)
    camph_array = camph_polytope.independent_endmember_occupancies
    n_published_members = 11
    n_mbrs = len(camph_array)
    n_missing_members = n_mbrs - n_published_members

    print('The published model has {0} independent endmembers, '
          'but the site solution space allows for '
          '{1}.'.format(n_published_members, n_mbrs))

    """
    The published model is missing one independent endmember,
    which means that parts of the site-occupancy space are inaccessible
    to the solution. We can use the function generate_complete_basis
    along with the site occupancies of the endmembers
    reported in the original paper to calculate one of the endmembers
    which complete the basis:
    """
    camph_partial_basis = np.array([[1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0],  # tr
                                    [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 1, 0],  # ts
                                    [0, 1, 0, 1, 0, 0.5, 0, 0.5, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 1, 0],  # parg
                                    [1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0],  # gl
                                    [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0],  # cumm
                                    [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],  # grun
                                    [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],  # a
                                    [1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],  # b
                                    [1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0],  # mrb
                                    [0, 0, 1, 1, 0, 0.5, 0, 0.5, 0, 0, 1, 0, 0, 0, 0.5, 0.5, 1, 0],  # kprg
                                    [1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0.5, 0.5, 0, 1]]) # tts

    print('The following endmember is one which was not contained '
          'within the original basis:')
    complete_basis = generate_complete_basis(camph_partial_basis, camph_array)
    last_endmember = complete_basis[-n_missing_members:]

    print(site_occupancies_to_strings(site_species, multiplicities,
                                      last_endmember))

    camph_incomplete_polytope = polytope_from_endmember_occupancies(camph_partial_basis)

    print('There are {0} endmembers in the incomplete polytope, '
          'and {1} in the complete polytope.'.format(len(camph_incomplete_polytope.endmember_occupancies),
                                                     len(camph_polytope.endmember_occupancies)))

    """
    One can also list the missing endmembers (change next line to 'if True:' to print these)
    """
    if False:
        for occupancies in camph_polytope.endmember_occupancies:
            mindiff = np.min(np.sum(np.abs(camph_incomplete_polytope.endmember_occupancies - occupancies), axis=1))
            if mindiff > 1.e-6:
                print(occupancies)

    """
    The SolutionPolytope object can also be used to create a set of
    equally-spaced "pseudocompounds". This functionality is contained
    within the class function "grid", which can partition the polytope
    based on site occupancies or on proportions of independent endmembers.
    The user can choose the resolution of the grid of pseudocompounds.
    """

    print('\n\nThe SolutionPolytope object can also be used to create '
          'a set of equally-spaced pseudocompounds.')
    # simple coupled binary solution:
    binary_polytope = polytope_from_charge_balance([[1, 2, 3]], 2,
                                                   return_fractions=True)

    grid = binary_polytope.grid(points_per_edge=41,
                                grid_type='site occupancies')

    print('Here are the site-occupancies for the first 10 (out of 41) '
          'pseudocompounds for a coupled binary system: [Mg2+, Al3+, Si4+]:')
    print(grid[0:10])
    print('The endmembers of this system are:')
    print(binary_polytope.independent_endmember_occupancies.astype(float))

    print('\nWe can also grid more complex polytopes. Here is a coarse grid '
          'for bridgmanite in the FMASO system '
          '([Mg,Fe,Fe3+,Al3+][Fe3+,Al3+,Si]):')
    # bridgmanite in FMASO system: Mg Fe Fe3+ Al3+ | Fe3+ Al3+ Si
    bridgmanite_polytope = polytope_from_charge_balance([[2, 2, 3, 3],
                                                         [3, 3, 4]], 6,
                                                        return_fractions=True)
    grid = bridgmanite_polytope.grid(points_per_edge=3,
                                     grid_type='site occupancies')
    print(grid)
    print('The endmembers of this system are:')
    print(bridgmanite_polytope.independent_endmember_occupancies.astype(float))

    """
    The grid function also allows the user to set limits on the grid.
    These limits can be declared as an array in the form b + A*x > 0,
    where the first element of every row is "b", and the other columns
    correspond to the n-1 independent endmember proportions
    (with an implicit equality constraint that the
    independent endmembers must sum to 1)
    or n site occupancies. The global limits can be returned from the polytope
    via its attributes independent_endmember_limits and site_occupancy_limits
    """

    print('\nThe function also allows the user to grid parts of the polytope.')
    print('Here we plot the pseudocompounds from three different griddings of '
          'a simple cubic polytope '
          '(for a solution with three sites, no charge balance constraints)')
    cubical_polytope = polytope_from_charge_balance([[2, 2], [2, 2], [2, 2]],
                                                    6, return_fractions=True)

    e_occ = cubical_polytope.independent_endmember_occupancies.astype(float)
    print('independent endmember site occupancies:')
    print(e_occ)  # ABB, BAA, AAA, AAB

    print('\nGlobal endmember proportion limits')
    print('(in form b + A*x > 0; last endmember not included):')
    print(cubical_polytope.independent_endmember_limits)

    print('\nGlobal site occupancy limits')
    print('(in form b + A*x < 0):')
    print(cubical_polytope.site_occupancy_limits)

    try:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        try:  # python2
            ax.set_aspect('equal')
        except NotImplementedError:
            pass

        gridded_proportions = cubical_polytope.grid(points_per_edge=3)
        gridded_occupancies = gridded_proportions.dot(e_occ)
        ax.scatter(gridded_occupancies[:, 0],
                   gridded_occupancies[:, 2],
                   gridded_occupancies[:, 4], label='global grid')

        # limits in terms of the first three independent endmembers
        endmember_limits = [[-0.1, 1., 0., 0.],  # ABB > 0.1
                            [-0.1, 0., 1., 0.],  # BAA > 0.1
                            [-0.1, 0., 0., 1.],  # AAA > 0.1
                            [0.3, -1., 0., 0.],  # ABB < 0.3
                            [0.3, 0., -1., 0.],  # BAA < 0.3
                            [0.3, 0., 0., -1.]]  # AAA < 0.3
        gridded_proportions = cubical_polytope.grid(points_per_edge=5,
                                                    grid_type='independent endmember proportions',
                                                    limits=endmember_limits)
        gridded_occupancies = gridded_proportions.dot(e_occ)
        ax.scatter(gridded_occupancies[:, 0],
                   gridded_occupancies[:, 2],
                   gridded_occupancies[:, 4],
                   label='endmember proportion refinement')

        # limits in terms of site occupancies
        site_limits = [[-0.1, 1., 0., 0., 0., 0., 0.],  # A on Site 1 > 0.1
                       [-0.1, 0., 0., 1., 0., 0., 0.],  # A on Site 2 > 0.1
                       [-0.1, 0., 0., 0., 0., 1.,  0.],  # A on Site 3 > 0.1
                       [0.3, -1., 0., 0., 0., 0., 0.],  # A on Site 1 < 0.3
                       [0.3, 0., 0., -1., 0., 0., 0.],  # A on Site 2 < 0.3
                       [0.3, 0., 0., 0., 0., -1., 0.]]  # A on Site 3 < 0.3
        gridded_occupancies = cubical_polytope.grid(points_per_edge=5,
                                                    grid_type='site occupancies',
                                                    limits=site_limits)
        ax.scatter(gridded_occupancies[:, 0],
                   gridded_occupancies[:, 2],
                   gridded_occupancies[:, 4],
                   label='site occupancy refinement')
        ax.set_xlabel('p(A,Site 1)')
        ax.set_ylabel('p(A,Site 2)')
        ax.set_zlabel('p(A,Site 3)')
        plt.legend()
        plt.show()
    except ValueError:
        pass

    """
    In this last example, we show how the polytope of an instance of the
    burnman.SolidSolution class can be found, and also how the independent
    endmember set of this solution can be changed
    """
    print('\nIn this last example, we show how the polytope '
          'of an instance of the burnman.SolidSolution class can be found, '
          'and also how the independent endmember set '
          'of this solution can be changed.')
    # Clinopyroxene from Jennings and Holland, 2015
    cpx = JH_2015.clinopyroxene()
    cpx_polytope = polytope_from_endmember_occupancies(cpx.solution_model.endmember_occupancies,
                                                       return_fractions=False)

    n_ind = len(cpx_polytope.independent_endmember_occupancies)
    all_site_occupancies = cpx_polytope.endmember_occupancies
    n_all = len(all_site_occupancies)
    print('There are {0} independent endmembers in '
          'Jennings and Holland (2015) clinopyroxene and '
          '{1} endmembers in total. The site-occupancies for all '
          '{1} endmembers are:'.format(n_ind, n_all))

    occs = site_occupancies_to_strings(cpx.solution_model.sites,
                                       cpx.solution_model.site_multiplicities,
                                       all_site_occupancies)
    for occ in occs:
        print(occ)

    print('\nTransforming to a new basis containing only Fe-diopside and '
          'Ca-tschermaks...')
    new_cpx = transform_solution_to_new_basis(cpx, [[1, 1, 0, 0, 0, 0, 0, -1],
                                                    [0, 0, 1, 0, 0, 0, 0, 0]],
                                              solution_name='fdi-cats')

    P = 1.e9
    T = 1000.
    cpx.set_composition([0.5, 0.5, 0.5, 0, 0, 0, 0, -0.5])
    new_cpx.set_composition([0.5, 0.5])
    cpx.set_state(P, T)
    new_cpx.set_state(P, T)

    print('Checking that a 50:50 solution of the two endmembers has '
          'the same properties in both solutions:')
    old_formula = burnman.processchemistry.formula_to_string(cpx.formula)
    new_formula = burnman.processchemistry.formula_to_string(new_cpx.formula)
    print('Chemical formula: {0}, {1}'.format(old_formula, new_formula))
    print('Gibbs free energy (J/mol): {0:.5f}, {1:.5f}'.format(cpx.gibbs, new_cpx.gibbs))
    print('Entropy (J/K/mol): {0:.5f}, {1:.5f}'.format(cpx.S, new_cpx.S))
    print('Volume (cm^3/mol): {0:.5f}, {1:.5f}'.format(cpx.V*1.e6, new_cpx.V*1.e6))
