# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman.minerals import JH_2015
from burnman.solutionbases import transform_solution_to_new_basis
from burnman.solutionbases import feasible_solution_in_component_space
from burnman.solutionbases import polytope_from_charge_balance, polytope_from_solution_model
from burnman.solutionbases import site_occupancies_to_strings, generate_complete_basis
from burnman.processanalyses import compute_and_set_phase_compositions, assemblage_affinity_misfit

grid_test = True
limits_and_grid_refinement_test = True
endmember_test = True
camph_test = True
solution_model_test = True


def float_array(arr):
    return np.array([list(map(float, v)) for v in arr])


if grid_test is True:
    # simple coupled binary solution:
    binary_polytope = polytope_from_charge_balance([[1, 2, 3]], 2,
                                                   return_fractions=True)

    grid = binary_polytope.grid(points_per_edge = 41, grid_type='site occupancies')
    print(grid[0:10])
    print(float_array(binary_polytope.independent_endmember_occupancies))
    print(len(grid))

    # simple tetrahedral solution:
    tetrahedral_polytope = polytope_from_charge_balance([[2, 2, 2, 2]], 2,
                                                        return_fractions=True)

    grid = tetrahedral_polytope.grid(points_per_edge = 41)
    print(grid[0:10])
    print(float_array(binary_polytope.independent_endmember_occupancies))
    print(len(grid))

    # simple cubical solution:
    cubical_polytope = polytope_from_charge_balance([[2,2], [1, 1], [0,0]], 3,
                                                    return_fractions=True)

    grid = cubical_polytope.grid(points_per_edge = 41)
    print(grid[0:10])
    print(float_array(binary_polytope.independent_endmember_occupancies))
    print(len(grid))


    # bridgmanite in FMASO system: Mg Fe Fe3+ Al3+ | Fe3+ Al3+ Si
    bridgmanite_polytope = polytope_from_charge_balance([[2,2,3,3], [3,3,4]], 6,
                                                        return_fractions=True)
    grid = bridgmanite_polytope.grid(points_per_edge = 5)
    print(grid)
    print(float_array(binary_polytope.independent_endmember_occupancies))
    print(len(grid))

if limits_and_grid_refinement_test is True:
    cubical_polytope = polytope_from_charge_balance([[2,2], [1, 1], [0,0]], 3,
                                                    return_fractions=True)

    e_occ = cubical_polytope.independent_endmember_occupancies.astype(float)
    print('Endmember occupancies:')
    print(e_occ)

    print('\nSolution model limits')
    print('(in form b + A*x > 0; last endmember not included):')
    print(cubical_polytope.independent_endmember_limits)

    # limits in terms of the first three independent endmembers
    endmember_limits = [[-0.1, 1., 0., 0.],
                        [-0.1, 0., 1., 0.],
                        [-0.1, 0., 0., 1.],
                        [0.3, -1., 0., 0.],
                        [0.3, 0., -1., 0.],
                        [0.3, 0., 0., -1.]]

    # limits in terms of site occupancies
    site_limits = [[-0.1, 1., 0., 0., 0., 0., 0.],
                   [-0.1, 0., 0., 1., 0., 0., 0.],
                   [-0.1, 0., 0., 0., 0., 1.,  0.],
                   [0.3, -1., 0., 0., 0., 0., 0.],
                   [0.3, 0., 0., -1., 0., 0., 0.],
                   [0.3, 0., 0., 0., 0., -1., 0.]]

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    try: # python2
        ax.set_aspect('equal')
    except:
        pass

    gridded_proportions = cubical_polytope.grid(points_per_edge=3)
    gridded_occupancies = gridded_proportions.dot(e_occ)
    ax.scatter(gridded_occupancies[:,0], gridded_occupancies[:,2], gridded_occupancies[:,4], label='global grid')

    gridded_proportions = cubical_polytope.grid(points_per_edge=5, grid_type='independent endmember proportions', limits=endmember_limits)
    gridded_occupancies = gridded_proportions.dot(e_occ)
    ax.scatter(gridded_occupancies[:,0], gridded_occupancies[:,2], gridded_occupancies[:,4], label='endmember proportion refinement')

    gridded_occupancies = cubical_polytope.grid(points_per_edge=5, grid_type='site occupancies', limits=site_limits)
    ax.scatter(gridded_occupancies[:,0], gridded_occupancies[:,2], gridded_occupancies[:,4], label='site occupancy refinement')
    plt.legend()
    plt.show()

if endmember_test is True:
    # bridgmanite in FMASO system: Mg Fe Fe3+ Al3+ | Fe3+ Al3+ Si
    bridgmanite_polytope = polytope_from_charge_balance([[2,2,3,3], [3,3,4]], 6,
                                                        return_fractions=False)
    print(bridgmanite_polytope.endmember_occupancies)

    # Two site
    simple_majorite_polytope = polytope_from_charge_balance([[2, 3, 4], [2, 3, 4]], 6,
                                                            return_fractions=False)
    print(simple_majorite_polytope.endmember_occupancies)
    print(simple_majorite_polytope.independent_endmember_occupancies)

    # Four site majorite, independent endmembers
    complex_majorite_polytope = polytope_from_charge_balance([[1, 1.5, 2], [1, 1.5, 2],
                                                              [1, 1.5, 2], [1, 1.5, 2]], 6,
                                                             return_fractions=False)
    print(complex_majorite_polytope.independent_endmember_occupancies)

    # Jamie's ?problematic? example

    disordered_px_polytope = polytope_from_charge_balance([[2, 2, 3], [2, 2, 3],[6,8]], 12,
                                                          return_fractions=False)
    print(disordered_px_polytope.independent_endmember_occupancies)

    # majorite from Holland et al. (2013)  Mg Fe Ca Na | Mg Fe Al Si
    majorite_polytope = polytope_from_charge_balance([[6, 6, 6, 3], [4, 4, 6, 8]], 12, return_fractions=False)
    print(majorite_polytope.independent_endmember_occupancies)

    # NCFMASO garnet
    # Multiplicity 3     | Multiplicity 2     | total charge
    # Ca, Fe2+, Mg2+, Na | Mg2+, Fe2+ Fe3+, Al3+, Si4+   | 12
    # SWAPPED SITES SO THAT NA ONLY IN ONE ENDMEMBER
    NCFMASO_majorite_polytope = polytope_from_charge_balance([[4, 4, 6, 6, 8],
                                                              [6, 6, 6, 3]], 12,
                                                             return_fractions=False)
    print(len(NCFMASO_majorite_polytope.independent_endmember_occupancies))
    print(NCFMASO_majorite_polytope.independent_endmember_occupancies)

    # NCFMASO garnet
    # Multiplicity 3     | Multiplicity 1     | Multiplicity 1     | total charge
    # Na, Ca, Fe2+, Mg2+ | Mg2+, Fe2+ Fe3+, Al3+, Si4+   | Mg2+, Fe2+ Fe3+, Al3+, Si4+   | 12
    NCFMASO_majorite_polytope = polytope_from_charge_balance([[3, 6, 6, 6],
                                                              [2, 2, 3, 3, 4],
                                                              [2, 2, 3, 3, 4]], 12,
                                                             return_fractions=False)
    print(len(NCFMASO_majorite_polytope.independent_endmember_occupancies))
    print(NCFMASO_majorite_polytope.independent_endmember_occupancies)

    # Spinel from Holland et al., 2018
    # Mg Fe Fe3 Al |  Mg, Fe, Fe3, Al, Cr, Ti
    spinel_polytope = polytope_from_charge_balance([[2, 2, 3, 3], [4, 4, 6, 6, 6, 8]], 8,
                                                   return_fractions=False)
    print(spinel_polytope.independent_endmember_occupancies)

    # Opx from Holland et al., 2018
    #       M1                       M2               T
    #       Mg  Fe  Al  Fe3 Cr  Ti | Mg  Fe  Ca  Na | Si  Al
    opx_polytope = polytope_from_charge_balance([[2, 2, 3, 3, 3, 4], [2, 2, 2, 1], [4, 3]], 8,
                                                return_fractions=False)
    print(opx_polytope.independent_endmember_occupancies)


if camph_test is True:
    # Camph from Diener and Powell, 2012
    # A*1    M13*3   M2*2            M4*2          T1*4
    # v Na | Mg Fe | Mg Fe Al Fe3  | Ca Mg Fe Na | Si Al
    camph_polytope = polytope_from_charge_balance([[0, 1],
                                                   [6, 6],
                                                   [4, 4, 6, 6],
                                                   [4, 4, 4, 2],
                                                   [16, 12]], 30., return_fractions=False)
    camph_array = camph_polytope.independent_endmember_occupancies
    print(camph_array)
    n_published_members = 9
    n_mbrs = len(camph_array)
    n_missing_members = n_mbrs - n_published_members


    camph_partial_basis = np.array([[1,0,1,0,1,  0,0,  0,1,0,0,0,  1,  0  ], # tr
                                    [1,0,1,0,0,  0,1,  0,1,0,0,0,  0.5,0.5], # ts
                                    [0,1,1,0,0.5,0,0.5,0,1,0,0,0,  0.5,0.5], # parg
                                    [1,0,1,0,0,  0,1,  0,0,0,0,1,  1,  0  ], # gl
                                    [1,0,1,0,1,  0,0,  0,0,1,0,0,  1,  0  ], # cumm
                                    [1,0,0,1,0,  1,0,  0,0,0,1,0,  1,  0  ], # grun
                                    [1,0,1,0,0,  1,0,  0,0,0,1,0,  1,  0  ], # a
                                    [1,0,0,1,1,  0,0,  0,0,0,1,0,  1,  0  ], # b
                                    [1,0,1,0,0,  0,0,  1,0,0,0,1,  1,  0  ]]) # mrb

    print('{0} endmembers: note that the published model has {1} endmembers'.format(n_mbrs, n_published_members))

    # Camph from Green et al., 2016; Holland et al., 2018
    # A*1      M13*3   M2*2               M4*2          T1*4    V*2
    # v Na K | Mg Fe | Mg Fe Al Fe3 Ti  | Ca Mg Fe Na | Si Al | OH O
    camph_polytope = polytope_from_charge_balance([[0, 1, 1],
                                                   [6, 6],
                                                   [4, 4, 6, 6, 8],
                                                   [4, 4, 4, 2],
                                                   [16, 12],
                                                   [-2, -4]], 28., return_fractions=False)
    camph_array = camph_polytope.independent_endmember_occupancies
    print(camph_array)
    n_published_members = 11
    n_mbrs = len(camph_array)
    n_missing_members = n_mbrs - n_published_members

    camph_partial_basis = np.array([[1,0,0,1,0,1,0,0,0,0,1,0,0,0,1,0,1,0], # tr
                                    [1,0,0,1,0,0,0,1,0,0,1,0,0,0,0.5,0.5,1,0], # ts
                                    [0,1,0,1,0,0.5,0,0.5,0,0,1,0,0,0,0.5,0.5,1,0], # parg
                                    [1,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0], # gl
                                    [1,0,0,1,0,1,0,0,0,0,0,1,0,0,1,0,1,0], # cumm
                                    [1,0,0,0,1,0,1,0,0,0,0,0,1,0,1,0,1,0], # grun
                                    [1,0,0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0], # a
                                    [1,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,1,0], # b
                                    [1,0,0,1,0,0,0,0,1,0,0,0,0,1,1,0,1,0], # mrb
                                    [0,0,1,1,0,0.5,0,0.5,0,0,1,0,0,0,0.5,0.5,1,0], # kprg
                                    [1,0,0,1,0,0,0,0,0,1,1,0,0,0,0.5,0.5,0,1]]) # tts

    print('{0} endmembers: note that the published model has {1} endmembers'.format(n_mbrs, n_published_members))
    print('The following vector(s) complete the basis:')
    print(generate_complete_basis(camph_partial_basis, camph_array)[-n_missing_members:])


if solution_model_test is True:
    # Clinopyroxene from Jennings and Holland
    cpx = JH_2015.clinopyroxene()

    cpx_polytope = polytope_from_solution_model(cpx.solution_model, return_fractions=False)

    print(cpx_polytope.independent_endmember_occupancies)

    endmember_sums = cpx_polytope.dependent_endmembers_as_independent_endmember_proportions

    all_site_occupancies = cpx_polytope.endmember_occupancies
    print(site_occupancies_to_strings(cpx.solution_model.sites,
                                      cpx.solution_model.site_multiplicities,
                                      all_site_occupancies))

    cpx.set_state(1.e5, 1000.)
    cpx.set_composition([0.5, 0.5, 0.5, 0, 0, 0, 0, -0.5])
    print(cpx.gibbs)
    new_cpx = transform_solution_to_new_basis(cpx, [[1, 1, 0, 0, 0, 0, 0, -1],
                                                    [0, 0, 1, 0, 0, 0, 0, 0]],
                                              solution_name='fdi-cats')

    new_cpx.set_state(1.e5, 1000.)
    new_cpx.set_composition([0.5, 0.5])
    print(new_cpx.gibbs)

    gt = JH_2015.garnet()
    new_gt = transform_solution_to_new_basis(gt, [[0, 1, -1, 1, 0]],
                                             endmember_names = ['skiagite'],
                                             solution_name='andr-sk')

    gt = JH_2015.garnet()
    new_gt = transform_solution_to_new_basis(gt, [[0, 0, 0, 1, 0],
                                                  [0, 1, -1, 1, 0]],
                                             endmember_names = ['andradite', 'skiagite'],
                                             solution_name='andr-sk')

    cpx = JH_2015.clinopyroxene()
    new_cpx = transform_solution_to_new_basis(cpx, [[1, 1, 0, 0, 0, 0, 0, -1],
                                                    [0, 0, 1, 0, 0, 0, 0, 0]],
                                              solution_name='fdi-cats')

    P = 1.e9
    T = 1000.
    cpx.set_composition([0.5, 0.5, 0.5, 0, 0, 0, 0, -0.5])
    new_cpx.set_composition([0.5, 0.5])
    cpx.set_state(P, T)
    new_cpx.set_state(P, T)
    print(burnman.processchemistry.formula_to_string(cpx.formula))
    print(burnman.processchemistry.formula_to_string(new_cpx.formula))
    print(cpx.gibbs, new_cpx.gibbs)
    print(cpx.S, new_cpx.S)
    print(cpx.V, new_cpx.V)
