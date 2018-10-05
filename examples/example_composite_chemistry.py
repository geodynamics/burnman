# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_composite_chemistry
------------------

This example illustrates chemistry calculations for a composite
Material.

*Specifically uses:*

* :class:`burnman.Material`
* :func:`burnman.processchemistry.stoichiometric_matrix_and_limits`
* :func:`burnman.processchemistry.reaction_vectors`

*Demonstrates:*

* Calculation of the stoichiometric matrix and associated limits
* Calculation of an independent set of reaction vectors

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
from burnman.minerals import SLB_2011, JH_2015

if __name__ == "__main__":
    
    # 1) Sometimes we have a set of mineral compositions within a composite, and we want to
    # estimate the phase proportions. The solution will not necessarily be unique, so this
    # function should be used with caution.

    print('1) Fitting phase proportions to a bulk composition\n')
    ol = SLB_2011.mg_fe_olivine()
    opx = SLB_2011.orthopyroxene()
    assemblage = burnman.Composite([ol, opx], name='olivine-orthopyroxene mix')    
    ol.set_composition([0.9, 0.1])
    opx.set_composition([0.9, 0.1, 0.0, 0.0])

    bulk_composition = {'Mg': 3.62, 'Fe': 0.38, 'Si': 2.4, 'O': 8.8}

    amounts, residuals, rnorm = assemblage.set_potential_phase_amounts_from_bulk(bulk_composition)
    print(assemblage)
    print('Total number of moles: {0:.3f}'.format(assemblage.n_moles))
    print('Residual: {0:0.4f}'.format(rnorm))
    
    

    # 2a) In other circumstances, we might have a bulk composition, and want to
    # obtain the proportions of the endmembers within each of the solutions,
    # as well as the phase proportions.

    # In this example, we take the composition of a natural grossular-rich
    # garnet (in weight percent oxides), and find the proportions of the
    # endmembers which best fit the composition.
    
    print('\n\n2a) Fitting endmember proportions to a natural composition\n')
    
    gt = JH_2015.garnet()
    assemblage = burnman.Composite([gt], name='Natural garnet (JH2015 endmembers)')
    c = burnman.composition.OrderedCounter({'SiO2':  38.96,
                                            'TiO2':  0.71,
                                            'Al2O3': 19.93,
                                            'Cr2O3': 0.,
                                            'Fe2O3': 3.43,
                                            'FeO':   3.25,
                                            'MnO':   0.03,
                                            'MgO':   1.31,
                                            'CaO':   32.52})
    garnet_composition = burnman.Composition(c, 'weight')

    garnet_composition.renormalize('atomic', 'O', 24.)
    amounts, residuals, rnorm = assemblage.set_potential_composition_from_bulk(garnet_composition.atomic_composition,
                                                                        exact_solution_required = False,
                                                                        use_solution_guesses=False)
    

    print(assemblage)
    print('Residual: {0:0.4f}'.format(rnorm))

    # 2b) It is possible for endmembers in multisite solutions to have
    # negative proportions. BurnMan deals with this by calculating the constraints
    # on solutions and using these during the inversion.

    # In this example, we invert a garnet composition Fe3Al0.8Mg0.6Si3.6O12
    # for proportions in the system (py, alm, gr, maj, jd_maj).
    # A little noise is added to simulate microprobe compositions.
    # Oxygen is added as an unfitted element
    # (most published analyses do not measure oxygen content)
    
    
    print('\n\n2b) Fitting endmember proportions to a natural composition\n'
          '(example with negative endmember proportions)\n')

    np.random.seed(2013)

    elements = ['Fe', 'Al', 'Mg', 'Si', 'O']
    perfect_composition = np.array([3., 0.8, 0.6, 3.6, 12.])
    noise = 0.2*(np.random.rand(5) - 0.5)
    b = perfect_composition + noise 

    bulk_composition = {e: b[i] for i, e in enumerate(elements)}
    excluded_endmembers = [[4]] # exclude jadeitic majorite

    gt = SLB_2011.garnet()
    assemblage = burnman.Composite([gt], name='SLB garnet')
    amounts, residuals, rnorm = assemblage.set_potential_composition_from_bulk(bulk_composition,
                                                                               unfitted_elements='O',
                                                                               excluded_endmembers = excluded_endmembers,
                                                                               exact_solution_required = False,
                                                                               use_solution_guesses = False)
    print(assemblage)
    print('Residual: {0:0.4f}'.format(rnorm))


    
    print('\n\n3) Fitting both phase amounts and compositions to a bulk composition\n'
          '(in this case, an upper mantle assemblage)\n')
    
    ol = SLB_2011.mg_fe_olivine()
    opx = SLB_2011.orthopyroxene()
    gt = SLB_2011.garnet()

    assemblage = burnman.Composite([ol, opx, gt], name='Upper mantle assemblage')

    bulk_composition = {'Ca': 0.1, 'Fe': 0.8,
                        'Mg': 8.0, 'Al': 2.2,
                        'Si': 7.9, 'O': 28.0}

    ol.guess = [0.9, 0.1]
    opx.guess = [0.8, 0.1, 0.1, 0.]
    gt.guess = [0.9, 0.1, 0., 0.0, 0.]

    excluded_endmembers = [[], [3], [3]] # no Ca-Tschermaks or Mg-majorite
    amounts, residuals, rnorm = assemblage.set_potential_composition_from_bulk(bulk_composition,
                                                                               excluded_endmembers = excluded_endmembers,
                                                                               exact_solution_required = True,
                                                                               use_solution_guesses = True)
    print(assemblage)
    print('Residual: {0:0.4f}'.format(rnorm))


