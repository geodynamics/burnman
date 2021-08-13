# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_solution_creation_and_manipulation
----------------

This example script demonstrates how the site-occupancy space of
a solid solution is mathematically equivalent to a polytope, and how
this insight can be used to semi-automate the creation and
manipulation of solution models.

*Uses:*

* :class:`burnman.polytope.MaterialPolytope`
* :func:`burnman.polytope.polytope_from_charge_balance`
* :func:`burnman.polytope.polytope_from_endmember_occupancies`
* :func:`burnman.polytope.transform_solution_to_new_basis`
* :func:`burnman.polytope.site_occupancies_to_strings`
* :func:`burnman.polytope.generate_complete_basis`
* :doc:`mineral_database`


*Demonstrates:*

* Generation of several different polytopes from charge balance constraints
* Interrogation of solution polytopes for endmember site distributions
* Generation of independent endmember sets from polytopes
* Completion of an independent endmember set from a partial set
* Transformation of endmember basis sets
"""
from __future__ import absolute_import
import numpy as np
import matplotlib.pyplot as plt

import burnman_path  # adds the local burnman directory to the path
import burnman
from burnman.polytope import polytope_from_composite
from burnman.polytope import polytope_from_charge_balance
from burnman.polytope import polytope_from_endmember_occupancies
from burnman.polytope import transform_solution_to_new_basis
from burnman.polytope import site_occupancies_to_strings
from burnman.polytope import generate_complete_basis
from burnman.polytope import independent_row_indices
from burnman.minerals import SLB_2011
from burnman.polytope import simplify_composite_with_composition

assert burnman_path  # silence pyflakes warning

if __name__ == "__main__":

    gt = SLB_2011.garnet()
    ol = SLB_2011.mg_fe_olivine()
    wad = SLB_2011.mg_fe_wadsleyite()
    opx = SLB_2011.orthopyroxene()
    stv = SLB_2011.stishovite()

    ol.set_composition([0.93, 0.07])
    wad.set_composition([0.93, 0.07])
    opx.set_composition([0.8, 0.1, 0.05, 0.05])
    gt.set_composition([0.8, 0.1, 0.05, 0.03, 0.02])

    assemblage = burnman.Composite([ol, opx, gt, stv], [0.7, 0.3, 0., 0.])
    composition = {'Na': 0.02, 'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9,
                   'Ca': 0.2, 'Al': 0.4, 'O': 6.81}

    assemblage2 = burnman.Composite([ol, wad, stv], [0.7, 0.3, 0.])
    composition2 = {'Fe': 0.2, 'Mg': 1.8, 'Si': 2., 'O': 6.}

    assemblage3 = burnman.Composite([ol, gt], [0.7, 0.3])
    composition3 = {'Fe': 0.2, 'Mg': 1.8, 'Si': 2., 'O': 6.}

    assemblage4 = assemblage3
    composition4 = {'Fe': 3., 'Mg': 1., 'Si': 3.9, 'O': 11.8}

    a, c = (assemblage4, composition4)
    new_assemblage = simplify_composite_with_composition(a, c)
    if new_assemblage is not a:
        print('assemblage simplified')
    old_polytope = polytope_from_composite(a, c, return_fractions=True)
    new_polytope = polytope_from_composite(new_assemblage, c,
                                           return_fractions=True)

    print(old_polytope.endmember_occupancies)
    print(new_polytope.endmember_occupancies)
