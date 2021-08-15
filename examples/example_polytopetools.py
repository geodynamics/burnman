# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_polytopetools
---------------------

This example script demonstrates how to use the linear algebra functions
included in the polytopetools module. In particular, it shows how a
user can automatically calculate the endmembers bounding a solution or a
composite polytope.

*Uses:*

* :class:`burnman.polytope.MaterialPolytope`
* :func:`burnman.polytope.generate_complete_basis`
* :func:`burnman.polytopetools.solution_polytope_from_charge_balance`
* :func:`burnman.polytopetools.solution_polytope_from_endmember_occupancies`
* :func:`burnman.processchemistry.site_occupancies_to_strings`
* :func:`burnman.polytopetools.composite_polytope_at_constrained_composition`
* :func:`burnman.polytopetools.simplify_composite_with_composition`
* :doc:`mineral_database`


*Demonstrates:*

* Generation of several different polytopes from charge balance constraints
* Interrogation of solution polytopes for endmember site distributions
* Generation of independent endmember sets from polytopes
* Completion of an independent endmember set from a partial set
* Simplifying composites at a given bulk composition
"""
from __future__ import absolute_import
import numpy as np

import burnman_path  # adds the local burnman directory to the path
import burnman
from burnman.minerals import SLB_2011
from burnman.polytopetools import solution_polytope_from_charge_balance
from burnman.polytopetools import solution_polytope_from_endmember_occupancies
from burnman.polytopetools import composite_polytope_at_constrained_composition as composite_polytope
from burnman.polytopetools import simplify_composite_with_composition
from burnman.processchemistry import site_occupancies_to_strings
from burnman.processchemistry import formula_to_string

assert burnman_path  # silence pyflakes warning

if __name__ == "__main__":
    """
    Part 1: Solution models

    First, we show how to generate a polytope from a set of species charges
    and total charge. The function we use here is
    solution_polytope_from_charge_balance.

    One attribute of the returned MaterialPolytope instance is
    endmember_occupancies, which is a 2D numpy array of the
    site-species occupancies of each endmember.

    The function site_occupancies_to_strings converts the array into a
    list of strings corresponding to the endmember formulae.

    A first example is bridgmanite in the FMASO system, with formula:
    [Mg,Fe,Fe3+,Al3+][Fe3+,Al3+,Si]O3
    """

    # bridgmanite in FMASO system: Mg Fe Fe3+ Al3+ | Fe3+ Al3+ Si
    print('The polytope tools in BurnMan can be used to calculate '
          'all potential endmembers in a solid solution.')
    print('\nEndmember occupancies for bridgmanite in the FMASO system')
    print('([Mg,Fe,Fe3+,Al3+][Fe3+,Al3+,Si]O3)')
    print('Each of the following lines represents a distinct endmember:')
    bdg_poly = solution_polytope_from_charge_balance([[2, 2, 3, 3],
                                                      [3, 3, 4]], 6,
                                                     return_fractions=False)
    endmember_occupancies = bdg_poly.endmember_occupancies
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
    py_maj_poly = solution_polytope_from_charge_balance([[2, 3, 4],
                                                         [2, 3, 4]], 6,
                                                        return_fractions=False)
    print(py_maj_poly.endmember_occupancies)

    """
    We can also calculate a basis set of endmembers (independent endmembers)
    from the polytope
    """
    print('A potential set of independent endmember site occupancies:')
    ind_occupancies = py_maj_poly.independent_endmember_occupancies
    print(ind_occupancies)

    print('Formulae corresponding to these independent endmember '
          'occupancies:')
    print(site_occupancies_to_strings([['Mg', 'Al', 'Si'],
                                       ['Mg', 'Al', 'Si']],
                                      [1, 1],
                                      ind_occupancies))

    print('The complete set of endmembers expressed as proportions of the '
          'independent endmember set:')
    print(py_maj_poly.endmembers_as_independent_endmember_amounts)

    print('\nWe can also create a polytope using the site-species occupancies '
          'of a set of independent endmembers of a solution. '
          'Here we do this for the garnet solid solution provided in the '
          'Stixrude and Lithgow-Bertelloni database.')

    gt = SLB_2011.garnet()

    print('Endmember compositions:')
    for f in gt.endmember_formulae:
        print(formula_to_string(f))

    gt_poly = solution_polytope_from_endmember_occupancies(
        gt.solution_model.endmember_occupancies, return_fractions=True)

    assert np.all(np.abs(gt.solution_model.endmember_occupancies
                         - gt_poly.independent_endmember_occupancies) < 1.e-5)

    print('Endmembers as independent endmember amounts:')
    print(gt_poly.endmembers_as_independent_endmember_amounts)

    strs = site_occupancies_to_strings(gt.solution_model.sites,
                                       gt.solution_model.site_multiplicities,
                                       gt_poly.endmember_occupancies)

    print('Site occupancies for all endmembers:')
    for string in strs:
        print(string)

    """
    Part 2: Composites
    """
    print('\n\nComposites are also constrained by linear equalities and '
          'inequalities. Here, we demonstrate how BurnMan can take a '
          'composite and a bulk composition, and simplify the composite so '
          'that it only contains endmembers which can have non-zero '
          'quantities.')

    gt = SLB_2011.garnet()
    ol = SLB_2011.mg_fe_olivine()

    assemblage = burnman.Composite([gt, ol])
    composition = {'Fe': 2, 'Mg': 18, 'Si': 15, 'O': 50}

    print('\nThe starting assemblage includes NCFMAS garnet and FMS olivine.')
    print(f'Starting composition: {formula_to_string(composition)}')

    new_assemblage = simplify_composite_with_composition(assemblage,
                                                         composition)
    if new_assemblage is not assemblage:
        print('The simplifying function has returned a modified assemblage.')
    old_polytope = composite_polytope(assemblage, composition,
                                      return_fractions=True)
    new_polytope = composite_polytope(new_assemblage, composition,
                                      return_fractions=True)

    print('\nThe old assemblage included endmembers with the '
          'following compositions:')
    print([formula_to_string(f) for f in assemblage.endmember_formulae])

    print('Which led to the following extreme vertices satifying '
          'the bulk compositional constraints:')
    print(old_polytope.endmember_occupancies)

    print('\nThe new assemblage has fewer endmembers, with compositions:')
    print([formula_to_string(f) for f in new_assemblage.endmember_formulae])

    print('The same vertices now look like this:')
    print(new_polytope.endmember_occupancies)
