from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np
from sympy import Matrix
from fractions import Fraction
import importlib

from burnman import Composite
from burnman.minerals import SLB_2011, JH_2015
from burnman.tools.polytope import solution_polytope_from_charge_balance
from burnman.tools.polytope import solution_polytope_from_endmember_occupancies
from burnman.tools.polytope import simplify_composite_with_composition
from burnman import MaterialPolytope


class polytope(BurnManTest):
    def test_polytope_from_charge_balance(self):
        bdg_poly = solution_polytope_from_charge_balance(
            [[2, 2, 3, 3], [3, 3, 4]], 6, return_fractions=False
        )
        self.assertTrue(bdg_poly.endmember_occupancies.shape[0] == 6)
        self.assertTrue(bdg_poly.endmember_occupancies.shape[1] == 7)

    def test_solution_polytope_from_endmember_occupancies(self):
        gt = SLB_2011.garnet()
        gt_poly = solution_polytope_from_endmember_occupancies(
            gt.solution_model.endmember_occupancies, return_fractions=True
        )

        self.assertTrue(
            np.all(
                np.abs(
                    gt.solution_model.endmember_occupancies
                    - gt_poly.independent_endmember_occupancies
                )
                < 1.0e-5
            )
        )

    def test_simplify_composite(self):
        gt = SLB_2011.garnet()
        ol = SLB_2011.mg_fe_olivine()
        assemblage = Composite([ol, gt], [0.7, 0.3])
        composition = {"Fe": 3.0, "Mg": 1.0, "Si": 3.9, "O": 11.8}

        new_assemblage = simplify_composite_with_composition(assemblage, composition)
        new_gt = new_assemblage.phases[1]
        self.assertTrue(new_gt.n_endmembers == 2)
        strings = list(sorted([e[1] for e in new_gt.endmembers]))
        self.assertEqual(strings[0], "[Fe]3[Mg][Si]")
        self.assertEqual(strings[1], "[Mg]3[Mg][Si]")

    def test_simplify_composite_and_composition(self):
        gt = SLB_2011.garnet()
        gt.set_composition([-0.1, 0.1, 0.0, 1.0, 0.0])
        ol = SLB_2011.mg_fe_olivine()
        assemblage = Composite([ol, gt], [0.7, 0.3])
        composition = {"Fe": 3.0, "Mg": 1.0, "Si": 3.9, "O": 11.8}

        new_assemblage = simplify_composite_with_composition(assemblage, composition)
        new_gt = new_assemblage.phases[1]
        self.assertTrue(new_gt.n_endmembers == 2)
        strings = list(sorted([e[1] for e in new_gt.endmembers]))
        self.assertEqual(strings[0], "[Fe]3[Mg][Si]")
        self.assertEqual(strings[1], "[Mg]3[Mg][Si]")
        self.assertArraysAlmostEqual([0.1, 0.9], new_gt.molar_fractions)

    def test_cddlib_versions(self):
        gt = JH_2015.garnet()
        endmember_occupancies = gt.solution_model.endmember_occupancies

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
            (
                np.zeros((len(equalities[0]) - 1, 1)),
                np.identity(len(equalities[0]) - 1),
            ),
            axis=1,
        )

        equalities = np.array([[Fraction(v) for v in r] for r in equalities])
        pos_constraints = np.array([[Fraction(v) for v in r] for r in pos_constraints])

        poly = MaterialPolytope(
            equalities,
            pos_constraints,
            independent_endmember_occupancies=endmember_occupancies,
        )

        try:
            _ = importlib.import_module("cdd.gmp")
            self.assertTrue(type(poly.raw_vertices[0][0]) is Fraction)
        except ImportError:
            self.assertTrue(type(poly.raw_vertices[0][0]) is float)


if __name__ == "__main__":
    unittest.main()
