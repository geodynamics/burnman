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
from burnman.classes.polytope import SimplexGrid


class polytope(BurnManTest):
    def test_simplex(self):
        simplex = SimplexGrid(3, 4)
        grid_list = np.array(simplex.grid(generate_type="list"))
        grid_array = simplex.grid(generate_type="array")

        self.assertArraysAlmostEqual(grid_list[-1], grid_array[-1])
        self.assertTrue(simplex.n_points() == len(grid_array))

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

        self.assertTrue(len(gt.solution_model.endmember_occupancies) == 5)
        self.assertTrue(len(gt_poly.independent_endmember_occupancies) == 5)

        self.assertTrue(len(gt.solution_model.endmember_occupancies[0]) == 9)
        self.assertTrue(len(gt_poly.endmember_occupancies[0]) == 9)
        self.assertTrue(len(gt_poly.independent_endmember_occupancies[0]) == 9)

        self.assertTrue(gt_poly.n_endmembers == 7)
        self.assertTrue(len(gt_poly.limits) == 14)

        self.assertTrue(len(gt_poly.independent_endmember_limits) == 6)
        self.assertTrue(len(gt_poly.independent_endmember_limits[0]) == 5)

        grid1 = gt_poly.grid(2, grid_type="independent endmember proportions")
        self.assertTrue(len(grid1) == 7)
        self.assertTrue(len(grid1[0]) == 5)
        grid2 = gt_poly.grid(2, grid_type="site occupancies")
        self.assertTrue(len(grid2) == 7)
        self.assertTrue(len(grid2[0]) == 9)

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
