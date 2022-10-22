from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np

from burnman import Composite
from burnman.minerals import SLB_2011
from burnman.tools.polytope import solution_polytope_from_charge_balance
from burnman.tools.polytope import solution_polytope_from_endmember_occupancies
from burnman.tools.polytope import simplify_composite_with_composition


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


if __name__ == "__main__":
    unittest.main()
