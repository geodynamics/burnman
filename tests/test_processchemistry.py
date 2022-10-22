from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np

from burnman.utils.chemistry import dictionarize_formula
from burnman.utils.chemistry import sum_formulae
from burnman.utils.chemistry import formula_mass
from burnman.utils.chemistry import convert_formula
from burnman.utils.chemistry import site_occupancies_to_strings


class processchemistry(BurnManTest):
    def test_dictionarize_formula(self):
        f = dictionarize_formula("Mg0.5FeAl2SiO12")
        self.assertArraysAlmostEqual(
            [f["Mg"], f["Fe"], f["Al"], f["Si"], f["O"]], [0.5, 1.0, 2.0, 1.0, 12.0]
        )

    def test_dictionarize_formula_repeated_element(self):
        f = dictionarize_formula("Mg0.5FeAl2SiO12O2")
        self.assertArraysAlmostEqual(
            [f["Mg"], f["Fe"], f["Al"], f["Si"], f["O"]], [0.5, 1.0, 2.0, 1.0, 14.0]
        )

    def test_sum_formulae(self):
        f1 = dictionarize_formula("Mg2SiO4")
        f2 = dictionarize_formula("Fe2SiO4")
        f = sum_formulae([f1, f2], [0.25, 0.75])
        self.assertArraysAlmostEqual(
            [f["Mg"], f["Fe"], f["Si"], f["O"]], [0.5, 1.5, 1.0, 4.0]
        )

    def test_formula_mass(self):
        f = dictionarize_formula("Mg2SiO4")
        self.assertFloatEqual(formula_mass(f), 0.1406931)

    def test_convert_formula(self):
        f = dictionarize_formula("Mg2SiO4")
        f_mass = convert_formula(f, to_type="mass", normalize=False)
        self.assertArraysAlmostEqual(
            [f_mass["Mg"], f_mass["Si"], f_mass["O"]], [0.04861, 0.0280855, 0.0639976]
        )
        f_mol = convert_formula(f_mass, to_type="molar", normalize=False)
        self.assertArraysAlmostEqual([f_mol["Mg"], f_mol["Si"], f_mol["O"]], [2, 1, 4])

    def test_convert_formula_normalize(self):
        f = dictionarize_formula("Mg2SiO4")
        f_mass = convert_formula(f, to_type="mass", normalize=True)
        self.assertArraysAlmostEqual(
            [f_mass["Mg"], f_mass["Si"], f_mass["O"]], [0.345504, 0.199622, 0.454874]
        )
        f_mol = convert_formula(f_mass, to_type="molar", normalize=True)
        self.assertArraysAlmostEqual(
            [f_mol["Mg"], f_mol["Si"], f_mol["O"]], [2.0 / 7.0, 1.0 / 7.0, 4.0 / 7.0]
        )

    def test_site_occupancies_to_strings(self):
        site_names = [["Mg", "Fe", "Ca"], ["Al", "Fef"]]
        site_multiplicities = [3.0, 2.0]
        site_occupancies = np.array(
            [[1.0, 0.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 0.0, 1.0]]
        )

        strings = site_occupancies_to_strings(
            site_names, site_multiplicities, site_occupancies
        )

        self.assertEqual(strings[0], "[Mg]3[Al]2")
        self.assertEqual(strings[1], "[Fe]3[Fef]2")


if __name__ == "__main__":
    unittest.main()
