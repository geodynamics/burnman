import unittest
from util import BurnManTest
import numpy as np

from burnman.classes.composition import Composition


class composition(BurnManTest):
    def test_molar_no_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        self.assertArraysAlmostEqual(
            list(c_dict.values()), c.molar_composition.values()
        )

    def test_molar_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=True)

        c_fractions = np.array(list(c_dict.values())) / sum(list(c_dict.values()))
        self.assertArraysAlmostEqual(c_fractions, c.molar_composition.values())

    def test_weight_no_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="weight", normalize=False)
        self.assertArraysAlmostEqual(c_dict.values(), c.weight_composition.values())

    def test_weight_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="weight", normalize=True)

        c_fractions = np.array(list(c_dict.values())) / sum(list(c_dict.values()))
        self.assertArraysAlmostEqual(c_fractions, c.weight_composition.values())

    def test_moles_to_atoms_no_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        self.assertFloatEqual(sum(c.atomic_composition.values()), 5.0)

    def test_moles_to_atoms_normalize(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=True)
        self.assertFloatEqual(c.atomic_composition["O"], 1.5)

    def test_add_component_molar(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c.add_components(c_dict, unit_type="molar")

        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())) * 2.0, c.molar_composition.values()
        )

    def test_add_component_weight(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="weight", normalize=False)
        c.add_components(c_dict, unit_type="weight")

        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())) * 2.0, c.weight_composition.values()
        )

    def test_add_component_mixed(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c.add_components(c.weight_composition, unit_type="weight")
        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())) * 2.0, c.molar_composition.values()
        )

    def test_add_component_mixed_inv(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="weight", normalize=False)
        c.add_components(c.molar_composition, unit_type="molar")

        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())) * 2.0, c.weight_composition.values()
        )

    def test_change_component_set(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c.change_component_set(["MgO2", "SiO"])

        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())), c.molar_composition.values()
        )

    def test_negative_composition(self):
        c_dict = {"MgO": -1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        self.assertArraysAlmostEqual(
            np.array(list(c_dict.values())), c.molar_composition.values()
        )

    def test_add_two_compositions(self):
        c_dict_1 = {"MgO": 1.0, "SiO2": 1.0}
        c_dict_2 = {"MgO": 2.0, "FeO": 3.0}
        c1 = Composition(c_dict_1, unit_type="molar", normalize=False)
        c2 = Composition(c_dict_2, unit_type="molar", normalize=False)

        c3 = c1 + c2

        self.assertArraysAlmostEqual(
            np.array([3.0, 1.0, 3.0]),
            [c3.molar_composition[key] for key in ["MgO", "SiO2", "FeO"]],
        )

        self.assertArraysAlmostEqual(
            np.array([1.0, 1.0]), [c1.molar_composition[key] for key in ["MgO", "SiO2"]]
        )
        self.assertArraysAlmostEqual(
            np.array([2.0, 3.0]), [c2.molar_composition[key] for key in ["MgO", "FeO"]]
        )

    def test_difference_two_compositions(self):
        c_dict_1 = {"MgO": 1.0, "SiO2": 1.0}
        c_dict_2 = {"MgO": 2.0, "FeO": 3.0}
        c1 = Composition(c_dict_1, unit_type="molar", normalize=False)
        c2 = Composition(c_dict_2, unit_type="molar", normalize=False)

        c3 = c1 - c2

        self.assertArraysAlmostEqual(
            np.array([-1.0, 1.0, -3.0]),
            [c3.molar_composition[key] for key in ["MgO", "SiO2", "FeO"]],
        )

        self.assertArraysAlmostEqual(
            np.array([1.0, 1.0]), [c1.molar_composition[key] for key in ["MgO", "SiO2"]]
        )
        self.assertArraysAlmostEqual(
            np.array([2.0, 3.0]), [c2.molar_composition[key] for key in ["MgO", "FeO"]]
        )

    def test_multiply_composition(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)

        c2 = c * 2.0

        self.assertArraysAlmostEqual(
            np.array([2.0, 2.0]), [c2.molar_composition[key] for key in ["MgO", "SiO2"]]
        )

        self.assertArraysAlmostEqual(
            np.array([1.0, 1.0]), [c.molar_composition[key] for key in ["MgO", "SiO2"]]
        )

    def test_divide_composition(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)

        c2 = c / 2.0

        self.assertArraysAlmostEqual(
            np.array([0.5, 0.5]), [c2.molar_composition[key] for key in ["MgO", "SiO2"]]
        )

        self.assertArraysAlmostEqual(
            np.array([1.0, 1.0]), [c.molar_composition[key] for key in ["MgO", "SiO2"]]
        )

    def test_add_in_place(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        d_dict = {"FeO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c += Composition(d_dict, unit_type="molar", normalize=False)

        self.assertArraysAlmostEqual(
            np.array([1.0, 2.0, 1.0]),
            [c.molar_composition[key] for key in ["MgO", "SiO2", "FeO"]],
        )

    def test_subtract_in_place(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        d_dict = {"FeO": 1.0, "SiO2": 2.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c -= Composition(d_dict, unit_type="molar", normalize=False)

        self.assertArraysAlmostEqual(
            np.array([1.0, -1.0, -1.0]),
            [c.molar_composition[key] for key in ["MgO", "SiO2", "FeO"]],
        )

    def test_multiply_in_place(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c *= 2.0

        self.assertArraysAlmostEqual(
            np.array([2.0, 2.0]), [c.molar_composition[key] for key in ["MgO", "SiO2"]]
        )

    def test_divide_in_place(self):
        c_dict = {"MgO": 1.0, "SiO2": 1.0}
        c = Composition(c_dict, unit_type="molar", normalize=False)
        c /= 2.0

        self.assertArraysAlmostEqual(
            np.array([0.5, 0.5]), [c.molar_composition[key] for key in ["MgO", "SiO2"]]
        )


if __name__ == "__main__":
    unittest.main()
