from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np

from burnman.composition import Composition



class composition(BurnManTest):

    def test_molar_no_normalize(self):

        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=False)
        self.assertArraysAlmostEqual(
            c_dict.values(), c.molar_composition.values())

    def test_molar_normalize(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=True)

        self.assertArraysAlmostEqual(np.array(c_dict.values())/sum(c_dict.values()),
                                     c.molar_composition.values())

    def test_weight_no_normalize(self):

        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='weight', normalize=False)
        self.assertArraysAlmostEqual(
            c_dict.values(), c.weight_composition.values())

    def test_weight_normalize(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='weight', normalize=True)

        self.assertArraysAlmostEqual(np.array(c_dict.values())/sum(c_dict.values()),
                                     c.weight_composition.values())

    def test_moles_to_atoms_no_normalize(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=False)

        self.assertFloatEqual(sum(c.atomic_composition.values()), 5.)

    def test_moles_to_atoms_normalize(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=True)
        self.assertFloatEqual(c.atomic_composition['O'], 0.6)

    def test_add_component_molar(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=False)
        c.add_components(c_dict, unit_type='molar')

        self.assertArraysAlmostEqual(np.array(c_dict.values())*2.,
                                     c.molar_composition.values())

    def test_add_component_weight(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='weight', normalize=False)
        c.add_components(c_dict, unit_type='weight')

        self.assertArraysAlmostEqual(np.array(c_dict.values())*2.,
                                     c.weight_composition.values())

    def test_add_component_mixed(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=False)
        c.add_components(c.weight_composition, unit_type='weight')

        self.assertArraysAlmostEqual(np.array(c_dict.values())*2.,
                                     c.molar_composition.values())

    def test_add_component_mixed_inv(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='weight', normalize=False)
        c.add_components(c.molar_composition, unit_type='molar')

        self.assertArraysAlmostEqual(np.array(c_dict.values())*2.,
                                     c.weight_composition.values())

    def test_change_component_set(self):
        c_dict = {'MgO': 1., 'SiO2': 1.}
        c = Composition(c_dict, unit_type='molar', normalize=False)
        c.change_component_set(['MgO2', 'SiO'])

        self.assertArraysAlmostEqual(np.array(c_dict.values()),
                                     c.molar_composition.values())


if __name__ == '__main__':
    unittest.main()
