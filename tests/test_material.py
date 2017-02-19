from __future__ import absolute_import
import unittest
import inspect
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals
from util import BurnManTest
import numpy as np

class test_material_name(BurnManTest):

    """ test Material.name and that we can edit and override it in Mineral"""

    class min_no_name(burnman.Mineral):

        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                'equation_of_state': 'slb3',
                'T_0': 300.,
                'P_0': 0.,
                'V_0': 11.24e-6,
                'K_0': 161.0e9,
                'Kprime_0': 3.8,
                'G_0': 131.0e9,
                'Gprime_0': 2.1,
                'molar_mass': .0403,
                'n': 2,
                'Debye_0': 773.,
                'grueneisen_0': 1.5,
                'q_0': 1.5,
                'eta_s_0': 2.8}
            burnman.Mineral.__init__(self)

    class min_with_name(burnman.Mineral):

        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                'name': 'name set in params',
                'equation_of_state': 'slb3',
                'T_0': 300.,
                'P_0': 0.,
                'V_0': 11.24e-6,
                'K_0': 161.0e9,
                'Kprime_0': 3.8,
                'G_0': 131.0e9,
                'Gprime_0': 2.1,
                'molar_mass': .0403,
                'n': 2,
                'Debye_0': 773.,
                'grueneisen_0': 1.5,
                'q_0': 1.5,
                'eta_s_0': 2.8}
            burnman.Mineral.__init__(self)

    class min_with_name_manually(burnman.Mineral):

        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                'equation_of_state': 'slb3',
                'T_0': 300.,
                'P_0': 0.,
                'V_0': 11.24e-6,
                'K_0': 161.0e9,
                'Kprime_0': 3.8,
                'G_0': 131.0e9,
                'Gprime_0': 2.1,
                'molar_mass': .0403,
                'n': 2,
                'Debye_0': 773.,
                'grueneisen_0': 1.5,
                'q_0': 1.5,
                'eta_s_0': 2.8}
            self.name = "manually set"
            burnman.Mineral.__init__(self)

    def test_mineral_without_name(self):
        m = self.min_no_name()
        self.assertEqual(m.name, "min_no_name")
        m.name = "bla"
        self.assertEqual(m.name, "bla")

    def test_mineral_with_name(self):

        m = self.min_with_name()
        self.assertEqual(m.name, "name set in params")
        m.name = "bla"
        self.assertEqual(m.name, "bla")

    def test_set_name_before_init(self):
        m = self.min_with_name_manually()
        self.assertEqual(m.name, "manually set")
        m.name = "bla"
        self.assertEqual(m.name, "bla")

    def test_evaluate_list(self):
        m = self.min_with_name()
        Ss = m.evaluate(['S'], [1.e5, 1.e5], [300., 300.])
        self.assertEqual(Ss[0].shape, (2,))

    def test_evaluate_1Darray(self):
        m = self.min_with_name()
        Ss = m.evaluate(['S'], np.array([1.e5, 1.e5]), np.array([300., 300.]))
        self.assertEqual(Ss[0].shape, (2,))
        
    def test_evaluate_2Darray(self):
        m = self.min_with_name()
        Ss = m.evaluate(['S'],
                        np.array([[1.e5, 1.e5, 1.e5], [1.e5, 1.e5, 1.e5]]),
                        np.array([[300., 300., 300], [300., 300., 300]]))
        self.assertEqual(Ss[0].shape, (2, 3))

        
if __name__ == '__main__':
    unittest.main()
