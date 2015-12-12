from __future__ import absolute_import
import unittest
import inspect
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from util import BurnManTest

# Instantiate every mineral in the mineral library
class test_material_name(BurnManTest):

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
        self.assertEquals(m.name, "min_no_name")
        m.name = "bla"
        self.assertEquals(m.name, "bla")

    def test_mineral_with_name(self):

        m = self.min_with_name()
        self.assertEquals(m.name, "name set in params")
        m.name = "bla"
        self.assertEquals(m.name, "bla")

    def test_set_name_before_init(self):
        m = self.min_with_name_manually()
        self.assertEquals(m.name, "manually set")
        m.name = "bla"
        self.assertEquals(m.name, "bla")


if __name__ == '__main__':
    unittest.main()
