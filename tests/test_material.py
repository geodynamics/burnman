from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np

import burnman


class test_material_name(BurnManTest):
    """test Material.name and that we can edit and override it in Mineral"""

    class min_no_name(burnman.Mineral):
        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                "equation_of_state": "slb3",
                "T_0": 300.0,
                "P_0": 0.0,
                "V_0": 11.24e-6,
                "K_0": 161.0e9,
                "Kprime_0": 3.8,
                "G_0": 131.0e9,
                "Gprime_0": 2.1,
                "molar_mass": 0.0403,
                "n": 2,
                "Debye_0": 773.0,
                "grueneisen_0": 1.5,
                "q_0": 1.5,
                "eta_s_0": 2.8,
            }
            burnman.Mineral.__init__(self)

    class min_with_name(burnman.Mineral):
        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                "name": "name set in params",
                "equation_of_state": "slb3",
                "T_0": 300.0,
                "P_0": 0.0,
                "V_0": 11.24e-6,
                "K_0": 161.0e9,
                "Kprime_0": 3.8,
                "G_0": 131.0e9,
                "Gprime_0": 2.1,
                "molar_mass": 0.0403,
                "n": 2,
                "Debye_0": 773.0,
                "grueneisen_0": 1.5,
                "q_0": 1.5,
                "eta_s_0": 2.8,
            }
            burnman.Mineral.__init__(self)

    class min_with_name_manually(burnman.Mineral):
        """
        Stixrude & Lithgow-Bertelloni 2005 and references therein
        """

        def __init__(self):
            self.params = {
                "equation_of_state": "slb3",
                "T_0": 300.0,
                "P_0": 0.0,
                "V_0": 11.24e-6,
                "K_0": 161.0e9,
                "Kprime_0": 3.8,
                "G_0": 131.0e9,
                "Gprime_0": 2.1,
                "molar_mass": 0.0403,
                "n": 2,
                "Debye_0": 773.0,
                "grueneisen_0": 1.5,
                "q_0": 1.5,
                "eta_s_0": 2.8,
            }
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
        Ss = m.evaluate(["S"], [1.0e5, 1.0e5], [300.0, 300.0])
        self.assertEqual(Ss[0].shape, (2,))

    def test_evaluate_1Darray(self):
        m = self.min_with_name()
        Ss = m.evaluate(["S"], np.array([1.0e5, 1.0e5]), np.array([300.0, 300.0]))
        self.assertEqual(Ss[0].shape, (2,))

    def test_evaluate_2Darray(self):
        m = self.min_with_name()
        Ss = m.evaluate(
            ["S"],
            np.array([[1.0e5, 1.0e5, 1.0e5], [1.0e5, 1.0e5, 1.0e5]]),
            np.array([[300.0, 300.0, 300], [300.0, 300.0, 300]]),
        )
        self.assertEqual(Ss[0].shape, (2, 3))

    def test_set_state_with_volume(self):
        m = self.min_with_name()
        P0 = 6.0e9
        T0 = 1000.0
        T1 = 298.15

        m.set_state(P0, T0)
        V = m.V
        m.set_state_with_volume(V, T1)
        P1 = m.pressure
        m.set_state(P0, T0)  # forget new state
        m.set_state(P1, T1)  # return to new state
        self.assertFloatEqual(V, m.V)


if __name__ == "__main__":
    unittest.main()
