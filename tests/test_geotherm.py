from __future__ import absolute_import
import unittest
from util import BurnManTest

import burnman


class mypericlase(burnman.Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """

    def __init__(self):
        self.params = {
            "equation_of_state": "slb3",
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


class geotherm(BurnManTest):
    def test_brown_shankland(self):
        T = burnman.geotherm.brown_shankland([600.0e3])
        self.assertArraysAlmostEqual(T, [1825.0])

    def test_anderson(self):
        T = burnman.geotherm.anderson([621.0e3])
        self.assertArraysAlmostEqual(T, [1805.0])

    def test_stacey_continental(self):
        T = burnman.geotherm.stacey_continental([11.0e3, 620.0e3])
        self.assertArraysAlmostEqual(T, [540.0, 2225])

    def test_stacey_oceanic(self):
        T = burnman.geotherm.stacey_oceanic([11.0e3, 620.0e3])
        self.assertArraysAlmostEqual(T, [550.0, 2225])

    def test_adiabat(self):
        rock = mypericlase()
        pressure = [100.0e9, 150.0e9]
        rock.set_method("slb3")
        T0 = 1500.0
        test_K_adiabat = burnman.geotherm.adiabatic(pressure, T0, rock)
        self.assertArraysAlmostEqual(test_K_adiabat, [1500, 1650.22034002])


if __name__ == "__main__":
    unittest.main()
