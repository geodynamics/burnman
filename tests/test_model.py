import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
from util import BurnManTest


class min1 (burnman.Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
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
            'eta_s_0': 2.8 }

class min2 (min1):
    def __init__(self):
        self.min1.__init__()
        self.params[''] = 0.0


class VRH_average(BurnManTest):

    def model1(self):
        rock = min1()
        rock.set_method('slb3')
        p = [40e9]
        T = [2000]
        return burnman.Model(rock, p, T, burnman.averaging_schemes.VoigtReussHill())


    def test_vs1(self):
        m = self.model1()
        self.assertArraysAlmostEqual(m.v_s(), [6359.44860613])

    def test_vp1(self):
        m = self.model1()
        self.assertArraysAlmostEqual(m.v_p(), [11043.57489092])

    def test_vphi1(self):
        m = self.model1()
        self.assertArraysAlmostEqual(m.v_phi(), [8248.46031729])

    def test_heatstuff1(self):
        m = self.model1()
        self.assertArraysAlmostEqual(m.heat_capacity_p(), [52.32168504])
        self.assertArraysAlmostEqual(m.thermal_expansivity(), [2.40018801e-05])
        self.assertArraysAlmostEqual(m.heat_capacity_v(), [49.342414])

        # reproduce by hand:
        min = m.rock
        min.set_state(m.p[0], m.T[0])
        self.assertArraysAlmostEqual(m.thermal_expansivity(), [min.thermal_expansivity()])
        self.assertArraysAlmostEqual(m.heat_capacity_v(), [min.heat_capacity_v()])
        self.assertArraysAlmostEqual(m.heat_capacity_p(), [min.heat_capacity_p()])







if __name__ == '__main__':
    unittest.main()
