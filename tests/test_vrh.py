import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals, voigt_reuss_hill

class VRH_average(unittest.TestCase):
    def test_one_object(self):
        v = voigt_reuss_hill.vhr_average([2.0],[0.123])
        self.assertAlmostEqual(0.123, v)

    def test_two_same(self):
        v = voigt_reuss_hill.vhr_average([2.0, 2.0],[0.456, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_one_no_volume(self):
        v = voigt_reuss_hill.vhr_average([0.0, 2.0],[0.123, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_mix(self):
        v = voigt_reuss_hill.vhr_average([1.0, 2.0],[0.1, 0.2])        
        self.assertAlmostEqual(0.15833333333333, v)

class VRH(unittest.TestCase):
    def test_1(self):
        phases = [minerals.periclase()]
        phases[0].set_method('slb') 
        molar_abundances = [1.0]
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, phases, molar_abundances)
        self.assertAlmostEqual(3790.847, rho, 2)
        self.assertAlmostEqual(10330.045, v_p, 2)
        self.assertAlmostEqual(6315.856, v_s, 2)
        self.assertAlmostEqual(7315.949, v_phi, 2)
        self.assertAlmostEqual(202.898, K_vrh, 2)
        self.assertAlmostEqual(151.217, mu_vrh, 2)

    def same(self, number):
        phases = [minerals.periclase()] * number
        for p in phases:
            p.set_method('slb')
        molar_abundances = [1.0/number] * number
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, phases, molar_abundances)
        self.assertAlmostEqual(3790.847, rho, 2)
        self.assertAlmostEqual(10330.045, v_p, 2)
        self.assertAlmostEqual(6315.856, v_s, 2)
        self.assertAlmostEqual(7315.949, v_phi, 2)
        self.assertAlmostEqual(202.898, K_vrh, 2)
        self.assertAlmostEqual(151.217, mu_vrh, 2)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        phases = [minerals.periclase(), minerals.fe_perovskite()]
        for p in phases:
            p.set_method('slb')
        molar_abundances = [1.0, 0.0]
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, phases, molar_abundances)
        self.assertAlmostEqual(3790.847, rho, 2)
        self.assertAlmostEqual(10330.045, v_p, 2)
        self.assertAlmostEqual(6315.856, v_s, 2)
        self.assertAlmostEqual(7315.949, v_phi, 2)
        self.assertAlmostEqual(202.898, K_vrh, 2)
        self.assertAlmostEqual(151.217, mu_vrh, 2)


if __name__ == '__main__':
    unittest.main()
