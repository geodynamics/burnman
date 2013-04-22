import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals, voigt_reuss_hill

class mypericlase (burnman.material):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'ref_V': 11.24e-6,
            'ref_K': 161.0e9,
            'K_prime': 3.8,
            'ref_mu': 131.0e9,
            'mu_prime': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'ref_Debye': 773.,
            'ref_grueneisen': 1.5,
            'q0': 1.5,
            'eta_0s': 2.8 }

class VRH_average(unittest.TestCase):
    def test_one_object(self):
        v = voigt_reuss_hill.vrh_average([2.0],[0.123])
        self.assertAlmostEqual(0.123, v)

    def test_two_same(self):
        v = voigt_reuss_hill.vrh_average([2.0, 2.0],[0.456, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_one_no_volume(self):
        v = voigt_reuss_hill.vrh_average([0.0, 2.0],[0.123, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_mix(self):
        v = voigt_reuss_hill.vrh_average([1.0, 2.0],[0.1, 0.2])        
        self.assertAlmostEqual(0.15833333333333, v)

class VRH(unittest.TestCase):
    def test_1(self):
        rock = burnman.composite ( ( (mypericlase(), 1.0),) )
        rock.set_method('slb3') 
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, rock)
        self.assertAlmostEqual(3791.392, rho, 2)
        self.assertAlmostEqual(10285.368, v_p, 2)
        self.assertAlmostEqual(6308.811, v_s, 2)
        self.assertAlmostEqual(7260.900, v_phi, 2)
        self.assertAlmostEqual(199.884, K_vrh/1.e9, 2)
        self.assertAlmostEqual(150.901, mu_vrh/1.e9, 2)

    def same(self, number):
        rock = burnman.composite (  [(mypericlase(), 1.0/number)]*number  )
        
        rock.set_method('slb3')
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, rock)
        self.assertAlmostEqual(3791.392, rho, 2)
        self.assertAlmostEqual(10285.368, v_p, 2)
        self.assertAlmostEqual(6308.811, v_s, 2)
        self.assertAlmostEqual(7260.900, v_phi, 2)
        self.assertAlmostEqual(199.884, K_vrh/1.e9, 2)
        self.assertAlmostEqual(150.901, mu_vrh/1.e9, 2)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.composite ( ( (minerals.periclase(), 1.0), 
                                     (minerals.fe_perovskite(), 0.0) ) )
        rock.set_method('slb3')
        rho, v_p, v_s, v_phi, K_vrh, mu_vrh = \
            voigt_reuss_hill.voigt_reuss_hill(10e9, 300, rock)
        self.assertAlmostEqual(3791.392, rho, 2)
        self.assertAlmostEqual(10285.368, v_p, 2)
        self.assertAlmostEqual(6308.811, v_s, 2)
        self.assertAlmostEqual(7260.900, v_phi, 2)
        self.assertAlmostEqual(199.884, K_vrh/1.e9, 2)
        self.assertAlmostEqual(150.901, mu_vrh/1.e9, 2)


if __name__ == '__main__':
    unittest.main()
