import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals
import burnman.averaging_schemes as avg

class mypericlase (burnman.mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """
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

class VRH_average(unittest.TestCase):
    def test_one_object(self):
        v = avg.voigt_reuss_hill_function([2.0],[0.123])
        self.assertAlmostEqual(0.123, v)

    def test_two_same(self):
        v = avg.voigt_reuss_hill_function([2.0, 2.0],[0.456, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_one_no_volume(self):
        v = avg.voigt_reuss_hill_function([0.0, 2.0],[0.123, 0.456])        
        self.assertAlmostEqual(0.456, v)

    def test_mix(self):
        v = avg.voigt_reuss_hill_function([1.0, 2.0],[0.1, 0.2])        
        self.assertAlmostEqual(0.15833333333333, v)

class VRH(unittest.TestCase):
    def test_1(self):
        rock = burnman.composite ( ( (mypericlase(), 1.0),) )
        rock.set_method('slb3') 
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            burnman.velocities_from_rock(rock, [10e9,], [300,])
        self.assertAlmostEqual(3791.392, rho[0], 2)
        self.assertAlmostEqual(10285.368, v_p[0], 2)
        self.assertAlmostEqual(6308.811, v_s[0], 2)
        self.assertAlmostEqual(7260.900, v_phi[0], 2)
        self.assertAlmostEqual(199.884, K_vrh[0]/1.e9, 2)
        self.assertAlmostEqual(150.901, G_vrh[0]/1.e9, 2)

    def same(self, number):
        rock = burnman.composite (  [(mypericlase(), 1.0/number)]*number  )
        
        rock.set_method('slb3')
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            burnman.velocities_from_rock(rock, [10e9,], [300,])
        self.assertAlmostEqual(3791.392, rho[0], 2)
        self.assertAlmostEqual(10285.368, v_p[0], 2)
        self.assertAlmostEqual(6308.811, v_s[0], 2)
        self.assertAlmostEqual(7260.900, v_phi[0], 2)
        self.assertAlmostEqual(199.884, K_vrh[0]/1.e9, 2)
        self.assertAlmostEqual(150.901, G_vrh[0]/1.e9, 2)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.composite ( [ (minerals.SLB_2005.periclase(), 1.0), 
                                     (minerals.SLB_2005.fe_perovskite(), 0.0) ] )
        rock.set_method('slb3')
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            burnman.velocities_from_rock(rock,[10e9,], [300,])
        self.assertAlmostEqual(3791.392, rho[0], 2)
        self.assertAlmostEqual(10285.368, v_p[0], 2)
        self.assertAlmostEqual(6308.811, v_s[0], 2)
        self.assertAlmostEqual(7260.900, v_phi[0], 2)
        self.assertAlmostEqual(199.884, K_vrh[0]/1.e9, 2)
        self.assertAlmostEqual(150.901, G_vrh[0]/1.e9, 2)


if __name__ == '__main__':
    unittest.main()
