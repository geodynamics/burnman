from __future__ import absolute_import
import unittest
import inspect
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import anisotropy
from util import BurnManTest
import numpy as np

class test_anisotropy(BurnManTest):
    def test_isotropic(self):
        l, G = (0.4e11, 0.24e11)
        elastic_constants = [l, G]
        rho = 2.735e3
        m = anisotropy.AnisotropicMaterial(elastic_constants, rho, 'isotropic')

        E = G*(3.*l + 2.*G)/(l+G)
        Vp = np.sqrt((m.bulk_modulus_reuss + 4./3.*m.shear_modulus_reuss)/rho)
        Vs = np.sqrt(m.shear_modulus_reuss/rho)
        
        array_1 = [m.bulk_modulus_reuss, m.bulk_modulus_vrh,
                   m.bulk_modulus_voigt, m.shear_modulus_reuss,
                   m.shear_modulus_vrh, m.shear_modulus_voigt,
                   m.universal_elastic_anisotropy + 6., m.isotropic_poisson_ratio,
                   E, Vp, Vs, Vs]

        
        d1 = [1., 0., 0.]
        d2 = [0., 1., 0.]
    
        beta_100 = m.linear_compressibility(direction=d1)
        E_100 = m.youngs_modulus(direction=d1)
        G_100_010 = m.shear_modulus(plane_normal=d1, shear_direction=d2)
        nu_100_010 = m.poissons_ratio(axial_direction=d1, lateral_direction=d2)
        wave_speeds, wave_directions = m.wave_velocities(propagation_direction=d1)
        Vp, Vs1, Vs2 = wave_speeds

        array_2 = [1./3./beta_100, 1./3./beta_100,
                   1./3./beta_100, G_100_010,
                   G_100_010, G,
                   6., nu_100_010,
                   E_100, Vp, Vs1, Vs2]

        self.assertArraysAlmostEqual(array_1, array_2)

if __name__ == '__main__':
    unittest.main()
