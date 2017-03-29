# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_anisotropy
-----------------

This example illustrates the basic functions required to convert 
an elastic stiffness tensor into elastic properties.

*Specifically uses:*

* :class:`burnman.AnisotropicMaterial`

*Demonstrates:*

* anisotropic functions

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import anisotropy 

if __name__ == "__main__":
    # Let's first look at an isotropic material
    # As an example, let's use the lambda (C12), mu (C44) and rho
    # for flow basalt given by Christensen et al. (1980)
    # Initial Reports of the Deep Sea Drilling Project 59: 515-17.

    elastic_constants = [0.4e11, 0.24e11]
    rho = 2.735e3
    basalt = anisotropy.AnisotropicMaterial(elastic_constants, rho, 'isotropic')

    print('Basalt isotropic elastic properties:\n')
    print('Bulk modulus bounds: {0:.3e} {1:.3e} {2:.3e}'.format(basalt.bulk_modulus_reuss,
                                                                basalt.bulk_modulus_vrh,
                                                                basalt.bulk_modulus_voigt))
    print('Shear modulus bounds: {0:.3e} {1:.3e} {2:.3e}'.format(basalt.shear_modulus_reuss,
                                                                 basalt.shear_modulus_vrh,
                                                                 basalt.shear_modulus_voigt))
    print('Universal elastic anisotropy: {0:.4f}\n'
          'Isotropic poisson ratio: {1:.4f}\n'.format(basalt.universal_elastic_anisotropy,
                                                      basalt.isotropic_poisson_ratio))
    
    d1 = [1., 0., 0.]
    d2 = [0., 1., 0.]
    
    beta_100 = basalt.linear_compressibility(direction=d1)
    E_100 = basalt.youngs_modulus(direction=d1)
    G_100_010 = basalt.shear_modulus(plane_normal=d1, shear_direction=d2)
    nu_100_010 = basalt.poissons_ratio(axial_direction=d1, lateral_direction=d2)
    wave_speeds, wave_directions = basalt.wave_velocities(propagation_direction=d1)
    Vp, Vs1, Vs2 = wave_speeds
    
    print('Linear compressibility along 100: {0:.3e}\n'
          'Young\'s modulus along 100: {1:.3e}\n'
          'Shear modulus on 100 plane in direction 010: {2:.3e}\n'
          'Poisson ratio for 100/010: {3}\n'
          'Vp, Vs1, Vs2 (km/s): '
          '{4:.2f}, {5:.2f}, {6:.2f}\n'.format(beta_100, E_100,
                                               G_100_010, nu_100_010,
                                               Vp/1.e3, Vs1/1.e3, Vs2/1.e3))
    
    
    # The next three lines initialise the AnisotropicMaterial
    talc_stiffness = [219.83e9,  59.66e9,  -4.82e9,  -0.82e9, -33.87e9, -1.04e9,
                      216.38e9, -3.67e9,   1.79e9, -16.51e9,  -0.62e9,
                      48.89e9,    4.12e9, -15.52e9,  -3.59e9,
                      26.54e9,    -3.6e9,  -6.41e9,
                      22.85e9,   -1.67e9,
                      78.29e9]
    rho = 2.75e3
    talc = anisotropy.AnisotropicMaterial(talc_stiffness, rho, 'triclinic')
    
    # Now we can calculate some isotropic properties of talc 
    # Print the bounds on bulk and shear modulus and
    # meaasures of elastic anisotropy and poisson ratio
    print('Talc elastic properties:\n')
    print('Bulk modulus bounds: {0:.3e} {1:.3e} {2:.3e}'.format(talc.bulk_modulus_reuss,
                                                     talc.bulk_modulus_vrh,
                                                     talc.bulk_modulus_voigt))
    print('Shear modulus bounds: {0:.3e} {1:.3e} {2:.3e}'.format(talc.shear_modulus_reuss,
                                                     talc.shear_modulus_vrh,
                                                     talc.shear_modulus_voigt))
    print('Universal elastic anisotropy: {0:.3f}\n'
          'Isotropic poisson ratio: {1:.3f}'.format(talc.universal_elastic_anisotropy,
                                                    talc.isotropic_poisson_ratio))

    # Let's make a pretty plot illustrating the anisotropy in talc
    talc.plot_velocities()


