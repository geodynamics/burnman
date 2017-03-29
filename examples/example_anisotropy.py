# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_anisotropy
-----------------

This example illustrates the basic functions required to convert 
an elastic stiffness tensor into elastic properties including:

See :cite:`Mainprice2011` Geological Society of London Special Publication 
and https://materialsproject.org/wiki/index.php/Elasticity_calculations
for mathematical descriptions of each function.

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
    # The next three lines initialise the AnisotropicMaterial
    talc_stiffness = [219.83e9, 59.66e9,  -4.82e9,  -0.82e9, -33.87e9, -1.04e9,
                      216.38e9, -3.67e9,   1.79e9, -16.51e9,  -0.62e9,
                      48.89e9,   4.12e9, -15.52e9,  -3.59e9,
                      26.54e9,   -3.6e9,  -6.41e9,
                      22.85e9,  -1.67e9,
                      78.29e9]
    rho = 2.75e3
    talc = anisotropy.AnisotropicMaterial(talc_stiffness, rho, 'triclinic')

    # Now we can calculate some isotropic properties of talc 
    # Print the bounds on bulk and shear modulus and
    # meaasures of elastic anisotropy and poisson ratio
    print('{0:.3e} {1:.3e} {2:.3e}'.format(talc.bulk_modulus_reuss,
                                           talc.bulk_modulus_vrh,
                                           talc.bulk_modulus_voigt))
    print('{0:.3e} {1:.3e} {2:.3e}'.format(talc.shear_modulus_reuss,
                                           talc.shear_modulus_vrh,
                                           talc.shear_modulus_voigt))
    print('{0:.3e} {1:.3e}'.format(talc.universal_elastic_anisotropy,
                                   talc.isotropic_poisson_ratio))

    # Let's make a pretty plot illustrating the anisotropy in talc
    talc.plot_velocities()

    # This plot uses the following functions:
    # talc.linear_compressibility(direction)
    # talc.youngs_modulus(direction)
    # talc.shear_modulus(plane_normal, shear_direction)
    # talc.poissons_ratio(longitudinal_direction, transverse_direction)
    # talc.wave_velocities(propagation_direction)
