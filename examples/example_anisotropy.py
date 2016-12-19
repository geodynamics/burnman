# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_anisotropy
-----------------

This example illustrates the basic functions required to convert 
an elastic stiffness tensor into elastic properties including:

See :cite:`Mainprice2011` Geological Society of London Special Publication 
for a mathematical description of each function.

*Specifically uses:*

* :class:`burnman.anisotropy.voigt_notation_to_stiffness_tensor`
* :class:`christoffel_tensor`
* :class:`compliance_tensor`
* :class:`volume_compressibility`
* :class:`linear_compressibility`
* :class:`youngs_modulus`
* :class:`shear_modulus`
* :class:`poissons_ratio`
* :class:`wave_velocities`

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
    
    talc_stiffness = [[219.83,   59.66,   -4.82,  -0.82,  -33.87, -1.04],
                      [59.66,  216.38,   -3.67,   1.79,  -16.51,  -0.62],
                      [-4.82,   -3.67,   48.89,   4.12,  -15.52,  -3.59],
                      [-0.82,    1.79,    4.12,  26.54,    -3.6,  -6.41],
                      [-33.87,  -16.51,  -15.52,   -3.6,   22.85, -1.67],
                      [-1.04,   -0.62,   -3.59,  -6.41,   -1.67,  78.29]]
    '''
    K = 2.
    mu = 1.
    a = K + 4./3.*mu
    b = K - 2./3.*mu
    talc_stiffness = np.array([[a, b, b, 0., 0., 0.],
                               [b, a, b, 0., 0., 0.],
                               [b, b, a, 0., 0., 0.],
                               [0.,     0.,  0., mu, 0., 0.],
                               [0.,     0.,  0., 0., mu, 0.],
                               [0.,     0.,  0., 0., 0., mu]])
    '''
    #talc_stiffness = anisotropy.voigt_notation_to_stiffness_tensor(talc_stiffness)
    
    talc_compliance = anisotropy.compliance_tensor(talc_stiffness)
    beta = anisotropy.volume_compressibility(talc_stiffness)
    density = 3.

    zeniths = np.linspace(np.pi/2., np.pi, 31)
    azimuths = np.linspace(0., 2.*np.pi, 91)
    Rs = np.sin(zeniths)/(1. - np.cos(zeniths))
    r, theta = np.meshgrid(Rs, azimuths)
    
    vps = np.empty_like(r)
    vs1s = np.empty_like(r)
    vs2s = np.empty_like(r)
    betas = np.empty_like(r)
    Es = np.empty_like(r)
    for i, az in enumerate(azimuths):
        for j, phi in enumerate(zeniths):
            d = np.array([np.cos(az)*np.sin(phi), np.sin(az)*np.sin(phi), -np.cos(phi)]) # change_hemispheres
            velocities = anisotropy.wave_velocities(talc_stiffness, d, density)
            betas[i][j] = anisotropy.linear_compressibility(talc_stiffness, d)
            Es[i][j] = anisotropy.youngs_modulus(talc_stiffness, d)
            vps[i][j] = velocities[0][0]
            vs1s[i][j] = velocities[0][1]
            vs2s[i][j] = velocities[0][2]
            
            
            
    fig = plt.figure()

    names = ['Vp', 'anisotropy', 'Vp/Vs1', 'linear beta', 'Youngs Modulus']
    items = [vps, (vs1s - vs2s)/(vs1s + vs2s), vps/vs1s, betas, Es]
    ax = []
    im = []
    for i, item in enumerate(items):
        ax.append(fig.add_subplot(2, 3, i+1, projection='polar'))
        ax[i].set_title(names[i])
        im.append(ax[i].contourf(theta, r, item, 100, cmap=plt.cm.jet_r, vmin=np.min(item), vmax=np.max(item)))
        fig.colorbar(im[i], ax=ax[i])


    plt.show()
