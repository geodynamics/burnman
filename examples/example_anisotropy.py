# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_anisotropy
------------------

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

    try:
        plt.style.use('ggplot')
        plt.rcParams['axes.facecolor'] = 'white'
        plt.rcParams['axes.edgecolor'] = 'black'
        plt.rcParams['figure.figsize'] = 16, 10 # inches
    except:
        pass

    
    # Let's first look at an isotropic material
    # As an example, let's use the lambda (C12), mu (C44) and rho
    # for flow basalt given by Christensen et al. (1980)
    # Initial Reports of the Deep Sea Drilling Project 59: 515-17.

    elastic_constants = [0.4e11, 0.24e11]
    rho = 2.735e3
    basalt = anisotropy.IsotropicMaterial(rho, elastic_constants)

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
    
    
    # Now let's look at the properties of an anisotropic material.
    # Here we choose talc, as it is the mineral used as an example
    # in Mainprice et al. (2011)
    talc_stiffness = [219.83e9,  59.66e9,  -4.82e9,  -0.82e9, -33.87e9, -1.04e9,
                      216.38e9, -3.67e9,   1.79e9, -16.51e9,  -0.62e9,
                      48.89e9,    4.12e9, -15.52e9,  -3.59e9,
                      26.54e9,    -3.6e9,  -6.41e9,
                      22.85e9,   -1.67e9,
                      78.29e9]
    rho = 2.75e3
    talc = anisotropy.TriclinicMaterial(rho, talc_stiffness)

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

    # Finally, let's make a pretty plot illustrating the anisotropy in talc
    def plot_anisotropic_seismic_properties(mineral):
        """
        Makes colour plots of:
        Compressional wave velocity: Vp
        Anisotropy: (Vs1 - Vs2)/(Vs1 + Vs2)
        Vp/Vs1
        linear compressibility: beta
        Youngs Modulus: E
        """
        
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
                velocities = mineral.wave_velocities(d)
                betas[i][j] = mineral.linear_compressibility(d)
                Es[i][j] = mineral.youngs_modulus(d)
                vps[i][j] = velocities[0][0]
                vs1s[i][j] = velocities[0][1]
                vs2s[i][j] = velocities[0][2]
                
        fig = plt.figure()
        names = ['Vp (km/s)', 'Vs1 (km/s)', 'Vp/Vs1', 'S-wave anisotropy (%)', 'Linear compressibility (GPa$^{-1}$)', 'Youngs Modulus (GPa)']
        items = [vps/1000., vs1s/1000., vps/vs1s, 200.*(vs1s - vs2s)/(vs1s + vs2s), betas*1.e9, Es/1.e9]
        ax = []
        im = []
        ndivs = 100
        for i, item in enumerate(items):
            ax.append(fig.add_subplot(2, 3, i+1, projection='polar'))
            ax[i].set_yticks([100])
            ax[i].set_title(names[i])

            vmin = np.min(item)
            vmax = np.max(item)
            spacing = np.power(10., np.floor(np.log10(vmax - vmin)))
            nt = int((vmax - vmin - vmax%spacing + vmin%spacing)/spacing)
            if nt == 1:
                spacing = spacing/4.
            elif nt < 4:
                spacing = spacing/2.
            elif nt > 8:
                spacing = spacing*2.
                
            tmin = vmin + (spacing - vmin%spacing)
            tmax = vmax - vmax%spacing
            nt = int((tmax - tmin)/spacing + 1)
            
            ticks = np.linspace(tmin, tmax, nt)
            im.append(ax[i].contourf(theta, r, item, ndivs, cmap=plt.cm.jet_r, vmin=vmin, vmax=vmax))
            lines = ax[i].contour(theta, r, item, ticks, colors=('black',), linewidths=(1,))
            
            cbar = fig.colorbar(im[i], ax=ax[i], ticks=ticks)
            cbar.add_lines(lines)

        plt.tight_layout()
        plt.savefig("output_figures/example_anisotropy.png")
        plt.show()

    plot_anisotropic_seismic_properties(talc)


