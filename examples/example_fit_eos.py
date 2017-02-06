# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_fit_eos
----------------

This example demonstrates BurnMan's functionality to fit PVT data to
an EoS of the user's choice. 

teaches:
- least squares fitting

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

if __name__ == "__main__":
    # First, let's read in the PVT equation of state data for MgO from Dewaele et al., (2000).
    T, Terr, Pta, P, Perr, V, Verr = np.loadtxt('../burnman/data/input_fitting/PVT_MgO_Dewaele_et_al_2000.dat', unpack=True)
    PTV = np.array([P*1.e9, T, burnman.tools.molar_volume_from_unit_cell_volume(V, 4.)]).T
    nul = 0.*PTV.T[0]

    Pcov = np.power(Perr*1.e9, 2.)
    Tcov = np.power(Terr, 2.)
    Vcov = np.power(burnman.tools.molar_volume_from_unit_cell_volume(Verr, 4.), 2.)
    
    PTV_covariance = np.array([[Pcov, nul, nul], [nul, Tcov, nul], [nul, nul, Vcov]]).T

    
    # Here's where we fit the data
    # We choose to fit to the SLB_2011 equation of state with four parameters (V_0, K_0, K'_0, grueneisen_0)
    per = burnman.minerals.SLB_2011.periclase()
    params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0']
    fitted_eos = burnman.tools.fit_PTV_data(per, params, PTV, PTV_covariance, verbose=False)

    
    # We're done! That wasn't too painful, was it?!

    
    # Print the optimized parameters
    print('Equation of state calculations')
    print('Optimized equation of state for periclase:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, params)
        
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, params)
    plt.show()

    # Now let's plot the residuals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    im = ax.scatter(PTV[:,0]/1.e9, PTV[:,1], c=fitted_eos.weighted_residuals, cmap=plt.cm.RdYlBu, s=50)
    fig.colorbar(im, ax=ax)
    max_abs_r = np.max(np.abs(fitted_eos.weighted_residuals))
    im.set_clim(-max_abs_r, max_abs_r)
    plt.show()
    
    # Here we plot our equation of state, along with the 95% confidence intervals for the volume
    pressures = np.linspace(1.e5, 60.e9, 101)
    volumes = np.empty_like(pressures)
    PTVs = np.empty((len(pressures), 3))
    
    for T in [298.15, 1000., 2000.]:
        for i, P in enumerate(pressures):
            per.set_state(P, T)
            PTVs[i] = [P, T, per.V]

        # Plot the 95% confidence bands 
        cp_bands = burnman.nonlinear_fitting.orthogonal_distance_confidence_prediction_bands(fitted_eos, PTVs, 0.95, [], 'V')
        plt.plot(cp_bands[0][:,0]/1.e9, cp_bands[0][:,2] * 1.e6, linestyle='--', color='r', label='95% confidence bands')
        plt.plot(cp_bands[1][:,0]/1.e9, cp_bands[1][:,2] * 1.e6, linestyle='--', color='r')
        
        plt.plot(PTVs[:,0] / 1.e9, PTVs[:,2] * 1.e6,
                 label='Optimized fit for periclase at {0:.0f} K'.format(T))
        
    plt.errorbar(PTV[:,0] / 1.e9, PTV[:,2] * 1.e6,
                 xerr=np.sqrt(PTV_covariance.T[0][0]) / 1.e9,
                 yerr=np.sqrt(PTV_covariance.T[2][2]) * 1.e6,
                 linestyle='None', marker='o', label='Dewaele et al. (2000)')

    plt.scatter(fitted_eos.data_mle[:,0] / 1.e9, fitted_eos.data_mle[:,2] * 1.e6, s=20, label='Maximum likelihood estimates')
    plt.ylabel("Volume (cm^3/mol)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Periclase EoS")
    plt.show()


    # We can also look at the uncertainty in other properties
    # For example, let's look at the uncertainty in P wave velocities, bulk modulus, thermal expansion and thermal pressure
    fig = plt.figure()
    for T in [298.15, 2000.]:
        for i, P in enumerate(pressures):
            per.set_state(P, T)
            PTVs[i] = [P, T, per.V]

        for i, (material_property, scaling, name) in enumerate([('p_wave_velocity', 1.e3, 'P wave velocity (km/s)'),
                                                                ('K_T', 1.e9, 'Bulk modulus (GPa)'),
                                                                ('alpha', 1., 'Thermal expansion (/K)'),
                                                                (['alpha', 'K_T'], 1.e6, 'Thermal pressure (MPa/K)')]):
            ax = fig.add_subplot(2, 2, i+1)
            
            # Plot the 95% confidence bands for the various material properties
            cp2_bands = burnman.nonlinear_fitting.confidence_prediction_bands(fitted_eos,
                                                                              burnman.tools.attribute_function(per, material_property),
                                                                              PTVs, 0.95, 'V')
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0] + cp2_bands[1])/2/scaling, label='Best fit at {0:.0f} K'.format(T))
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0])/scaling, linestyle='--', color='r', label='95% confidence band')
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[1])/scaling, linestyle='--', color='r')
            plt.ylabel(name)
            plt.xlabel('Pressure (GPa)')
            
    plt.legend(loc='upper right')
    plt.show()

    
    # Finally, let's look at the effect of uncertainties in the equation of state on gibbs free energy at high pressure
    fig = plt.figure()
    for T in [298.15, 2000.]:
        for i, P in enumerate(pressures):
            per.set_state(P, T)
            PTVs[i] = [P, T, per.V]

        ax = fig.add_subplot(1, 1, 1)
        scaling = 1.e3
        
        # Plot the 95% confidence bands for the gibbs free energy
        cp2_bands = burnman.nonlinear_fitting.confidence_prediction_bands(fitted_eos,
                                                                          burnman.tools.attribute_function(per, 'gibbs'),
                                                                          PTVs, 0.95, 'V')
        ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0] - cp2_bands[1])/2/scaling, label='95% confidence half width at {0:.0f} K'.format(T))
        plt.ylabel('Gibbs free energy uncertainty (kJ/mol)')
        plt.xlabel('Pressure (GPa)')
            
    plt.legend(loc='lower right')
    plt.show()

    
