# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
EoS calculation example
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

def round_to_n(x, xerr, n):
    return round(x, -int(np.floor(np.log10(np.abs(xerr)))) + (n - 1))

if __name__ == "__main__":
    # First, let's read in the PVT equation of state data for MgO from Dewaele et al., (2000).
    T, Terr, Pta, P, Perr, V, Verr = np.loadtxt('../burnman/data/input_fitting/PVT_MgO_Dewaele_et_al_2000.dat', unpack=True)
    PTV = np.array([P*1.e9, T, burnman.tools.molar_volume_from_unit_cell_volume(V, 4.)]).T
    nul = 0.*PTV.T[0]
    PTV_covariance = np.array([[Perr*1.e9, nul, nul], [nul, Terr, nul], [nul, nul, burnman.tools.molar_volume_from_unit_cell_volume(Verr, 4.)]]).T

    
    # Here's where we fit the data
    # We choose to fit to the SLB_2011 equation of state with four parameters (V_0, K_0, K'_0, grueneisen_0)
    per = burnman.minerals.SLB_2011.periclase()
    params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0']
    fitted_eos = burnman.tools.fit_PTV_data(per, params, PTV, PTV_covariance, verbose=False)

    
    # We're done! That wasn't too painful, was it?!

    
    # Print the optimized parameters
    print('Equation of state calculations')
    print('Optimized equation of state for periclase:')
    for i, p in enumerate(params):
        p_rnd = round_to_n(fitted_eos.popt[i], np.sqrt(fitted_eos.pcov[i][i]), 1)
        c_rnd = round_to_n(np.sqrt(fitted_eos.pcov[i][i]), np.sqrt(fitted_eos.pcov[i][i]), 1)
        scale = np.power(10., np.floor(np.log10(p_rnd)))
        nd = np.floor(np.log10(p_rnd)) - np.floor(np.log10(c_rnd))
        print ('{0:s}: ({1:{4}{5}f} +/- {2:{4}{5}f}) x {3:.0e}'.format(p, p_rnd/scale, c_rnd/scale, scale, 0, (nd)/10.))

        
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, params)
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
        cp_bands = burnman.nonlinear_fitting.orthogonal_distance_confidence_prediction_bands(fitted_eos, PTVs, 0.95, [])
        plt.plot(cp_bands[0][:,0]/1.e9, cp_bands[0][:,2] * 1.e6, linestyle='--', color='r', label='95% confidence bands')
        plt.plot(cp_bands[1][:,0]/1.e9, cp_bands[1][:,2] * 1.e6, linestyle='--', color='r')
        
        plt.plot(PTVs[:,0] / 1.e9, PTVs[:,2] * 1.e6,
                 label='Optimized fit for periclase at {0:.0f} K'.format(T))
        
    plt.errorbar(PTV[:,0] / 1.e9, PTV[:,2] * 1.e6,
                 xerr=PTV_covariance.T[0][0] / 1.e9,
                 yerr=PTV_covariance.T[2][2] * 1.e6,
                 linestyle='None', marker='o', label='Dewaele et al. (2000)')

    plt.ylabel("Volume (cm^3/mol)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Periclase EoS")
    plt.show()


    # We can also look at the uncertainty in other properties
    # For example, let's look at the uncertainty in P wave velocities
    for T in [298.15, 2000.]:
        for i, P in enumerate(pressures):
            per.set_state(P, T)
            PTVs[i] = [P, T, per.V]
            
        # Plot the 95% confidence bands for the bulk modulus
        cp2_bands = burnman.nonlinear_fitting.confidence_prediction_bands(fitted_eos,
                                                                          burnman.tools.attribute_function(per, 'p_wave_velocity'),
                                                                          PTVs, 0.95)
        plt.plot(PTVs[:,0]/1.e9, (cp2_bands[0] + cp2_bands[1])/2.e3, label='Best fit at {0:.0f} K'.format(T))
        plt.plot(PTVs[:,0]/1.e9, (cp2_bands[0])/1.e3, linestyle='--', color='r', label='95% confidence band')
        plt.plot(PTVs[:,0]/1.e9, (cp2_bands[1])/1.e3, linestyle='--', color='r')
    plt.ylabel("P wave velocity (km/s)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="lower right")
    plt.title("Periclase EoS; uncertainty in P wave velocities")
    plt.show()

    
