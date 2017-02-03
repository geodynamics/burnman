# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_tools
----------------

This example demonstrates BurnMan's tools, which are currently
- EoS consistency
- equation of state fitting
- equilibrium temperature and pressure calculations
- Hugoniot calculation

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
    
    # First, let's check the EoS consistency of SLB_2011 periclase
    burnman.tools.check_eos_consistency(burnman.minerals.SLB_2011.periclase(), P=10.e9, T=3000., verbose=True)
    print('')
    
    # Next, let's create the Mg2SiO4 phase diagram
    forsterite = burnman.minerals.HP_2011_ds62.fo()
    mg_wadsleyite = burnman.minerals.HP_2011_ds62.mwd()
    mg_ringwoodite = burnman.minerals.HP_2011_ds62.mrw()
    periclase = burnman.minerals.HP_2011_ds62.per()
    mg_perovskite = burnman.minerals.HP_2011_ds62.mpv()

    temperatures = np.linspace(1000., 2000., 21)
    pressures_fo_wd = np.empty_like(temperatures)
    pressures_wd_rw = np.empty_like(temperatures)
    pressures_rw_perpv = np.empty_like(temperatures)

    # Here's one example where we find the equilibrium temperature:
    P = 14.e9
    T = burnman.tools.equilibrium_temperature([forsterite, mg_wadsleyite],
                                              [1.0, -1.0], P)
    print('Endmember equilibrium calculations')
    print('fo -> wad equilibrium at', P / 1.e9,
          "GPa is reached at", round_to_n(T, T, 4), "K")
    print('')

    # Now let's make the whole diagram using equilibrium_pressure
    for i, T in enumerate(temperatures):
        P = burnman.tools.equilibrium_pressure([forsterite, mg_wadsleyite],
                                               [1.0, -1.0],
                                               T)
        pressures_fo_wd[i] = P

        P = burnman.tools.equilibrium_pressure([mg_wadsleyite, mg_ringwoodite],
                                               [1.0, -1.0],
                                               T)
        pressures_wd_rw[i] = P

        P = burnman.tools.equilibrium_pressure(
            [mg_ringwoodite, periclase, mg_perovskite],
            [1.0, -1.0, -1.0],
            T)
        pressures_rw_perpv[i] = P

    plt.plot(temperatures, pressures_fo_wd / 1.e9, label='fo -> wd')
    plt.plot(temperatures, pressures_wd_rw / 1.e9, label='wd -> rw')
    plt.plot(temperatures, pressures_rw_perpv / 1.e9, label='rw -> per + pv')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pressure (GPa)")
    plt.legend(loc="lower right")
    plt.title("Mg2SiO4 phase diagram")
    plt.ylim(0., 30.)
    plt.show()

    # Now let's fit some EoS data
    # Let's just take a bit of data from Andrault et al. (2003) on stishovite
    PV = np.array([[0.0001, 46.5126, 0.0061],
                   [1.168, 46.3429, 0.0053],
                   [2.299, 46.1756, 0.0043],
                   [3.137, 46.0550, 0.0051],
                   [4.252, 45.8969, 0.0045],
                   [5.037, 45.7902, 0.0053],
                   [5.851, 45.6721, 0.0038],
                   [6.613, 45.5715, 0.0050],
                   [7.504, 45.4536, 0.0041],
                   [8.264, 45.3609, 0.0056],
                   [9.635, 45.1885, 0.0042],
                   [11.69, 44.947, 0.002],
                   [17.67, 44.264, 0.002],
                   [22.38, 43.776, 0.003],
                   [29.38, 43.073, 0.009],
                   [37.71, 42.278, 0.008],
                   [46.03, 41.544, 0.017],
                   [52.73, 40.999, 0.009],
                   [26.32, 43.164, 0.006],
                   [30.98, 42.772, 0.005],
                   [34.21, 42.407, 0.003],
                   [38.45, 42.093, 0.004],
                   [43.37, 41.610, 0.004],
                   [47.49, 41.280, 0.007]])
    
    PTV = np.array([PV.T[0]*1.e9,
                    PV.T[0]*0. + 298.15,
                    burnman.tools.molar_volume_from_unit_cell_volume(PV.T[1], 2.)]).T

    nul = 0.*PTV.T[0]
    PTV_covariances = np.array([[0.03*PTV.T[0], nul, nul], [nul, nul, nul], [nul, nul, burnman.tools.molar_volume_from_unit_cell_volume(PV.T[2], 2.)]]).T

    # Here's where we fit the data
    # The mineral parameters are automatically updated during fitting
    stv = burnman.minerals.HP_2011_ds62.stv()
    params = ['V_0', 'K_0', 'Kprime_0']
    fitted_eos = burnman.tools.fit_PTV_data(stv, params, PTV, PTV_covariances)

    # Print the optimized parameters
    print('Equation of state calculations')
    print('Optimized equation of state for stishovite:')
    for i, p in enumerate(params):
        print (p + ':', round_to_n(fitted_eos.popt[i], np.sqrt(fitted_eos.pcov[i][i]), 1),
               '+/-', round_to_n(np.sqrt(fitted_eos.pcov[i][i]), np.sqrt(fitted_eos.pcov[i][i]), 1))
    
    # Finally, let's plot our equation of state
    T = 298.15
    pressures = np.linspace(1.e5, 60.e9, 101)
    volumes = np.empty_like(pressures)

    PTVs = np.empty((len(pressures), 3))
    for i, P in enumerate(pressures):
        stv.set_state(P, T)
        PTVs[i] = [P, T, stv.V]
        
    cp_bands = fitted_eos.orthogonal_distance_confidence_prediction_bands(PTVs, 0.95, [])

    plt.plot(cp_bands[0][:,0]/1.e9, cp_bands[0][:,2] * 1.e6, linestyle='--', color='r', label='95% confidence bands')
    plt.plot(cp_bands[1][:,0]/1.e9, cp_bands[1][:,2] * 1.e6, linestyle='--', color='r')
    plt.plot(cp_bands[2][:,0]/1.e9, cp_bands[2][:,2] * 1.e6, linestyle='--', color='b', label='95% prediction bands')
    plt.plot(cp_bands[3][:,0]/1.e9, cp_bands[3][:,2] * 1.e6, linestyle='--', color='b')
        
    plt.plot(PTVs[:,0] / 1.e9, PTVs[:,2] * 1.e6,
             label='Optimized fit for stishovite')
    plt.errorbar(PTV[:,0] / 1.e9, PTV[:,2] * 1.e6,
                 xerr=PTV_covariances.T[0][0] / 1.e9,
                 yerr=PTV_covariances.T[2][2] * 1.e6,
                 linestyle='None', marker='o', label='Andrault et al. (2003)')

    plt.ylabel("Volume (cm^3/mol)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Stishovite EoS (room temperature)")
    plt.show()


    cp2_bands = fitted_eos.secondary_function_confidence_prediction_bands(PTVs, 0.95)
    plt.plot(PTVs[:,0]/1.e9, (cp2_bands[0] + cp2_bands[1])/2.e9, linestyle='--', color='b', label='Best fit')
    plt.plot(PTVs[:,0]/1.e9, (cp2_bands[0])/1.e9, linestyle='--', color='r', label='95% confidence band')
    plt.plot(PTVs[:,0]/1.e9, (cp2_bands[1])/1.e9, linestyle='--', color='r')
    plt.ylabel("Bulk modulus (GPa)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Stishovite EoS; uncertainty in bulk modulus (room temperature)")
    plt.show()

    
    # Here's a calculation of the Hugoniot of periclase up to 120 GPa
    print('')
    print('Hugoniot calculations')
    pressures = np.linspace(1.e5, 120.e9, 101)
    temperatures, volumes = burnman.tools.hugoniot(
        periclase, 1.e5, 298.15, pressures)
    plt.plot(pressures / 1.e9, temperatures, label='298.15 K')
    print('Room temperature Hugoniot temperature at',
          pressures[-1] / 1.e9, 'GPa:', int(temperatures[-1] + 0.5), 'K')

    temperatures, volumes = burnman.tools.hugoniot(
        periclase, 1.e5, 1000., pressures)
    plt.plot(pressures / 1.e9, temperatures, label='1000 K')
    print('1000 K Hugoniot temperature at',
          pressures[-1] / 1.e9, 'GPa:', int(temperatures[-1] + 0.5), 'K')

    plt.legend(loc="upper left")
    plt.ylabel("Temperature (K)")
    plt.xlabel("Pressure (GPa)")
    plt.title("Periclase Hugoniot")
    plt.show()
