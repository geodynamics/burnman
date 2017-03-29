# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_tools
----------------

This example demonstrates BurnMan's tools, which are currently
- EoS consistency
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
    
    # Here's a calculation of the Hugoniot of periclase up to 120 GPa
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
