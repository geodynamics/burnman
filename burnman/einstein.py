# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, 2014, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import constants

"""
Functions for the Einstein model of a solid.
"""


def thermal_energy(T, einstein_T, n):
    """
    calculate the thermal energy of a substance.  Takes the temperature,
    the Einstein temperature, and n, the number of atoms per molecule.
    Returns thermal energy in J/mol
    """
    if T == 0:
        return 3.*n*constants.gas_constant*einstein_T*0.5 # zero point energy
    x = einstein_T/T
    E_th = 3.*n*constants.gas_constant*einstein_T*( 0.5 + 1. / (np.exp( x ) - 1.0) ) # include the zero point energy
    return E_th

def heat_capacity_v(T,einstein_T,n):
    """
    Heat capacity at constant volume.  In J/K/mol
    """
    if T ==0:
        return 0
    x = einstein_T/T
    C_v = 3.0*n*constants.gas_constant* ( x * x * np.exp( x ) / np.power( np.exp( x ) - 1.0, 2.0 ) )
    return C_v




if __name__ == "__main__":

    import matplotlib.pyplot as plt

    temperatures = np.linspace(0,5000, 200)
    vibrational_energy = np.empty_like(temperatures)
    heat_capacity = np.empty_like(temperatures)
    Einstein_T = 1000.
    for i in range(len(temperatures)):
      vibrational_energy[i] = thermal_energy(temperatures[i], Einstein_T, 1.0)
      heat_capacity[i] = heat_capacity_v(temperatures[i], Einstein_T, 1.0)

    plt.subplot(121)
    plt.plot(temperatures, vibrational_energy)
    plt.subplot(122)
    plt.plot(temperatures, heat_capacity)
    plt.show()



