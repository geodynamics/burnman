# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_compare_all_methods
---------------------------

This example demonstrates how to call each of the individual calculation
methodologies that exist within BurnMan. See below for current options. This
example calculates seismic velocity profiles for the same set of minerals and
a plot of :math:`V_s, V_\phi` and :math:`\\rho` is produce for the user to compare each of the
different methods.

*Specifically uses:*

* :doc:`eos`


*Demonstrates:*

* Each method for calculating velocity profiles currently included within BurnMan

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
from burnman import minerals

if __name__ == "__main__":
    # Input composition.

    amount_perovskite = 0.95
    rock = burnman.Composite([minerals.Murakami_etal_2012.fe_perovskite(),
                              minerals.Murakami_etal_2012.fe_periclase_LS()],
                             [amount_perovskite, 1.0 - amount_perovskite])

    #(min pressure, max pressure, pressure step)
    seis_p = np.arange(25e9, 125e9, 5e9)

    # Input adiabat potential temperature
    T0 = 1500.0

    # Now we'll calculate the models by forcing the rock to use a method. The
    # preset equation of state for the Murakami_etal_2012 minerals is 'slb2'

    """ 'slb2' (finite-strain 2nd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus,
        matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shear modulus,
        matas et al. 2007)
    or 'bm2' (birch-murnaghan 2nd order, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case))
    or 'bm3' (birch-murnaghan 3rd order, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case))"""

    methods = ['bm3', 'bm2', 'mgd3', 'mgd2', 'slb3', 'slb2']
    colors = ['r', 'k', 'g', 'b', 'y', 'm']
    markers = ['+', 'x', '>', '^', '<', 'v']

    plt.figure(figsize=(12, 10))

    for m in range(len(methods)):
        rock.set_method(methods[m])
        temperature = burnman.geotherm.adiabatic(seis_p, T0, rock)

        print("Calculations are done for:")
        rock.debug_print()

        mat_rho_1, mat_vs_1, mat_vphi_1 = \
            rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

        # Now let's plot the comparison. You can conversely just output to a data   file
        #(see example_woutput.py)
        # plot Vs
        plt.subplot(2, 2, 1)
        plt.plot(
            seis_p / 1.e9, mat_vs_1 / 1.e3, color=colors[m], linestyle='-', marker=markers[m],
            markerfacecolor=colors[m], markersize=4)

        # plot Vphi
        plt.subplot(2, 2, 2)
        plt.plot(
            seis_p / 1.e9, mat_vphi_1 / 1.e3, color=colors[m], linestyle='-', marker=markers[m],
            markerfacecolor=colors[m], markersize=4)

        # plot density
        plt.subplot(2, 2, 3)
        plt.plot(
            seis_p / 1.e9, mat_rho_1 / 1.e3, color=colors[m], linestyle='-', marker=markers[m],
            markerfacecolor=colors[m], markersize=4)

        # plot temperature
        plt.subplot(2, 2, 4)
        plt.plot(
            seis_p / 1.e9, temperature, color=colors[m], linestyle='-', marker=markers[m],
            markerfacecolor=colors[m], markersize=4, label=methods[m])
        plt.legend(loc='upper left')

    plt.subplot(2, 2, 1)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Vs (km/s)')
    plt.subplot(2, 2, 2)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Vphi (km/s)')
    plt.subplot(2, 2, 3)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Density (kg/m^3)')
    plt.subplot(2, 2, 4)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Temperature (K)')
    plt.savefig("output_figures/example_compare_all_methods.png")
    plt.show()
