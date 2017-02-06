# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_fit_data
----------------

This example demonstrates BurnMan's functionality to fit thermoelastic data to
both 2nd and 3rd orders using the EoS of the user's choice at 300 K. User's
must create a file with :math:`P, T` and :math:`V_s`. See input_minphys/ for example input
files.

requires:
- compute seismic velocities

teaches:
- averaging

"""
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman

if __name__ == "__main__":

    # First, read in the data from file and convert to SI units.
    PTp_data = np.loadtxt('../burnman/data/input_minphys/Murakami_perovskite.txt')
    PTp_data[:,0] = PTp_data[:,0]*1.e9
    PTp_data[:,2] = PTp_data[:,2]*1.e3

    # Make the test mineral
    mg_perovskite_test = burnman.Mineral()
    mg_perovskite_test.params = {'V_0': 24.45e-6,
                                 'K_0': 281.e9,
                                 'Kprime_0': 4.1,
                                 'molar_mass': .10,
                                 'G_0': 200.e9,
                                 'Gprime_0': 2.}


    def best_fit():
        return burnman.tools.fit_PTp_data(mineral = mg_perovskite_test,
                                          p_flags = 'shear_wave_velocity',
                                          fit_params = ['G_0', 'Gprime_0'],
                                          PTp = PTp_data,
                                          verbose = False)

    pressures = np.linspace(1.e5, 150.e9, 101)
    temperatures = pressures*0. + 300.

    # Fit to the second order Birch-Murnaghan EoS
    mg_perovskite_test.set_method("bm2")
    fitted_eos = best_fit()
    print('2nd order fit:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    model_vs_2nd_order_correct = mg_perovskite_test.evaluate(['shear_wave_velocity'],
                                                             pressures, temperatures)[0]
    mg_perovskite_test.set_method("bm3")
    model_vs_2nd_order_incorrect = mg_perovskite_test.evaluate(['shear_wave_velocity'],
                                                               pressures, temperatures)[0]
    print('')
    

    # Fit to the third order Birch-Murnaghan EoS

    mg_perovskite_test.set_method("bm3")
    fitted_eos = best_fit()
    print('3rd order fit:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    model_vs_3rd_order_correct = mg_perovskite_test.evaluate(['shear_wave_velocity'],
                                                             pressures, temperatures)[0]
    mg_perovskite_test.set_method("bm2")
    model_vs_3rd_order_incorrect = mg_perovskite_test.evaluate(['shear_wave_velocity'],
                                                               pressures, temperatures)[0]
    print('')


    plt.plot(pressures / 1.e9, model_vs_2nd_order_correct / 1000., color='r',
             linestyle='-', linewidth=2, label="Correct 2nd order fit")
    plt.plot(pressures / 1.e9, model_vs_2nd_order_incorrect / 1000., color='r',
             linestyle='-.', linewidth=2, label="Incorrect 2nd order fit")
    plt.plot(pressures / 1.e9, model_vs_3rd_order_correct / 1000., color='b',
             linestyle='-', linewidth=2, label="Correct 3rd order fit")
    plt.plot(pressures / 1.e9, model_vs_3rd_order_incorrect / 1000., color='b',
             linestyle='-.', linewidth=2, label="Incorrect 3rd order fit")
    plt.scatter(PTp_data[:,0] / 1.e9, PTp_data[:,2] / 1.e3)
    plt.ylim([6.55, 8])
    plt.xlim([25., 135.])
    plt.ylabel("Shear velocity (km/s)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="lower right", prop={'size': 12}, frameon=False)
    plt.savefig("output_figures/example_fit_data.png")
    plt.show()
