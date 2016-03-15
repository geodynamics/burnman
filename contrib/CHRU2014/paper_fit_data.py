# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_fit_data
--------------

This script reproduces :cite:`Cottaar2014` Figure 4.

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
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import scipy.optimize as opt
import burnman
import misc.colors as colors

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


figsize = (6, 5)
prop = {'size': 12}
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
plt.rc('font', family='sans-serif')
figure = plt.figure(dpi=100, figsize=figsize)


def calc_shear_velocities(G_0, Gprime_0, mineral, pressures):

    mineral.params['G_0'] = G_0
    mineral.params['Gprime_0'] = Gprime_0

    shear_velocities = np.empty_like(pressures)
    for i in range(len(pressures)):
        mineral.set_state(pressures[i], 0.0)
                          # set state with dummy temperature
        shear_velocities[i] = mineral.v_s

    return shear_velocities


def error(guess, test_mineral, pressures, obs_vs):
    vs = calc_shear_velocities(guess[0], guess[1], test_mineral, pressures)

    vs_l2 = [(vs[i] - obs_vs[i]) * (vs[i] - obs_vs[i])
             for i in range(len(obs_vs))]
    l2_error = sum(vs_l2)

    return l2_error


if __name__ == "__main__":

    mg_perovskite_data = np.loadtxt("Murakami_perovskite.txt")
    obs_pressures = mg_perovskite_data[:, 0] * 1.e9
    obs_vs = mg_perovskite_data[:, 2] * 1000.

    pressures = np.linspace(25.e9, 135.e9, 100)

    # make the mineral to fit
    guess = [200.e9, 2.0]
    mg_perovskite_test = burnman.Mineral()
    mg_perovskite_test.params['V_0'] = 24.45e-6
    mg_perovskite_test.params['K_0'] = 281.e9
    mg_perovskite_test.params['Kprime_0'] = 4.1
    mg_perovskite_test.params['molar_mass'] = .10227

    # first, do the second-order fit
    mg_perovskite_test.set_method("bm2")
    func = lambda x: error(x, mg_perovskite_test, obs_pressures, obs_vs)
    sol = opt.fmin(func, guess)
    print("2nd order fit: G = ", sol[0] / 1.e9, "GPa\tG' = ", sol[1])
    model_vs_2nd_order_correct = calc_shear_velocities(
        sol[0], sol[1], mg_perovskite_test, pressures)
    mg_perovskite_test.set_method("bm3")
    model_vs_2nd_order_incorrect = calc_shear_velocities(
        sol[0], sol[1], mg_perovskite_test, pressures)

    # now do third-order fit
    mg_perovskite_test.set_method("bm3")
    func = lambda x: error(x, mg_perovskite_test, obs_pressures, obs_vs)
    sol = opt.fmin(func, guess)
    print("3rd order fit: G = ", sol[0] / 1.e9, "GPa\tG' = ", sol[1])
    model_vs_3rd_order_correct = calc_shear_velocities(
        sol[0], sol[1], mg_perovskite_test, pressures)
    mg_perovskite_test.set_method("bm2")
    model_vs_3rd_order_incorrect = calc_shear_velocities(
        sol[0], sol[1], mg_perovskite_test, pressures)

    plt.plot(
        pressures / 1.e9, model_vs_2nd_order_correct / 1000., color=colors.color(3), linestyle='-',
        marker='x', markevery=7, linewidth=1.5, label="Correct 2nd order extrapolation")
    plt.plot(
        pressures / 1.e9, model_vs_2nd_order_incorrect / 1000., color=colors.color(3), linestyle='--',
        marker='x', markevery=7, linewidth=1.5, label="2nd order fit, 3rd order extrapolation")
    plt.plot(
        pressures / 1.e9, model_vs_3rd_order_correct / 1000., color=colors.color(1),
        linestyle='-', linewidth=1.5, label="Correct 3rd order extrapolation")
    plt.plot(
        pressures / 1.e9, model_vs_3rd_order_incorrect / 1000., color=colors.color(1),
        linestyle='--', linewidth=1.5, label="3rd order fit, 2nd order extrapolation")
    plt.scatter(obs_pressures / 1.e9, obs_vs /
                1000., zorder=1000, marker='o', c='w')
    plt.ylim([6.7, 8])
    plt.xlim([25., 135.])
    if "RUNNING_TESTS" not in globals():
        plt.ylabel(
            r'Shear velocity ${V}_{\mathlarger{\mathlarger{\mathlarger{s}}}}$ (km/s)')
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="lower right", prop=prop)
    if "RUNNING_TESTS" not in globals():
        plt.savefig("example_fit_data.pdf", bbox_inches='tight')
    plt.show()
