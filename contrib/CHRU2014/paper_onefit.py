# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_onefit
------------

This script reproduces :cite:`Cottaar2014`, Figure 7.
It shows an example for a  best fit for a pyrolitic model within mineralogical error bars.
"""
from __future__ import absolute_import
from __future__ import print_function

import os.path
import sys
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import numpy.random
import burnman
import pickle
from burnman import minerals
from misc.helper_solid_solution import HelperSolidSolution
import matplotlib.cm
import matplotlib.colors
from scipy import interpolate
from scipy.stats import norm
import matplotlib.mlab as mlab
import misc.colors as colors


def make_rock():

    # approximate four component pyrolite model
    x_pv = 0.67
    x_fp = 0.33
    pv_fe_num = 0.07
    fp_fe_num = 0.2

    mg_perovskite = minerals.SLB_2011_ZSB_2013.mg_perovskite()
    fe_perovskite = minerals.SLB_2011_ZSB_2013.fe_perovskite()
    wuestite = minerals.SLB_2011_ZSB_2013.wuestite()
    periclase = minerals.SLB_2011_ZSB_2013.periclase()

    perovskite = HelperSolidSolution(
        [mg_perovskite, fe_perovskite], [1.0 - pv_fe_num, pv_fe_num])
    ferropericlase = HelperSolidSolution(
        [periclase, wuestite], [1.0 - fp_fe_num, fp_fe_num])

    pyrolite = burnman.Composite([perovskite, ferropericlase], [x_pv, x_fp])
    pyrolite.set_method('slb3')
    anchor_temperature = 1935.0

    return pyrolite, anchor_temperature


def output_rock(rock, file_handle):
    for ph in rock.staticphases:
        if(isinstance(ph.mineral, HelperSolidSolution)):
            for mineral in ph.mineral.endmembers:
                file_handle.write('\t' + mineral.to_string() + '\n')
                for key in mineral.params:
                    file_handle.write(
                        '\t\t' + key + ': ' + str(mineral.params[key]) + '\n')
        else:
            file_handle.write('\t' + ph.mineral.to_string() + '\n')
            for key in ph.mineral.params:
                file_handle.write(
                    '\t\t' + key + ': ' + str(ph.mineral.params[key]) + '\n')


def realization_to_array(rock, anchor_t):
    arr = [anchor_t]
    names = ['anchor_T']
    for ph in rock.staticphases:
        if(isinstance(ph.mineral, burnman.minerals_base.helper_solid_solution)):
            for mineral in ph.mineral.endmembers:
                for key in mineral.params:
                    if key != 'equation_of_state' and key != 'F_0' and key != 'T_0' and key != 'P_0':
                        arr.append(mineral.params[key])
                        names.append(mineral.to_string() + '.' + key)
        else:
            for key in ph.mineral.params:
                if key != 'equation_of_state' and key != 'F_0' and key != 'T_0' and key != 'P_0':
                    arr.append(ph.mineral.params[key])
                    names.append(mph.mineral.to_string() + '.' + key)
    return arr, names


def array_to_rock(arr, names):
    rock, _ = make_rock()
    anchor_t = arr[0]
    idx = 1
    for phase in rock.phases:
        if isinstance(phase, HelperSolidSolution):
            for mineral in phase.endmembers:
                while mineral.to_string() in names[idx]:
                    key = names[idx].split('.')[-1]
                    if key != 'equation_of_state' and key != 'F_0' and key != 'T_0' and key != 'P_0':
                        assert(mineral.to_string() in names[idx])
                        mineral.params[key] = arr[idx]
                        idx += 1
        else:
            raise Exception("unknown type")
    return rock, anchor_t


if __name__ == "__main__":
    # set up the seismic model
    seismic_model = burnman.seismic.PREM()
    npts = 10
    depths = np.linspace(850e3, 2700e3, npts)
    pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    pressures_sampled = np.linspace(
        pressure[0], pressure[-1], 20 * len(pressure))

    names = [
        'anchor_T', "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n",
        "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0", 'err', 'err_rho', 'err_vphi', 'err_vs']

    # those are the literature values:
    mymaplit = {'anchor_T': 2000,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0": 1.74,  # r:1.74
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0": 250.5e9,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0": 172.9e9,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0": 1.09,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0": 4.01,  # r:4.01
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0": 1.44,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0": 24.45e-6,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0": 9.059e2,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass": 0.1,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n": 5,
                "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0": 2.13,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0": 1.4,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0": 2.72e11,  # b: 2.637e11, r:2.72e11
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0": 1.33e11,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0": 1.1,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0": 4.1,  # r: 4.1
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0": 1.57,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0": 2.549e-05,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0": 8.71e2,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass": 0.1319,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n": 5,
                "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0": 2.3,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0": 2.1,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0": 161e9,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0": 1.310e11,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0": 1.700,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0": 3.8,  # b: 3.718 r:3.8
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0": 1.36,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0": 1.124e-05,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0": 767,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass": 0.0403,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n": 2,
                "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0": 2.8,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0": 1.4,  # r: 1.4
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0": 1.790e11,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0": 59.0e9,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0": 1.7,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0": 4.9,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0": 1.53,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0": 1.226e-05,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0": 4.54e2,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass": 0.0718,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n": 2,
                "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0": -0.1,
                'err': 0,
                'err_rho': 0,
                'err_vphi': 0,
                'err_vs': 0}

    # those are the values we use:
    mymap = {'anchor_T': 2000,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0": 1.779,  # r:1.74
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0": 250.5e9,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0": 172.9e9,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0": 1.09,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0": 3.917,  # r:4.01
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0": 1.44,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0": 24.45e-6,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0": 9.059e2,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass": 0.1,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n": 5,
             "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0": 2.13,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0": 1.4,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0": 2.637e11,  # b: 2.637e11, r:2.72e11
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0": 1.33e11,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0": 1.1,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0": 3.428,  # r: 4.1
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0": 1.57,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0": 2.549e-05,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0": 8.71e2,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass": 0.1319,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n": 5,
             "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0": 2.3,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0": 2.1,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0": 161e9,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0": 1.310e11,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0": 1.700,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0": 3.718,  # b: 3.718 r:3.8
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0": 1.36,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0": 1.124e-05,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0": 767,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass": 0.0403,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n": 2,
             "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0": 2.8,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0": 1.4,  # r: 1.4
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0": 1.790e11,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0": 59.0e9,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0": 1.7,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0": 4.9,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0": 1.53,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0": 1.226e-05,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0": 4.54e2,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass": 0.0718,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n": 2,
             "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0": -0.1,
             'err': 0,
             'err_rho': 0,
             'err_vphi': 0,
             'err_vs': 0}

    # make table:
    rows = ["V_0", "K_0", "Kprime_0", "G_0", "Gprime_0",
            "molar_mass", "n", "Debye_0", "grueneisen_0", "q_0", "eta_s_0"]
    for row in rows:
        val = []
        for n in names:
            if "." + row in n:
                val.append(mymap[n])
        print(row, "& %g && %g && %g && %g & \\" %
              (val[0], val[1], val[2], val[3]))

    dashstyle2 = (7, 3)
    dashstyle3 = (3, 2)

    fit = []
    lit = []
    for n in names:
        fit.append(mymap[n])
        lit.append(mymaplit[n])

    rock, anchor_t = array_to_rock(fit, names)
    temperature = burnman.geotherm.adiabatic(pressure, anchor_t, rock)
    rock.set_averaging_scheme(
        burnman.averaging_schemes.HashinShtrikmanAverage())
    rho, vp, vs, vphi, K, G = \
        rock.evaluate(
            ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], pressure, temperature)

    err_vs, err_vphi, err_rho = burnman.compare_l2(depths / np.mean(depths),
                                                   [vs / np.mean(seis_vs),
                                                    vphi / np.mean(seis_vphi),
                                                    rho / np.mean(seis_rho)],
                                                   [seis_vs / np.mean(seis_vs),
                                                    seis_vphi /
                                                    np.mean(seis_vphi),
                                                    seis_rho / np.mean(seis_rho)])
    error = np.sum([err_rho, err_vphi, err_vs])

    print("errors:", error, err_rho, err_vphi, err_vs)

    figsize = (6, 5)
    prop = {'size': 12}
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
    plt.rc('font', family='sans-serif')
    figure = plt.figure(dpi=100, figsize=figsize)

    # plot v_s
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, linestyle="-",
             color='k', linewidth=2.0, label='PREM')

    # plot v_phi
    plt.plot(pressure / 1.e9, seis_vphi / 1.e3,
             linestyle="-", color='k', linewidth=2.0)

    # plot density
    plt.plot(pressure / 1.e9, seis_rho / 1.e3,
             linestyle="-", color='k', linewidth=2.0)

    plt.plot(
        pressure / 1.e9, vphi / 1.e3, linestyle="-", color=colors.color(3),
        linewidth=1.0, marker='s', markerfacecolor=colors.color(3), label="vphi")
    plt.plot(pressure / 1.e9, vs / 1.e3, linestyle="-", color=colors.color(4),
             linewidth=1.0, marker='v', markerfacecolor=colors.color(4), label="vs")
    plt.plot(pressure / 1.e9, rho / 1.e3, linestyle="-", color=colors.color(2),
             linewidth=1.0, marker='o', markerfacecolor=colors.color(2), label="rho")

    rock, anchor_t = array_to_rock(lit, names)
    temperature = burnman.geotherm.adiabatic(pressure, anchor_t, rock)
    rock.set_averaging_scheme(
        burnman.averaging_schemes.HashinShtrikmanAverage())
    rho, vp, vs, vphi, K, G = \
        rock.evaluate(
            ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], pressure, temperature)
    plt.plot(pressure / 1.e9, vs / 1.e3, dashes=dashstyle2,
             color=colors.color(4), linewidth=1.0)
    plt.plot(pressure / 1.e9, vphi / 1.e3, dashes=dashstyle2,
             color=colors.color(3), linewidth=1.0, label="literature")
    plt.plot(pressure / 1.e9, rho / 1.e3, dashes=dashstyle2,
             color=colors.color(2), linewidth=1.0)

    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Velocities (km/s) and Density (kg/m$^3$)")
    plt.legend(bbox_to_anchor=(1.0, 0.9), prop={'size': 12})
    plt.xlim(25, 135)
    # plt.ylim(6,11)
    if "RUNNING_TESTS" not in globals():
        plt.savefig("onefit.pdf", bbox_inches='tight')
    print("wrote onefit.pdf")
    # plt.show()
