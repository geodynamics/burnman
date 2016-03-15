from __future__ import absolute_import
from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import os.path
import sys
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))
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

import signal
import sys


def signal_handler(signal, frame):
    print('You pressed Ctrl+C!')
    sys.exit(0)
signal.signal(signal.SIGINT, signal_handler)


def normal(loc=0.0, scale=1.0):
    if scale <= 0.0:
        return 0.0
    else:
        return numpy.random.normal(loc, scale)


def realize_mineral(mineral):
    # some of the minerals are missing an uncertainty for this.  Assign a
    # characteristic nominal value for those
    if mineral.uncertainties['err_Gprime_0'] == 0.0:
        mineral.uncertainties['err_Gprime_0'] = 0.1

    # sample the uncertainties for all the relevant parameters.  Assume that
    # molar mass, V0, and n are well known
    mineral.params['K_0'] = mineral.params['K_0'] + \
        normal(scale=mineral.uncertainties['err_K_0'])
    mineral.params['Kprime_0'] = mineral.params['Kprime_0'] + \
        normal(scale=mineral.uncertainties['err_Kprime_0'])
    mineral.params['G_0'] = mineral.params['G_0'] + \
        normal(scale=mineral.uncertainties['err_G_0'])
    mineral.params['Gprime_0'] = mineral.params['Gprime_0'] + \
        normal(scale=mineral.uncertainties['err_Gprime_0'])
    mineral.params['Debye_0'] = mineral.params['Debye_0'] + \
        normal(scale=mineral.uncertainties['err_Debye_0'])
    mineral.params['grueneisen_0'] = mineral.params['grueneisen_0'] + \
        normal(scale=mineral.uncertainties['err_grueneisen_0'])
    mineral.params['q_0'] = mineral.params['q_0'] + \
        normal(scale=mineral.uncertainties['err_q_0'])
    mineral.params['eta_s_0'] = mineral.params['eta_s_0'] + \
        normal(scale=mineral.uncertainties['err_eta_s_0'])
    return mineral


def realize_pyrolite():

    # approximate four component pyrolite model
    x_pv = 0.67
    x_fp = 0.33
    pv_fe_num = 0.07
    fp_fe_num = 0.2

    mg_perovskite = minerals.SLB_2011_ZSB_2013.mg_perovskite()
    realize_mineral(mg_perovskite)
    fe_perovskite = minerals.SLB_2011_ZSB_2013.fe_perovskite()
    realize_mineral(fe_perovskite)
    wuestite = minerals.SLB_2011_ZSB_2013.wuestite()
    realize_mineral(wuestite)
    periclase = minerals.SLB_2011_ZSB_2013.periclase()
    realize_mineral(periclase)

    perovskite = HelperSolidSolution(
        [mg_perovskite, fe_perovskite], [1.0 - pv_fe_num, pv_fe_num])
    ferropericlase = HelperSolidSolution(
        [periclase, wuestite], [1.0 - fp_fe_num, fp_fe_num])

    pyrolite = burnman.Composite([perovskite, ferropericlase], [x_pv, x_fp])
    pyrolite.set_method('slb3')

    anchor_temperature = normal(loc=1935.0, scale=200.0)

    return pyrolite, anchor_temperature


def output_rock(rock, file_handle):
    for ph in rock.phases:
        if(isinstance(ph, HelperSolidSolution)):
            for min in ph.endmembers:
                file_handle.write('\t' + min.to_string() + '\n')
                for key in min.params:
                    file_handle.write(
                        '\t\t' + key + ': ' + str(min.params[key]) + '\n')
        else:
            file_handle.write('\t' + ph.to_string() + '\n')
            for key in ph.params:
                file_handle.write(
                    '\t\t' + key + ': ' + str(ph.params[key]) + '\n')


def realization_to_array(rock, anchor_t):
    arr = [anchor_t]
    names = ['anchor_T']
    for ph in rock.phases:
        if isinstance(ph, HelperSolidSolution):
            for min in ph.endmembers:
                for key in min.params:
                    if key != 'equation_of_state':
                        arr.append(min.params[key])
                        names.append(min.to_string() + '.' + key)
        else:
            for key in ph.params:
                if key != 'equation_of_state':
                    arr.append(ph.mineral.params[key])
                    names.append(ph.mineral.to_string() + '.' + key)
    return arr, names


def array_to_rock(arr, names):
    rock, _ = realize_pyrolite()
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


# set up the seismic model
seismic_model = burnman.seismic.PREM()
npts = 10
depths = np.linspace(850e3, 2700e3, npts)
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
    ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

n_realizations = 10000
min_error = np.inf

pressures_sampled = np.linspace(pressure[0], pressure[-1], 20 * len(pressure))
fname = 'output_pyrolite_uncertainty.txt'


whattodo = ""
dbname = ""

goodfits = []
names = []

if len(sys.argv) >= 2:
    whattodo = sys.argv[1]
if len(sys.argv) >= 3:
    dbname = sys.argv[2]

if whattodo == "plotgood" and len(sys.argv) > 2:
    files = sys.argv[2:]
    print("files:", files)
    names = pickle.load(open(files[0] + ".names", "rb"))
    erridx = names.index("err")
    print(erridx)

    allfits = []

    for f in files:
        a = pickle.load(open(f, "rb"))
        allfits.extend(a)
        b = a
        # b = [i for i in a if i[erridx]<3e-5] -- filter, need to adjust error
        # value
        print("adding %d out of %d" % (len(b), len(a)))
        goodfits.extend(b)

    minerr = min([f[erridx] for f in allfits])
    print("min error is %f" % minerr)

    num = len(goodfits)
    print("we have %d good entries" % num)

    i = 0
    idx = 0
    figsize = (20, 15)
    font = {'family': 'normal',
            'weight': 'normal',
            'size': 8}

    matplotlib.rc('font', **font)
    prop = {'size': 12}
    # plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
    plt.rc('font', family='sans-serif')
    figure = plt.figure(dpi=150, figsize=figsize)
    plt.subplots_adjust(hspace=0.3)
    for name in names:
        if name.endswith(".n") or name.endswith(".V_0") or name.endswith(".molar_mass") or name.endswith(".F_0") or name.endswith(".P_0"):
            i += 1
            continue
        plt.subplot(5, 8, idx)
        idx += 1
        shortname = name.replace("'burnman.minerals.SLB_2011_ZSB_2013", "").replace(
            "'", "").replace("perovskite", "p")

        trace = []
        for entry in allfits:
            trace.append(entry[i])
        # n, bins, patches = plt.hist(np.array(trace), 20, normed=1,
        # facecolor='blue', alpha=0.75)
        hist, bins = np.histogram(np.array(trace), bins=50, density=True)
        (mu, sigma) = norm.fit(np.array(trace))
        y = mlab.normpdf(bins, mu, sigma)
        if sigma > 1e-10 and not shortname.startswith("err"):
            l = plt.plot(bins, y, 'b--', linewidth=1)

        trace = []
        if shortname.startswith("err"):
            shortname += "(log)"
            for entry in goodfits:
                trace.append(np.log(entry[i]) / np.log(10))
        else:
            for entry in goodfits:
                trace.append(entry[i])
        hist, bins = np.histogram(trace)

        n, bins, patches = plt.hist(
            np.array(trace), 20, facecolor='green', alpha=0.75, normed=True)
        (mu, sigma) = norm.fit(np.array(trace))
        y = mlab.normpdf(bins, mu, sigma)
        if sigma > 1e-10 and not shortname.startswith("err"):
            l = plt.plot(bins, y, 'r--', linewidth=1)

        plt.title("%s\nmean %.3e sd: %.3e" %
                  (shortname, mu, sigma), fontsize=8)

        i += 1
    plt.savefig('good.png')
    print("Writing good.png")
    # plt.show()

    figsize = (8, 6)
    figure = plt.figure(dpi=150, figsize=figsize)

    for fit in goodfits:
        # print(fit)
        # print(names)

        rock, anchor_t = array_to_rock(fit, names)
        temperature = burnman.geotherm.adiabatic(pressure, anchor_t, rock)

        rock.set_averaging_scheme(
            burnman.averaging_schemes.HashinShtrikmanAverage())
        rho, vs, vphi = rock.evaluate(
            ['rho', 'v_s', 'v_phi'], pressure, temperature)

        print(".")

        plt.plot(pressure / 1.e9, vs / 1.e3,
                 linestyle="-", color='r', linewidth=1.0)
        plt.plot(pressure / 1.e9, vphi / 1.e3,
                 linestyle="-", color='b', linewidth=1.0)
        plt.plot(pressure / 1.e9, rho / 1.e3,
                 linestyle="-", color='g', linewidth=1.0)

    print("done!")

    # plot v_s
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, linestyle="--",
             color='k', linewidth=2.0, label='PREM')

    # plot v_phi
    plt.plot(pressure / 1.e9, seis_vphi / 1.e3,
             linestyle="--", color='k', linewidth=2.0, label='PREM')

    # plot density
    plt.plot(pressure / 1.e9, seis_rho / 1.e3, linestyle="--",
             color='k', linewidth=2.0, label='PREM')

    plt.savefig('goodones.png')
    print("Writing goodones.png")
    plt.show()

elif whattodo == "plotone":
    figsize = (6, 5)
    figure = plt.figure(dpi=100, figsize=figsize)

    names = [
        'anchor_T', "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n",
        "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0", 'err', 'err_rho', 'err_vphi', 'err_vs']

    mymapbestfitnotused = {'anchor_T': 2000,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0": 1.779,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0": 2.500e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0": 1.728e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0": 1.098,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0": 3.917,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0": 1.442,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0": 2.445e-05,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0": 9.057e2,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass": 0.1,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n": 5,
                           "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0": 2.104,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0": 1.401,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0": 2.637e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0": 1.329e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0": 1.084,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0": 3.428,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0": 1.568,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0": 2.549e-05,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0": 8.707e2,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass": 0.1319,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n": 5,
                           "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0": 2.335,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0": 2.108,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0": 1.610e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0": 1.310e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0": 1.700,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0": 3.718,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0": 1.359,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0": 1.124e-05,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0": 7.672,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass": 0.0403,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n": 2,
                           "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0": 2.804,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0": 1.405,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0": 1.790e11,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0": 5.905e10,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0": 1.681,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0": 4.884,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0": 1.532,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0": 1.226e-05,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0": 4.543e2,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass": 0.0718,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n": 2,
                           "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0": -7.048e-2,
                           'err': 0,
                           'err_rho': 0,
                           'err_vphi': 0,
                           'err_vs': 0}
    print(
        "goal: 5.35427067017e-06 2.72810809096e-07 3.67937164518e-06 1.4020882159e-06")

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
    rho, vs, vphi = rock.evaluate(
        ['rho', 'v_s', 'v_phi'], pressure, temperature)

    err_vs, err_vphi, err_rho = burnman.compare_l2(
        depths, [vs, vphi, rho], [seis_vs, seis_vphi, seis_rho])
    error = np.sum(
        [err_rho / np.mean(seis_rho), err_vphi / np.mean(seis_vphi), err_vs / np.mean(seis_vs)])

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
    rho, vs, vphi = rock.evaluate(
        ['rho', 'v_s', 'v_phi'], pressure, temperature)

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
    plt.savefig("onefit.pdf", bbox_inches='tight')
    print("wrote onefit.pdf")
    # plt.show()


elif whattodo == "run" and len(sys.argv) > 2:
    outfile = open(fname, 'w')
    outfile.write("#pressure\t Vs \t Vp \t rho \n")
    best_fit_file = open('output_pyrolite_closest_fit.txt', 'w')

    for i in range(n_realizations):
        if (i > 0 and i % 25 == 0):
        # save good fits
            print("saving %d fits to %s" % (len(goodfits), dbname))
            pickle.dump(goodfits, open(dbname + ".tmp", "wb"))
            os.rename(dbname + ".tmp", dbname)
            pickle.dump(names, open(dbname + ".names", "wb"))

        print("realization", i + 1)
        try:
            # create the ith model
            pyrolite, anchor_temperature = realize_pyrolite()
            temperature = burnman.geotherm.adiabatic(
                pressure, anchor_temperature, pyrolite)

            # calculate the seismic observables
            pyrolite.set_averaging_scheme(
                burnman.averaging_schemes.HashinShtrikmanAverage())
            rho, vs, vphi = pyrolite.evaluate(
                ['rho', 'v_s', 'v_phi'], pressure, temperature)

            # estimate the misfit with the seismic model
            err_vs, err_vphi, err_rho = burnman.compare_l2(
                depths, [vs, vphi, rho], [seis_vs, seis_vphi, seis_rho])

            error = np.sum(
                [err_rho / np.mean(seis_rho), err_vphi / np.mean(seis_vphi), err_vs / np.mean(seis_vs)])

            if error < min_error:
                min_error = error
                print("new best realization with error", error)
                best_fit_file.write('Current best fit : ' + str(error) + '\n')
                output_rock(pyrolite, best_fit_file)

            a, names = realization_to_array(pyrolite, anchor_temperature)
            a.extend([error, err_rho, err_vphi, err_vs])
            names.extend(["err", "err_rho", "err_vphi", "err_vs"])
            goodfits.append(a)

            # interpolate to a higher resolution line
            frho = interpolate.interp1d(pressure, rho)
            fs = interpolate.interp1d(pressure, vs)
            fphi = interpolate.interp1d(pressure, vphi)

            pressure_list = pressures_sampled
            density_list = frho(pressures_sampled)
            vs_list = fs(pressures_sampled)
            vphi_list = fphi(pressures_sampled)

            data = list(zip(pressure_list, vs_list, vphi_list, density_list))
            np.savetxt(outfile, data, fmt='%.10e', delimiter='\t')

        except ValueError:
            print("failed, skipping")

    outfile.close()
    best_fit_file.close()

elif whattodo == "error":
    values = [
        1957.020221991886, 1.6590112209181886, 249335164670.39246, 170883524675.03842, 0.8922515920546608, 4.083536182853109, 1.4680357687136616, 2.445e-05, 907.6618871363347, 0.1, 5, 1.4575168081960164, 1.3379195339709193, 260344929478.3809, 138077598973.27307, 0.17942226498091196, 1.3948903373340595, 1.436924855529012, 2.549e-05, 881.2532665499875, 0.1319, 5, 3.1204661890247394, 2.1411938868468483, 164407523972.7836,
        131594720803.07439, 1.855224221011796, 3.867545309505681, 1.2953203656315155, 1.124e-05, 769.8199298156555, 0.0403, 2, 2.8860489779521985, 1.4263617489128713, 177341125271.45096, 59131041052.46985, 2.352310980469468, 5.1279202520952545, 1.6021924873676925, 1.226e-05, 440.13042122457716, 0.0718, 2, -1.6065263588976038, 7.5954915681374134e-05, 9.6441602176002807e-07, 4.4326026287552629e-05, 3.0664473372061482e-05]
    names = [
        'anchor_T', "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.mg_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.n",
        "'burnman.minerals.SLB_2011_ZSB_2013.fe_perovskite'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.n", "'burnman.minerals.SLB_2011_ZSB_2013.periclase'.eta_s_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Gprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.K_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.G_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.q_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Kprime_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.grueneisen_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.V_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.Debye_0", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.molar_mass", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.n", "'burnman.minerals.SLB_2011_ZSB_2013.wuestite'.eta_s_0", 'err', 'err_rho', 'err_vphi', 'err_vs']

    rock, anchor_t = array_to_rock(values, names)
    temperature = burnman.geotherm.adiabatic(pressure, anchor_t, rock)

    rock.set_averaging_scheme(
        burnman.averaging_schemes.HashinShtrikmanAverage())
    rho, vs, vphi = rock.evaluate(
        ['rho', 'v_s', 'v_phi'], pressure, temperature)

    err_vs, err_vphi, err_rho = burnman.compare_l2(
        depths, [vs, vphi, rho], [seis_vs, seis_vphi, seis_rho])
    error = np.sum([err_rho, err_vphi, err_vs])

    print(error, err_rho, err_vphi, err_vs)


elif whattodo == "plot":
    infile = open(fname, 'r')
    data = np.loadtxt(fname, skiprows=1)
    pressure_list = data[:, 0]
    density_list = data[:, 3]
    vs_list = data[:, 1]
    vphi_list = data[:, 2]
    infile.close()

    density_hist, rho_xedge, rho_yedge = np.histogram2d(
        pressure_list, density_list, bins=len(pressures_sampled), normed=True)
    vs_hist, vs_xedge, vs_yedge = np.histogram2d(
        pressure_list, vs_list, bins=len(pressures_sampled), normed=True)
    vphi_hist, vphi_xedge, vphi_yedge = np.histogram2d(
        pressure_list, vphi_list, bins=len(pressures_sampled), normed=True)

    vs_xedge /= 1.e9
    vphi_xedge /= 1.e9
    rho_xedge /= 1.e9
    vs_yedge /= 1.e3
    vphi_yedge /= 1.e3
    rho_yedge /= 1.e3

    left_edge = min(vs_xedge[0], vphi_xedge[0], rho_xedge[0])
    right_edge = max(vs_xedge[-1], vphi_xedge[-1], rho_xedge[-1])
    bottom_edge = 4.3
    top_edge = 11.3
    aspect_ratio = (right_edge - left_edge) / (top_edge - bottom_edge)
    gamma = 0.8  # Mess with this to change intensity of colormaps near the edges

    # do some setup for the figure
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
    plt.rc('font', family='sans-serif')
    plt.subplots_adjust(wspace=0.3)

    plt.subplot(111, aspect='equal')
    plt.xlim(left_edge, right_edge)
    plt.ylim(bottom_edge, top_edge)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Wave Speed (km/s)')

    # plot v_s
    vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
    c = matplotlib.colors.LinearSegmentedColormap.from_list(
        'vphi', [(0, '#ffffff'), (0.2, '#eff3ff'), (0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c')], gamma=gamma)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, linestyle="--",
             color='k', linewidth=2.0, label='PREM')

    # plot v_phi
    vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
    c = matplotlib.colors.LinearSegmentedColormap.from_list(
        'vphi', [(0, '#ffffff'), (0.2, '#fee5d9'), (0.4, '#fcae91'), (0.6, '#fb6a4a'), (0.8, '#de2d26'), (1.0, '#a50f15')], gamma=gamma)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_vphi / 1.e3,
             linestyle="--", color='k', linewidth=2.0, label='PREM')

    # plot density
    density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
    c = matplotlib.colors.LinearSegmentedColormap.from_list(
        'vphi', [(0, '#ffffff'), (0.2, '#edf8e9'), (0.4, '#bae4b3'), (0.6, '#74c476'), (0.8, '#31a354'), (1.0, '#006d2c')], gamma=gamma)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        density_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_rho / 1.e3, linestyle="--",
             color='k', linewidth=2.0, label='PREM')

    # save and show the image
    fig = plt.gcf()
    fig.set_size_inches(6.0, 6.0)
    if "RUNNING_TESTS" not in globals():
        fig.savefig("pyrolite_uncertainty.pdf", bbox_inches='tight', dpi=100)
        print("Writing pyrolite_uncertainty.pdf")
    plt.show()

else:
    print("Options:")
    print(
        "  run <dbname>    -- run realizations and write into given database name")
    print("  plot <dbname>   -- plot given database")
    print(
        "  plotgood <dbname1> <dbname2> ...   -- aggregate databases and plot")
    print("  plotone  -- plot a single hardcoded nice result")
    print(
        "  error    -- testing, compute errors of a single hardcoded realization")
