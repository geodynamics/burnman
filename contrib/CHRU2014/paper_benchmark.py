# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_benchmark
---------------

This script reproduces the benchmark in :cite:`Cottaar2014`, Figure 3.
"""
from __future__ import absolute_import


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))
sys.path.insert(1, os.path.abspath('.'))
import burnman

figsize = (6, 5)
prop = {'size': 12}
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
figure = plt.figure(dpi=100, figsize=figsize)


def check_slb_fig7_txt():
    """
    Calculates all values for forsterite and benchmarks with values from Stixrude and Lithgow-Bertelloni (personal communication)
    """
    forsterite = burnman.Mineral()
    forsterite.params = {'name': 'forsterite',
                         'V_0': 43.603e-6,
                         'K_0': 127.955e9,
                         'Kprime_0': 4.232,
                         'G_0': 81.6e9,
                         'Gprime_0': 1.4,
                         'molar_mass': .140695,
                         'n': 7.0,
                         'Debye_0': 809.183,
                         'grueneisen_0': .993,
                         'q_0': 2.093,
                         'eta_s_0': 2.364}
    forsterite.set_method('slb3')

    data = np.loadtxt("slb_benchmark.txt", skiprows=1)

    temperature = np.array(data[:, 2])
    pressure = np.array(data[:, 0])
    rho = np.array(data[:, 3])
    frho = np.empty_like(rho)
    rho_comp = np.empty_like(rho)
    Kt = np.array(data[:, 4])
    fKt = np.empty_like(Kt)
    Kt_comp = np.empty_like(Kt)
    Ks = np.array(data[:, 5])
    fKs = np.empty_like(Ks)
    Ks_comp = np.empty_like(Ks)
    G = np.array(data[:, 6])
    fG = np.empty_like(G)
    G_comp = np.empty_like(G)
    VB = np.array(data[:, 7])
    fVB = np.empty_like(VB)
    VB_comp = np.empty_like(VB)
    VS = np.array(data[:, 8])
    fVS = np.empty_like(VS)
    VS_comp = np.empty_like(VS)
    VP = np.array(data[:, 9])
    fVP = np.empty_like(VP)
    VP_comp = np.empty_like(VP)
    vol = np.array(data[:, 10])
    fvol = np.empty_like(vol)
    vol_comp = np.empty_like(vol)
    alpha = np.array(data[:, 11])
    falpha = np.empty_like(alpha)
    alpha_comp = np.empty_like(alpha)
    Cp = np.array(data[:, 12])
    fCp = np.empty_like(Cp)
    Cp_comp = np.empty_like(Cp)
    gr = np.array(data[:, 13])
    gr_comp = np.empty_like(gr)

    for i in range(len(temperature)):
        forsterite.set_state(pressure[i], temperature[i])
        rho_comp[i] = 100. * (forsterite.density / 1000. - rho[i]) / rho[i]
        Kt_comp[i] = 100. * (
            forsterite.isothermal_bulk_modulus / 1.e9 - Kt[i]) / Kt[i]
        Ks_comp[i] = 100. * (
            forsterite.adiabatic_bulk_modulus / 1.e9 - Ks[i]) / Ks[i]
        G_comp[i] = 100. * (forsterite.shear_modulus / 1.e9 - G[i]) / G[i]
        VB_comp[i] = 100. * (forsterite.v_phi / 1000. - VB[i]) / VB[i]
        VS_comp[i] = 100. * (forsterite.v_s / 1000. - VS[i]) / VS[i]
        VP_comp[i] = 100. * (forsterite.v_p / 1000. - VP[i]) / VP[i]
        vol_comp[i] = 100. * (forsterite.molar_volume * 1.e6 - vol[i]) / vol[i]
        alpha_comp[i] = 100. * (
            forsterite.thermal_expansivity / 1.e-5 - alpha[i]) / (alpha[-1])
        Cp_comp[i] = 100. * (forsterite.heat_capacity_p /
                             forsterite.params['molar_mass'] / 1000. - Cp[i]) / (Cp[-1])
        gr_comp[i] = (forsterite.grueneisen_parameter - gr[i]) / gr[i]

    plt.plot(temperature, rho_comp, label=r'$\rho$')
    plt.plot(temperature, Kt_comp, label=r'$K_S$')
    plt.plot(temperature, Ks_comp, label=r'$K_T$')
    plt.plot(temperature, G_comp, label=r'$G$')
    plt.plot(temperature, VS_comp, label=r'$V_S$')
    plt.plot(temperature, VP_comp, label=r'$V_P$')
    plt.plot(temperature, VB_comp, label=r'$V_\phi$')
    plt.plot(temperature, vol_comp, label=r'$V$')
    plt.plot(temperature, alpha_comp, label=r'$\alpha$')
    plt.plot(temperature, Cp_comp, label=r'$c_P$')
    plt.plot(temperature, gr_comp, label=r'$\gamma$')

    plt.xlim([0, 2200])
    plt.ylim([-0.002, 0.002])
    plt.yticks([-0.002, -0.001, 0, 0.001, 0.002])
    plt.xticks([0, 800, 1600, 2200])
    plt.xlabel("Temperature (K)")
    plt.ylabel("Difference (\%)")
    plt.legend(loc="lower center", prop=prop, ncol=4)
    if "RUNNING_TESTS" not in globals():
        plt.savefig("benchmark1.pdf", bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    check_slb_fig7_txt()
