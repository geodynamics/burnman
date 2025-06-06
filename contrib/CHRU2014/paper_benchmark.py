# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_benchmark
---------------

This script reproduces the benchmark in :cite:`Cottaar2014`, Figure 3.
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

if not os.path.exists("burnman") and os.path.exists("../../burnman"):
    sys.path.insert(1, os.path.abspath("../.."))
sys.path.insert(1, os.path.abspath("."))
import burnman

figsize = (6, 5)
prop = {"size": 12}
plt.rc("text", usetex=True)
plt.rc("font", family="sans-serif")
figure = plt.figure(dpi=100, figsize=figsize)


def check_slb_fig7_txt():
    """
    Calculates all values for forsterite and benchmarks with values from Stixrude and Lithgow-Bertelloni (personal communication)
    """
    forsterite = burnman.Mineral()
    forsterite.params = {
        "name": "forsterite",
        "V_0": 43.603e-6,
        "K_0": 127.955e9,
        "Kprime_0": 4.232,
        "G_0": 81.6e9,
        "Gprime_0": 1.4,
        "molar_mass": 0.140695,
        "n": 7.0,
        "Debye_0": 809.183,
        "grueneisen_0": 0.993,
        "q_0": 2.093,
        "eta_s_0": 2.364,
    }
    forsterite.set_method("slb3")

    data = np.loadtxt("slb_benchmark.txt", skiprows=1)

    temperature = np.array(data[:, 2])
    pressure = np.array(data[:, 0])
    rho = np.array(data[:, 3])
    rho_comp = np.empty_like(rho)
    Kt = np.array(data[:, 4])
    Kt_comp = np.empty_like(Kt)
    Ks = np.array(data[:, 5])
    Ks_comp = np.empty_like(Ks)
    G = np.array(data[:, 6])
    G_comp = np.empty_like(G)
    VB = np.array(data[:, 7])
    VB_comp = np.empty_like(VB)
    VS = np.array(data[:, 8])
    VS_comp = np.empty_like(VS)
    VP = np.array(data[:, 9])
    VP_comp = np.empty_like(VP)
    vol = np.array(data[:, 10])
    vol_comp = np.empty_like(vol)
    alpha = np.array(data[:, 11])
    alpha_comp = np.empty_like(alpha)
    Cp = np.array(data[:, 12])
    Cp_comp = np.empty_like(Cp)
    gr = np.array(data[:, 13])
    gr_comp = np.empty_like(gr)

    for i in range(len(temperature)):
        forsterite.set_state(pressure[i], temperature[i])
        rho_comp[i] = 100.0 * (forsterite.density / 1000.0 - rho[i]) / rho[i]
        Kt_comp[i] = (
            100.0 * (forsterite.isothermal_bulk_modulus_reuss / 1.0e9 - Kt[i]) / Kt[i]
        )
        Ks_comp[i] = (
            100.0 * (forsterite.isentropic_bulk_modulus_reuss / 1.0e9 - Ks[i]) / Ks[i]
        )
        G_comp[i] = 100.0 * (forsterite.shear_modulus / 1.0e9 - G[i]) / G[i]
        VB_comp[i] = 100.0 * (forsterite.v_phi / 1000.0 - VB[i]) / VB[i]
        VS_comp[i] = 100.0 * (forsterite.v_s / 1000.0 - VS[i]) / VS[i]
        VP_comp[i] = 100.0 * (forsterite.v_p / 1000.0 - VP[i]) / VP[i]
        vol_comp[i] = 100.0 * (forsterite.molar_volume * 1.0e6 - vol[i]) / vol[i]
        alpha_comp[i] = (
            100.0 * (forsterite.thermal_expansivity / 1.0e-5 - alpha[i]) / (alpha[-1])
        )
        Cp_comp[i] = (
            100.0
            * (
                forsterite.molar_heat_capacity_p
                / forsterite.params["molar_mass"]
                / 1000.0
                - Cp[i]
            )
            / (Cp[-1])
        )
        gr_comp[i] = (forsterite.grueneisen_parameter - gr[i]) / gr[i]

    plt.plot(temperature, rho_comp, label=r"$\\rho$")
    plt.plot(temperature, Kt_comp, label=r"$K_S$")
    plt.plot(temperature, Ks_comp, label=r"$K_T$")
    plt.plot(temperature, G_comp, label=r"$G$")
    plt.plot(temperature, VS_comp, label=r"$V_S$")
    plt.plot(temperature, VP_comp, label=r"$V_P$")
    plt.plot(temperature, VB_comp, label=r"$V_\\phi$")
    plt.plot(temperature, vol_comp, label=r"$V$")
    plt.plot(temperature, alpha_comp, label=r"$\\alpha$")
    plt.plot(temperature, Cp_comp, label=r"$c_P$")
    plt.plot(temperature, gr_comp, label=r"$\\gamma$")

    plt.xlim([0, 2200])
    plt.ylim([-0.002, 0.002])
    plt.yticks([-0.002, -0.001, 0, 0.001, 0.002])
    plt.xticks([0, 800, 1600, 2200])
    plt.xlabel("Temperature (K)")
    plt.ylabel("Difference (\\%)")
    plt.legend(loc="lower center", prop=prop, ncol=4)
    if "RUNNING_TESTS" not in globals():
        plt.savefig("benchmark1.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    check_slb_fig7_txt()
