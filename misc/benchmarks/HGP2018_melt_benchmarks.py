# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

HGP2018_melt_benchmarks
-----------------------

This script tests the Holland et al. (2018) melt model in
the CMS and MS systems.

"""
import numpy as np
import matplotlib.pyplot as plt
from burnman import Composite
from burnman.minerals import HGP_2018_ds633
from burnman import equilibrate

if __name__ == "__main__":
    di = HGP_2018_ds633.di()
    liq = HGP_2018_ds633.CMS_melt()

    liq.set_composition([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])
    composition = {"Mg": 1.0, "Ca": 1.0, "Si": 2.0, "O": 6.0}

    assemblage = Composite([di, liq])

    equality_constraints = [
        ["P", np.linspace(1.0e5, 5.0e9, 101)],
        ["phase_fraction", (di, 0.0)],
    ]
    sols, prm = equilibrate(composition, assemblage, equality_constraints)

    plt.plot(
        [sol.assemblage.pressure / 1.0e9 for sol in sols],
        [sol.assemblage.temperature for sol in sols],
        label="Solution model diopside melting",
    )

    liq = HGP_2018_ds633.diL()
    assemblage = Composite([di, liq])

    equality_constraints = [
        ["P", np.linspace(1.0e5, 5.0e9, 101)],
        ["phase_fraction", (di, 0.0)],
    ]
    sols, prm = equilibrate(composition, assemblage, equality_constraints)

    plt.plot(
        [sol.assemblage.pressure / 1.0e9 for sol in sols],
        [sol.assemblage.temperature for sol in sols],
        linestyle=":",
        label="Raw dataset diopside melting",
    )

    plt.legend()
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Temperature (K))")
    plt.show()

    liq = HGP_2018_ds633.MS_melt()
    per = HGP_2018_ds633.per()
    fo = HGP_2018_ds633.fo()
    pren = HGP_2018_ds633.pren()
    crst = HGP_2018_ds633.crst()

    # peritectics / eutectics
    liq.set_composition([0.1, 0.9])
    composition = {"Mg": 2.0, "Si": 1.5, "O": 5.0}
    assemblage = Composite([fo, pren, liq])
    equality_constraints = [["P", 1.0e5], ["phase_fraction", (liq, 0.0)]]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)

    T_fo_pren = assemblage.temperature
    xSi_fo_pren = liq.formula["Si"] / (liq.formula["Si"] + liq.formula["Mg"])

    composition = {"Mg": 1.0, "Si": 2.0, "O": 5.0}
    assemblage = Composite([pren, crst, liq])
    equality_constraints = [["P", 1.0e5], ["phase_fraction", (liq, 0.0)]]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)

    T_pren_crst = assemblage.temperature
    xSi_pren_crst = liq.formula["Si"] / (liq.formula["Si"] + liq.formula["Mg"])

    print("fo-pren peritectic")
    print(f"x(SiO2): {xSi_fo_pren:.3f}, T: {T_fo_pren:.2f} K")

    print("pren-crst eutectic")
    print(f"x(SiO2): {xSi_pren_crst:.3f}, T: {T_pren_crst:.2f} K")

    plt.plot([1.0 / 3.0, 1.0 / 3.0], [1750.0, 2250.0], color="black")
    plt.plot([1.0 / 3.0, xSi_fo_pren], [T_fo_pren, T_fo_pren], color="black")
    plt.plot([0.5, 0.5], [1750.0, T_fo_pren], color="black")
    plt.plot([0.5, 1.0], [T_pren_crst, T_pren_crst], color="black")

    for xSi0, xSi1, phase in [
        [1.0 / 3.0, xSi_fo_pren, fo],
        [xSi_fo_pren, xSi_pren_crst, pren],
        [xSi_pren_crst, 0.999, crst],
    ]:
        xSis = np.linspace(xSi0, xSi1, 31)
        Ts = np.empty_like(xSis)
        for i, xSi in enumerate(xSis):
            composition = {"Mg": 1.0 - xSi, "Si": xSi, "O": (1.0 - xSi) + xSi * 2.0}
            assemblage = Composite([phase, liq])
            assemblage.set_state(1.0e5, 2000.0)
            equality_constraints = [["P", 1.0e5], ["phase_fraction", (phase, 0.0)]]

            sol, prm = equilibrate(composition, assemblage, equality_constraints)

            Ts[i] = sol.assemblage.temperature

        plt.plot(xSis, Ts, color="black")

    plt.xlim(0.0, 1.0)
    plt.ylim(1750.0, 2250.0)
    plt.xlabel("xSiO$_2$/(xMgO + xSiO$_2$) (mole fraction)")
    plt.ylabel("Temperature (K)")
    plt.show()
