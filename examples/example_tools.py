# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
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
import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman.tools.chemistry import equilibrium_pressure
from burnman.tools.chemistry import equilibrium_temperature
from burnman.tools.chemistry import invariant_point
from burnman.tools.chemistry import hugoniot
from burnman.tools.eos import check_eos_consistency


def round_to_n(x, xerr, n):
    return round(x, -int(np.floor(np.log10(np.abs(xerr)))) + (n - 1))


if __name__ == "__main__":
    # First, let's check the EoS consistency of SLB_2011 periclase
    check_eos_consistency(
        burnman.minerals.SLB_2011.periclase(), P=10.0e9, T=3000.0, verbose=True
    )
    print("")

    # Next, let's create the Mg2SiO4 phase diagram
    forsterite = burnman.minerals.HP_2011_ds62.fo()
    liquid = burnman.minerals.HP_2011_ds62.foL()
    mg_wadsleyite = burnman.minerals.HP_2011_ds62.mwd()
    mg_ringwoodite = burnman.minerals.HP_2011_ds62.mrw()
    periclase = burnman.minerals.HP_2011_ds62.per()
    mg_perovskite = burnman.minerals.HP_2011_ds62.mpv()

    fo_wad_liq = invariant_point(
        [forsterite, mg_wadsleyite],
        [1.0, -1.0],
        [forsterite, liquid],
        [1.0, -1],
        pressure_temperature_initial_guess=[14.0e9, 2000.0],
    )

    # Here's one example where we find the equilibrium temperature:
    P = 14.0e9
    T = equilibrium_temperature([forsterite, mg_wadsleyite], [1.0, -1.0], P)
    print("Endmember equilibrium calculations")
    print(
        "fo -> wad equilibrium at",
        P / 1.0e9,
        "GPa is reached at",
        round_to_n(T, T, 4),
        "K",
    )
    print("")

    # Here we find the equilibrium temperature for forsterite melting
    pressures = np.linspace(1.0e5, fo_wad_liq[0], 11)

    T_fo_melt = np.empty_like(pressures)
    for i, P in enumerate(pressures):
        T_fo_melt[i] = equilibrium_temperature(
            [forsterite, liquid], [1.0, -1.0], P, 2000.0
        )

    # Now let's make the whole diagram using equilibrium_pressure
    temperatures = np.linspace(1000.0, fo_wad_liq[1], 21)
    pressures_fo_wd = np.empty_like(temperatures)
    pressures_wd_rw = np.empty_like(temperatures)
    pressures_rw_perpv = np.empty_like(temperatures)

    for i, T in enumerate(temperatures):
        P = equilibrium_pressure([forsterite, mg_wadsleyite], [1.0, -1.0], T)
        pressures_fo_wd[i] = P

        P = equilibrium_pressure([mg_wadsleyite, mg_ringwoodite], [1.0, -1.0], T)
        pressures_wd_rw[i] = P

        P = equilibrium_pressure(
            [mg_ringwoodite, periclase, mg_perovskite], [1.0, -1.0, -1.0], T
        )
        pressures_rw_perpv[i] = P

    plt.plot(T_fo_melt, pressures / 1.0e9, label="fo -> melt")
    plt.plot(temperatures, pressures_fo_wd / 1.0e9, label="fo -> wd")
    plt.plot(temperatures, pressures_wd_rw / 1.0e9, label="wd -> rw")
    plt.plot(temperatures, pressures_rw_perpv / 1.0e9, label="rw -> per + pv")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pressure (GPa)")
    plt.legend()
    plt.title("Mg2SiO4 phase diagram")
    plt.ylim(0.0, 30.0)
    plt.show()

    # Here's a calculation of the Hugoniot of periclase up to 120 GPa
    print("Hugoniot calculations")
    pressures = np.linspace(1.0e5, 120.0e9, 101)
    temperatures, volumes = hugoniot(periclase, 1.0e5, 298.15, pressures)
    plt.plot(pressures / 1.0e9, temperatures, label="298.15 K")
    print(
        "Room temperature Hugoniot temperature at",
        pressures[-1] / 1.0e9,
        "GPa:",
        int(temperatures[-1] + 0.5),
        "K",
    )

    temperatures, volumes = hugoniot(periclase, 1.0e5, 1000.0, pressures)
    plt.plot(pressures / 1.0e9, temperatures, label="1000 K")
    print(
        "1000 K Hugoniot temperature at",
        pressures[-1] / 1.0e9,
        "GPa:",
        int(temperatures[-1] + 0.5),
        "K",
    )

    plt.legend(loc="upper left")
    plt.ylabel("Temperature (K)")
    plt.xlabel("Pressure (GPa)")
    plt.title("Periclase Hugoniot")
    plt.show()
