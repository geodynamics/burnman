# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_olivine_binary
----------------------

This example demonstrates how BurnMan may be used to calculate the
equilibrium binary phase diagram for the three olivine polymorphs
(olivine, wadsleyite and ringwoodite).

The calculations use the equilibrate function. Unlike the
examples in example\\_equilibrate.py, which are constrained to a
fixed bulk composition, the bulk composition is allowed to vary
along the vector [n_Mg - n_Fe].

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Composite`
* :func:`burnman.equilibrate`
"""
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman import equilibrate
from burnman.minerals import SLB_2011


if __name__ == "__main__":
    # Initialize the minerals we will use in this example.
    ol = SLB_2011.mg_fe_olivine()
    wad = SLB_2011.mg_fe_wadsleyite()
    rw = SLB_2011.mg_fe_ringwoodite()

    # Set the starting guess compositions for each of the solutions
    ol.set_composition([0.90, 0.10])
    wad.set_composition([0.90, 0.10])
    rw.set_composition([0.80, 0.20])

    # Initialize the figure that will be used to plot the binary diagram
    fig = plt.figure()
    ax = [fig.add_subplot(1, 1, 1)]

    # Loop over three temperatures
    for T, color in [(1200.0, "blue"), (1600.0, "purple"), (2000.0, "red")]:
        # First, we find the compositions of the three phases
        # at the univariant.
        composition = {"Fe": 0.2, "Mg": 1.8, "Si": 1.0, "O": 4.0}
        assemblage = burnman.Composite([ol, wad, rw], [1.0, 0.0, 0.0])
        equality_constraints = [
            ("T", T),
            ("phase_fraction", (ol, 0.0)),
            ("phase_fraction", (rw, 0.0)),
        ]
        free_compositional_vectors = [{"Mg": 1.0, "Fe": -1.0}]

        sol, prm = equilibrate(
            composition,
            assemblage,
            equality_constraints,
            free_compositional_vectors,
            verbose=False,
        )

        if not sol.success:
            raise Exception(
                "Could not find solution for the univariant using "
                "provided starting guesses."
            )

        # We interrogate the stored copy of the assemblage for the pressure and
        # the composition of all three phases
        P_univariant = sol.assemblage.pressure
        x_fa_univariant = sol.assemblage.phases[0].molar_fractions[1]
        x_fwd_univariant = sol.assemblage.phases[1].molar_fractions[1]
        x_frw_univariant = sol.assemblage.phases[2].molar_fractions[1]

        print(f"Univariant pressure at {T:.0f} K: {P_univariant/1.e9:.3f} GPa")
        # Plot the line connecting the three phases
        ax[0].plot(
            [x_fa_univariant, x_frw_univariant],
            [P_univariant / 1.0e9, P_univariant / 1.0e9],
            color=color,
        )

        # Now solve for the stable sections of the three binary loops
        i = 0
        for d in [
            [ol, wad, np.linspace(x_fa_univariant, 0.001, 20)],
            [ol, rw, np.linspace(x_fa_univariant, 0.999, 20)],
            [wad, rw, np.linspace(x_fwd_univariant, 0.001, 20)],
        ]:
            m1, m2, x_fe_m1 = d
            assemblage = burnman.Composite([m1, m2], [1.0, 0.0])

            # Reset the compositions of the two phases to have compositions
            # close to those at the univariant point
            m1.set_composition([1.0 - x_fwd_univariant, x_fwd_univariant])
            m2.set_composition([1.0 - x_fwd_univariant, x_fwd_univariant])

            # Also set the pressure and temperature
            assemblage.set_state(P_univariant, T)

            # Here our equality constraints are temperature,
            # the phase fraction of the second phase,
            # and we loop over the composition of the first phase.
            equality_constraints = [
                ("T", T),
                (
                    "phase_composition",
                    (m1, [["Mg_A", "Fe_A"], [0.0, 1.0], [1.0, 1.0], x_fe_m1]),
                ),
                ("phase_fraction", (m2, 0.0)),
            ]

            sols, prm = equilibrate(
                composition,
                assemblage,
                equality_constraints,
                free_compositional_vectors,
                verbose=False,
            )

            # Process the solutions
            out = np.array(
                [
                    [
                        sol.assemblage.pressure,
                        sol.assemblage.phases[0].molar_fractions[1],
                        sol.assemblage.phases[1].molar_fractions[1],
                    ]
                    for sol in sols
                    if sol.success
                ]
            )
            pressures, x_fe_m1s, x_fe_m2s = out.T

            if i == 0:
                ax[0].plot(x_fe_m1s, pressures / 1.0e9, color=color, label=f"{T} K")
            else:
                ax[0].plot(x_fe_m1s, pressures / 1.0e9, color=color)
            ax[0].plot(x_fe_m2s, pressures / 1.0e9, color=color)
            ax[0].fill_betweenx(
                pressures / 1.0e9, x_fe_m1s, x_fe_m2s, color=color, alpha=0.2
            )
            i += 1

    ax[0].text(0.1, 6.0, "olivine", horizontalalignment="left")
    ax[0].text(
        0.01,
        14.2,
        "wadsleyite",
        horizontalalignment="left",
        bbox=dict(facecolor="white", edgecolor="white", boxstyle="round,pad=0.2"),
    )
    ax[0].text(0.9, 15.0, "ringwoodite", horizontalalignment="right")

    ax[0].set_xlabel("p(Fe$_2$SiO$_4$)")
    ax[0].set_ylabel("Pressure (GPa)")
    ax[0].legend()
    plt.show()
