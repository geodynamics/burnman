# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_spintransition_thermal
------------------------------

This example illustrates how to create a non-ideal solution model
for (Mg,Fe\\ :sup:`HS`\\ ,Fe\\ :sup:`LS`\\ )O ferropericlase that has a gradual
spin transition at finite temperature.

In the first part of the example, we create a Solution class from scratch
that incorporates three endmembers: MgO, FeO in the high spin state,
and FeO in the low spin state. We then embed this solution in a
RelaxedSolution class that allows us to define the isochemical high spin
to low spin transition as a free relaxation vector that is allowed to vary freely
to minimize the Gibbs free energy.

In the second part of the example, we illustrate how to use the
ferropericlase solution model directly from the Stixrude and Lithgow-Bertelloni (2024)
database, again embedding it in a RelaxedSolution to incorporate the spin transition.
We then show how to calculate thermodynamic and elastic properties
that incorporate the effects of the spin transition over a range of
pressures and temperatures along an oceanic geotherm.

In this example, we implicitly apply the Bragg-Williams approximation
(i.e., we assume that there is no short-range order by only incorporating
interactions that are a function of the average occupancy of species
on each distinct site). Furthermore, the implemented mixing model
explicitly precludes long range order.

*Specifically uses:*

* :class:`burnman.Mineral`
* :class:`burnman.Solution`
* :class:`burnman.RelaxedSolution`

*Demonstrates:*

* implementation of gradual spin transition in (Mg,Fe)O at a user-defined
  pressure and temperature
"""

import numpy as np
import matplotlib.pyplot as plt

from burnman import Solution, RelaxedSolution, seismic, geotherm
from burnman.classes.solutionmodel import SymmetricRegularSolution
from burnman.minerals import SLB_2024

if __name__ == "__main__":
    """
    First, we create three Mineral objects: one for periclase,
    one for high spin wuestite and another for low spin wuestite.

    These endmembers are taken directly from
    Stixrude and Lithgow-Bertelloni (2024).
    """

    periclase = SLB_2024.periclase()
    high_spin_wuestite = SLB_2024.wustite()
    low_spin_wuestite = SLB_2024.wustite_low_spin()

    """
    Now, we create a standard solution class for ferropericlase
    based on the full model from Stixrude and Lithgow-Bertelloni (2024).
    """

    class ferropericlase(Solution):
        def __init__(self, molar_fractions=None):
            self.name = "ferropericlase"
            self.solution_model = SymmetricRegularSolution(
                endmembers=[
                    [periclase, "[Mg]2[Mg]2"],
                    [high_spin_wuestite, "[Fe]2[Fe]2"],
                    [low_spin_wuestite, "[Fels]2[Fels]2"],
                ],
                energy_interaction=[
                    [44000.0, -87120.47],
                    [-60219.09],
                ],
                volume_interaction=[
                    [4.4e-07, 0.0],
                    [0.0],
                ],
            )

            Solution.__init__(self, molar_fractions=molar_fractions)

    """
    Now, by itself this solution class allows the user to set the proportions
    of MgO, FeO(HS) and FeO(LS) arbitrarily. However, we want to
    define a method that allows us to set the equilibrium proportions of
    FeO(HS) and FeO(LS) at a given pressure, temperature and bulk composition.
    To do this, we embed this solution in another class called a
    RelaxedSolution, which allows us to define relaxation vectors
    corresponding to the spin transition reaction:
    """
    fper = ferropericlase()
    fper_relaxed = RelaxedSolution(
        ferropericlase(),
        relaxation_vectors=np.array([[0.0, -1.0, 1.0]]),
        unrelaxed_vectors=np.array([[1.0, 0.0, 0.0], [0.0, 0.5, 0.5]]),
    )

    # Now we loop over a series of pressures at three different temperatures,
    # calculating the equilibrium composition of the solution at each.
    # We fix the bulk composition of the solution to be (Mg0.8Fe0.2)O.
    X_Fe = 0.2

    pressures = np.linspace(10.0e9, 150.0e9, 101)
    volumes = np.empty_like(pressures)
    volumes_HS = np.empty_like(pressures)
    volumes_LS = np.empty_like(pressures)
    p_LS = np.empty_like(pressures)

    fig = plt.figure(figsize=(8, 4))
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    for T, color in [[300.0, "blue"], [1800.0, "purple"], [3300.0, "red"]]:
        for i, P in enumerate(pressures):
            # Calculate and store the equilibrium volume and proportion of
            # iron in the low spin state
            fper_relaxed.set_state(P, T)
            fper_relaxed.set_composition([1.0 - X_Fe, X_Fe])

            volumes[i] = fper_relaxed.V / 4.0

            # Within the relaxed solution, the unrelaxed solution is accessible
            # via the .unrelaxed attribute. This allows us to extract the
            # molar fractions of the unrelaxed endmembers, and thus
            # the high spin and low spin proportions of iron.
            p_LS[i] = fper_relaxed.unrelaxed.molar_fractions[2] / (
                fper_relaxed.unrelaxed.molar_fractions[1]
                + fper_relaxed.unrelaxed.molar_fractions[2]
            )

            # Also calculate and store the volumes if all iron were in the
            # high or low spin state
            fper.set_state(P, T)
            fper.set_composition([1.0 - X_Fe, X_Fe, 0.0])
            volumes_HS[i] = fper.V / 4.0
            fper.set_composition([1.0 - X_Fe, 0.0, X_Fe])
            volumes_LS[i] = fper.V / 4.0

        # Do some plotting
        ax[0].fill_between(
            pressures / 1.0e9,
            volumes_HS * 1.0e6,
            volumes_LS * 1.0e6,
            alpha=0.15,
            color=color,
            label=f"{T} K, volume range",
        )
        ax[0].plot(
            pressures / 1.0e9,
            volumes * 1.0e6,
            c=color,
            linewidth=2,
            label=f"{T} K, equilibrium volume",
        )

        ax[1].plot(pressures / 1.0e9, 1.0 - p_LS, c=color, label=f"{T} K")

    # Add legends and axis titles to the plot
    for i in range(2):
        ax[i].legend()
        ax[i].set_xlabel("Pressure (GPa)")

    ax[0].set_ylabel("Volume (cm$^3$/mol)")
    ax[1].set_ylabel("High spin fraction")

    ax[0].set_ylim(7, 13)

    # Tidy the plot and show it
    fig.set_layout_engine("tight")
    plt.show()

    # The RelaxedSolution class is particularly powerful because
    # it allows us to calculate not only the equilibrium composition
    # at given P and T, but also the relaxed thermodynamic and elastic
    # properties that incorporate the effects of the spin transition.

    # This includes the direct effects on properties such as volume,
    # and entropy, as well as the effects on properties related to
    # the second derivatives of the Gibbs free energy, such as
    # seismic velocities and heat capacity. These second derivatives
    # are affected not only by the equilibrium proportions of high spin
    # and low spin iron, but also by the rate of change of these proportions
    # with respect to pressure and temperature.

    # For this part, we could use the same ferropericlase solution model
    # defined above, but here we illustrate the same approach
    # using the ferropericlase model taken directly from the
    # Stixrude and Lithgow-Bertelloni (2024) database.

    # Create a high spin ferropericlase solution
    fper_HS = SLB_2024.ferropericlase()
    fper_HS.set_composition([0.6, 0.4, 0.0, 0.0, 0.0])

    # Create a low spin ferropericlase solution
    fper_LS = SLB_2024.ferropericlase()
    fper_LS.set_composition([0.6, 0.0, 0.4, 0.0, 0.0])

    # Create the relaxed solution
    fper_relaxed = RelaxedSolution(
        SLB_2024.ferropericlase(),
        relaxation_vectors=np.array([[0.0, -1.0, 1.0, 0.0, 0.0]]),
        unrelaxed_vectors=np.array(
            [[1.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.5, 0.5, 0.0, 0.0]]
        ),
    )
    fper_relaxed.set_composition([0.6, 0.4])

    # Now, we evaluate the seismic velocities and density of
    # high spin, low spin, and relaxed ferropericlase over a range
    # of pressures and temperatures corresponding to Stacey's oceanic
    # geotherm.
    pressures = np.linspace(0.0e9, 135e9, 136)
    depths = seismic.PREM().depth(pressures)
    temperatures = geotherm.stacey_oceanic(depths)
    Vp_h, Vs_h, rho_h = fper_HS.evaluate(["v_p", "v_s", "rho"], pressures, temperatures)
    Vp_l, Vs_l, rho_l = fper_LS.evaluate(["v_p", "v_s", "rho"], pressures, temperatures)
    Vp_r, Vs_r, rho_r, mf_r = fper_relaxed.evaluate(
        ["v_p", "v_s", "rho", "molar_fractions"], pressures, temperatures
    )

    # Also calculate unrelaxed properties at the relaxed molar fractions for comparison
    fper = SLB_2024.ferropericlase()
    Vp_u, Vs_u, rho_u = fper.evaluate(
        ["v_p", "v_s", "rho"], pressures, temperatures, mf_r
    )

    # Finally, we plot the results
    fig = plt.figure(figsize=(10, 6))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]
    # Plot temperature profile
    ax[0].plot(depths / 1.0e3, temperatures, c="black", label="PREM geotherm")
    ax[0].set_ylabel("Temperature (K)")
    ax[0].set_xlabel("Depth (km)")
    ax[0].legend()
    # Plot pressure profile
    ax[1].plot(depths / 1.0e3, pressures / 1.0e9, c="black", label="PREM geotherm")
    ax[1].set_ylabel("Pressure (GPa)")
    ax[1].set_xlabel("Depth (km)")
    ax[1].legend()
    # Plot seismic velocities and density
    HS_colour = "blue"
    LS_colour = "red"
    MS_colour = "purple"
    relaxed_style = "-"
    unrelaxed_style = ":"
    HS_fractions = mf_r[:, 1] / (mf_r[:, 1] + mf_r[:, 2])
    ax[2].plot(depths / 1.0e3, 1.0 * np.ones_like(depths), c=HS_colour)
    ax[2].plot(depths / 1.0e3, 0.0 * np.ones_like(depths), c=LS_colour)
    ax[2].plot(depths / 1.0e3, HS_fractions, linewidth=2, c=MS_colour)
    ax[2].set_ylabel("High spin fraction")
    ax[2].set_xlabel("Depth (km)")

    ax[3].plot(depths / 1.0e3, Vp_h / 1.0e3, c=HS_colour)
    ax[3].plot(depths / 1.0e3, Vp_l / 1.0e3, c=LS_colour)
    ax[3].plot(
        depths / 1.0e3, Vp_r / 1.0e3, c=MS_colour, linewidth=2, linestyle=relaxed_style
    )
    ax[3].plot(
        depths / 1.0e3,
        Vp_u / 1.0e3,
        c=MS_colour,
        linewidth=2,
        linestyle=unrelaxed_style,
    )
    ax[3].set_ylabel("$V_P$ (km/s)")

    ax[4].plot(depths / 1.0e3, Vs_h / 1.0e3, c=HS_colour)
    ax[4].plot(depths / 1.0e3, Vs_l / 1.0e3, c=LS_colour)
    ax[4].plot(
        depths / 1.0e3, Vs_r / 1.0e3, c=MS_colour, linewidth=2, linestyle=relaxed_style
    )
    ax[4].plot(
        depths / 1.0e3,
        Vs_u / 1.0e3,
        c=MS_colour,
        linewidth=2,
        linestyle=unrelaxed_style,
    )
    ax[4].set_ylabel("$V_S$ (km/s)")

    ax[5].plot(depths / 1.0e3, rho_h, c=HS_colour, label="HS")
    ax[5].plot(depths / 1.0e3, rho_l, c=LS_colour, label="LS")
    ax[5].plot(
        depths / 1.0e3,
        rho_r,
        c=MS_colour,
        linewidth=2,
        linestyle=relaxed_style,
        label="relaxed",
    )
    ax[5].plot(
        depths / 1.0e3,
        rho_u,
        c=MS_colour,
        linewidth=2,
        linestyle=unrelaxed_style,
        label="unrelaxed",
    )
    ax[5].set_ylabel("Density (kg/m$^3$)")
    ax[5].set_xlabel("Depth (km)")
    ax[5].legend()

    fig.set_layout_engine("tight")
    plt.show()
