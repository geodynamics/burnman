# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_assemblage_relaxation
-----------------------------

In this example, we demonstrate how to use the `RelaxedComposite` class
to compute the properties of a mineral assemblage allowing for
different levels of compositional and phase relaxation.

The specific example used here is a binary assemblage of olivine and
wadsleyite at 1600 K, over the pressure range where both phases are
stable. We compute and compare the isothermal bulk modulus, thermal
expansivity and specific heat capacity of the assemblage in three cases:
1) No relaxation (standard Composite)
2) Partial relaxation (Fe-Mg exchange is allowed between olivine and wadsleyite)
3) Full relaxation (Fe-Mg exchange + phase proportions adjust to equilibrium)

One scientific motivation for this example is to illustrate how
different relaxation mechanisms can impact the thermodynamic
properties of mineral assemblages, which is important for:
1) Interpreting geophysical observations of Earth's interior.
   For example, should we expect a significant decrease in seismic
   velocities across the olivine-wadsleyite transition zone if
   reactions can take place over seismic timescales?
2) Understanding mantle convection and dynamics.
   For example, might the rate of growth of wadsleyite from olivine
   affect patterns of mantle convection?

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Composite`
* :func:`burnman.equilibrate`
"""
import numpy as np
import matplotlib.pyplot as plt

from burnman import equilibrate
from burnman.minerals import SLB_2011
from burnman import Composite, RelaxedComposite

if __name__ == "__main__":
    # Initialize the minerals we will use in this example.
    ol = SLB_2011.mg_fe_olivine()
    wad = SLB_2011.mg_fe_wadsleyite()

    # Create the assemblage with arbitrary proportions
    assemblage = Composite([ol, wad], [1.0, 0.0])
    ol.set_composition([0.90, 0.10])
    wad.set_composition([0.90, 0.10])

    # Find the wad-in and ol-in pressures at 1600 K
    T = 1600.0
    composition = assemblage.formula.copy()
    equality_constraints = [("T", T), ("phase_fraction", (wad, 0.0))]
    equilibrate(composition, assemblage, equality_constraints)
    P_wad_in = assemblage.pressure

    equality_constraints = [("T", T), ("phase_fraction", (ol, 0.0))]
    equilibrate(composition, assemblage, equality_constraints)
    P_ol_in = assemblage.pressure

    # Now get the unrelaxed properties over the two phase region
    # at their equilibrated states
    pressures = np.linspace(P_wad_in + 100.0, P_ol_in - 100.0, 101)
    equality_constraints = [("P", pressures), ("T", T)]

    sols, prm = equilibrate(composition, assemblage, equality_constraints)
    KS_u = np.array([sol.assemblage.K_S for sol in sols])
    alpha_u = np.array([sol.assemblage.alpha for sol in sols])
    Cp_u = np.array([sol.assemblage.C_p for sol in sols])
    M_u = np.array([sol.assemblage.molar_mass for sol in sols])

    # Create relaxed composites allowing for different levels of relaxation
    assemblage_relaxed = RelaxedComposite(assemblage, assemblage.reaction_basis)
    assemblage_partially_relaxed = RelaxedComposite(
        assemblage, [[1.0, -1.0, -1.0, 1.0]]
    )

    # Get the relaxed properties over the two phase region, again at their
    # equilibrated states
    temperatures = T * np.ones_like(pressures)
    KS_r, alpha_r, Cp_r, M_r = assemblage_relaxed.evaluate(
        ["K_S", "alpha", "C_p", "molar_mass"], pressures, temperatures
    )
    KS_pr, alpha_pr, Cp_pr, M_pr = assemblage_partially_relaxed.evaluate(
        ["K_S", "alpha", "C_p", "molar_mass"], pressures, temperatures
    )

    # Print some results to the screen
    print("Following results are for an olivine-wadsleyite assemblage at 1600 K.")
    print(
        "Order of results: Unrelaxed, Partial relaxation (Fe-Mg exchange), Full relaxation\n"
    )

    print(f"Start of ol-wad transition zone at {P_ol_in/1.0e9:.2f} GPa")
    print(
        f"Isentropic bulk modulus: {KS_u[0]/1.0e9:.2f} GPa, {KS_pr[0]/1.0e9:.2f} GPa, {KS_r[0]/1.0e9:.2f} GPa"
    )
    print(
        f"Thermal expansivity: {alpha_u[0]*1.0e6:.2f}e-6 /K, {alpha_pr[0]*1.0e6:.2f}e-6 /K, {alpha_r[0]*1.0e6:.2f}e-6 /K"
    )
    print(
        f"Specific heat capacity: {Cp_u[0]/M_u[0]:.1f} J/K/kg, {Cp_pr[0]/M_pr[0]:.1f} J/K/kg, {Cp_r[0]/M_r[0]:.1f} J/K/kg"
    )
    print()

    print(
        f"Middle of ol-wad transition zone at {(P_ol_in + P_wad_in)/2.0/1.0e9:.2f} GPa"
    )
    mid_index = len(pressures) // 2
    print(
        f"Isentropic bulk modulus: {KS_u[mid_index]/1.0e9:.2f} GPa, {KS_pr[mid_index]/1.0e9:.2f} GPa, {KS_r[mid_index]/1.0e9:.2f} GPa"
    )
    print(
        f"Thermal expansivity: {alpha_u[mid_index]*1.0e6:.2f}e-6 /K, {alpha_pr[mid_index]*1.0e6:.2f}e-6 /K, {alpha_r[mid_index]*1.0e6:.2f}e-6 /K"
    )
    print(
        f"Specific heat capacity: {Cp_u[mid_index]/M_u[mid_index]:.1f} J/K/kg, {Cp_pr[mid_index]/M_pr[mid_index]:.1f} J/K/kg, {Cp_r[mid_index]/M_r[mid_index]:.1f} J/K/kg"
    )
    print()

    print(f"End of ol-wad transition zone at {P_wad_in/1.0e9:.2f} GPa")
    print(
        f"Isentropic bulk modulus: {KS_u[-1]/1.0e9:.2f} GPa, {KS_pr[-1]/1.0e9:.2f} GPa, {KS_r[-1]/1.0e9:.2f} GPa"
    )
    print(
        f"Thermal expansivity: {alpha_u[-1]*1.0e6:.2f}e-6 /K, {alpha_pr[-1]*1.0e6:.2f}e-6 /K, {alpha_r[-1]*1.0e6:.2f}e-6 /K"
    )
    print(
        f"Specific heat capacity: {Cp_u[-1]/M_u[-1]:.1f} J/K/kg, {Cp_pr[-1]/M_pr[-1]:.1f} J/K/kg, {Cp_r[-1]/M_r[-1]:.1f} J/K/kg"
    )

    fig = plt.figure(figsize=(15, 5))
    ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]
    ax[0].plot(pressures / 1.0e9, KS_u / 1.0e9, label="Unrelaxed")
    ax[0].plot(pressures / 1.0e9, KS_pr / 1.0e9, linestyle="--")
    ax[0].plot(pressures / 1.0e9, KS_r / 1.0e9)
    ax[0].set_xlabel("Pressure (GPa)")
    ax[0].set_ylabel("K$_S$ (GPa)")

    ax[1].plot(pressures / 1.0e9, alpha_u * 1.0e6, label="Unrelaxed")
    ax[1].plot(pressures / 1.0e9, alpha_pr * 1.0e6, linestyle="--")
    ax[1].plot(pressures / 1.0e9, alpha_r * 1.0e6)
    ax[1].set_xlabel("Pressure (GPa)")
    ax[1].set_ylabel("alpha (1e-6/K)")

    ax[2].plot(pressures / 1.0e9, Cp_u / M_u, label="Unrelaxed")
    ax[2].plot(
        pressures / 1.0e9, Cp_pr / M_pr, linestyle="--", label="Fe-Mg exchange only"
    )
    ax[2].plot(pressures / 1.0e9, Cp_r / M_r, label="Fe-Mg exchange + phase change")
    ax[2].set_xlabel("Pressure (GPa)")
    ax[2].set_ylabel("C$_p$ (J/K/kg)")
    ax[2].legend()
    plt.show()

    # The results show that allowing for full relaxation
    # of the assemblage leads to significant reductions in the
    # isentropic bulk modulus and increases in thermal expansivity
    # and specific heat capacity across the olivine-wadsleyite
    # transition. In contrast, allowing only for Fe-Mg exchange
    # between the two phases has a much smaller effect on these properties.
