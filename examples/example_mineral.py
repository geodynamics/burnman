# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_mineral
---------------

This example shows how to create mineral objects in BurnMan,
and how to output their thermodynamic and thermoelastic quantities.

Mineral objects are the building blocks for more complex objects in BurnMan.
These objects are intended to represent minerals (or melts, or fluids)
of fixed composition, with a well defined equation of state that defines
the relationship between the current state (pressure and temperature) of
the mineral and its thermodynamic potentials and derivatives
(such as volume and entropy).

Mineral objects are initialized with a dictionary containing all of the
parameters required by the desired equation of state. BurnMan contains
implementations of may equations of state (:doc:`eos`).


*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Mineral`


*Demonstrates:*

* Different ways to define an endmember
* How to set state
* How to output thermodynamic and thermoelastic properties

"""
from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman.minerals import HGP_2018_ds633, SLB_2011
from burnman import CombinedMineral

from burnman.utils.chemistry import dictionarize_formula, formula_mass
from burnman.utils.chemistry import formula_to_string


if __name__ == "__main__":
    """
    Below, we demonstrate the creation of a forsterite object for the
    Stixrude and Lithgow-Bertelloni (2011) equation of state which uses
    a 3rd order expansion for the shear modulus (equation_of_state = 'slb3').
    """

    forsterite_formula = dictionarize_formula("Mg2SiO4")
    forsterite_params = {
        "name": "Forsterite",
        "formula": forsterite_formula,
        "equation_of_state": "slb3",
        "F_0": -2055403.0,
        "V_0": 4.3603e-05,
        "K_0": 1.279555e11,
        "Kprime_0": 4.21796,
        "Debye_0": 809.1703,
        "grueneisen_0": 0.99282,
        "q_0": 2.10672,
        "G_0": 81599990000.0,
        "Gprime_0": 1.46257,
        "eta_s_0": 2.29972,
        "n": sum(forsterite_formula.values()),
        "molar_mass": formula_mass(forsterite_formula),
    }

    forsterite = burnman.Mineral(params=forsterite_params)

    """
    BurnMan also contains the Stixrude and Lithgow-Bertelloni (2011)
    dataset, so we can also create an object corresponding to this mineral
    directly.
    """

    forsterite = SLB_2011.forsterite()

    """
    In petrology, we are often interested in phases for which we have little
    experimental or theoretical information.
    One common example is when we want to approximate the properties
    of an ordered phase relative to its disordered counterparts.
    In many cases, a reasonable procedure is to make a mechanical mixture
    of the known phases, such that they are the correct composition of
    the unknown phase, and then apply a linear correction to the Gibbs
    energy of the phase (i.e. Delta G = Delta E - T Delta S + P Delta V).
    In BurnMan, we do this using the CombinedMineral class.
    In the lines below, we create an ordered orthopyroxene with composition
    MgFeSi2O6 from a 50:50 mixture of enstatite and ferrosilite.
    We make this compound 6 kJ/mol more stable than the mechanical mixture.
    """

    fe_mg_orthopyroxene = CombinedMineral(
        name="ordered ferroenstatite",
        mineral_list=[HGP_2018_ds633.en(), HGP_2018_ds633.fs()],
        molar_amounts=[0.5, 0.5],
        free_energy_adjustment=[-6.0e3, 0.0, 0.0],
    )
    print(
        f"Formula of CombinedMineral {fe_mg_orthopyroxene.name}: "
        f"{formula_to_string(fe_mg_orthopyroxene.formula)}\n"
    )

    """
    In the following lines, we will compare the properties of two
    implementations of the forsterite mineral.
    """

    fo_SLB = SLB_2011.forsterite()
    fo_HP = HGP_2018_ds633.fo()

    # Overwrite names
    fo_SLB.name = "Forsterite (SLB2011)"
    fo_HP.name = "Forsterite (HGP2018)"

    # Create list of minerals
    minerals = [fo_SLB, fo_HP]

    """
    Once an object has been initialised, the user can set its state
    (pressure in Pa and temperature in K):
    """

    P = 1.0e5
    T = 1000.0

    for m in minerals:
        m.set_state(1.0e9, 1000)
        print(
            f"{m.name}\n"
            f"  Gibbs: {m.gibbs:.2f} J/mol\n"
            f"  Entropy: {m.molar_entropy:.2f} J/K/mol\n"
            f"  Volume: {m.molar_volume*1.e6:.2f} cm3/mol"
        )

    """
    It is common for users to want to return properties over a range of states.
    BurnMan makes this easy through the "evaluate" method of the Material class.
    """

    temperatures = np.linspace(300.0, 1300.0, 101)
    pressures = 1.0e5 * np.ones_like(temperatures)
    Cp_HP, V_HP = fo_HP.evaluate(
        ["molar_heat_capacity_p", "molar_volume"], pressures, temperatures
    )
    Cp_SLB, V_SLB = fo_SLB.evaluate(
        ["molar_heat_capacity_p", "molar_volume"], pressures, temperatures
    )

    # The following lines do the plotting
    fig = plt.figure(figsize=(8, 4))
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    ax[0].plot(temperatures, Cp_HP, label="HGP2018, 1 bar")
    ax[0].plot(temperatures, Cp_SLB, label="SLB2011, 1 bar")

    ax[1].plot(temperatures, V_HP * 1.0e6, label="HGP2018, 1 bar")
    ax[1].plot(temperatures, V_SLB * 1.0e6, label="SLB2011, 1 bar")

    for i in range(2):
        ax[i].set_xlabel("Temperature (K)")
        ax[i].legend()

    ax[0].set_ylabel("Molar isobaric heat capacity (J/K/mol)")
    ax[1].set_ylabel("Molar volume (cm$^3$/mol)")
    fig.set_layout_engine("tight")
    plt.show()
