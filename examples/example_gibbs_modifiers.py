# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_gibbs_modifiers
-----------------------

This example script demonstrates the modifications to
the gibbs free energy (and derivatives) that can be applied
as masks over the results from the equations of state.

These modifications currently take the forms:

* Landau corrections (implementations of Putnis (1992)
  and Holland and Powell (2011)
* Bragg-Williams corrections
  (implementation of Holland and Powell (1996))
* Linear (a simple delta_E + delta_V*P - delta_S*T
* Magnetic (Chin, Hertzman and Sundman (1987))

*Uses:*

* :doc:`mineral_database`


*Demonstrates:*

* creating a mineral with excess contributions
* calculating thermodynamic properties

"""
from __future__ import absolute_import

# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import numpy as np
import matplotlib.pyplot as plt


# Here we import the relevant modules from BurnMan.  The burnman
# module imports several of the most important functionalities of
# the library, including the ability to make composites, and compute
# thermoelastic properties of them.  The minerals module includes
# the mineral physical parameters for the predefined minerals in
# BurnMan
import burnman
from burnman import minerals


if __name__ == "__main__":
    # Here we show the interesting features of Landau transitions
    # First, lets look at the P wave velocity in stishovite as it transforms
    # to the CaCl2 structure at high pressure
    stv = minerals.SLB_2011.stishovite()

    T = 1500.0
    pressures = np.linspace(60.0e9, 80.0e9, 101)
    temperatures = pressures * 0.0 + 1500.0
    v_ps = stv.evaluate(["p_wave_velocity"], pressures, temperatures)[0]

    plt.plot(pressures / 1.0e9, v_ps / 1.0e3, label="stishovite")
    plt.xlabel("P (GPa)")
    plt.ylabel("$V_p$ (km/s)")
    plt.legend(loc="lower right")
    plt.show()

    # Landau transitions also cause spikes in heat capacity
    # Here we show an example of troilite, as implemented by
    # Evans et al. (2010) and incorporated into the dataset
    # of Holland and Powell (2011)

    # Here we show you how to create a mineral with a
    # Landau transition.
    # A special feature of burnman is that you can have
    # more than one Landau (or any other type of)
    # contribution.
    # Here's a copy of lot (low-temperature troilite) from
    # Holland and Powell (2011), with the Landau transition
    # of tro also included.
    from burnman.utils.chemistry import dictionarize_formula, formula_mass

    class lot(burnman.Mineral):
        def __init__(self):
            formula = "Fe1.0S1.0"
            formula = dictionarize_formula(formula)
            self.params = {
                "name": "lot",
                "formula": formula,
                "equation_of_state": "hp_tmt",
                "H_0": -102160.0,
                "S_0": 60.0,
                "V_0": 1.818e-05,
                "Cp": [50.2, 0.011052, -940000.0, 0.0],
                "a_0": 4.93e-05,
                "K_0": 65800000000.0,
                "Kprime_0": 4.17,
                "Kdprime_0": -6.3e-11,
                "n": sum(formula.values()),
                "molar_mass": formula_mass(formula),
            }
            self.property_modifiers = [
                [
                    "landau_hp",
                    {
                        "P_0": 100000.0,
                        "T_0": 298.15,
                        "Tc_0": 420.0,
                        "S_D": 10.0,
                        "V_D": 0.0,
                    },
                ],
                [
                    "landau_hp",
                    {
                        "P_0": 100000.0,
                        "T_0": 298.15,
                        "Tc_0": 598.0,
                        "S_D": 12.0,
                        "V_D": 4.1e-7,
                    },
                ],
            ]
            burnman.Mineral.__init__(self)

    troilite = lot()
    lot = minerals.HP_2011_ds62.lot()
    tro = minerals.HP_2011_ds62.tro()

    temperatures = np.linspace(300.0, 1300.0, 101)
    pressures = temperatures * 0.0 + 1.0e5

    C_ps_troilite = troilite.evaluate(["C_p"], pressures, temperatures)[0]
    C_ps_lot = lot.evaluate(["C_p"], pressures, temperatures)[0]
    C_ps_tro = tro.evaluate(["C_p"], pressures, temperatures)[0]

    plt.plot(temperatures, C_ps_lot, "r--", label="low temperature (HP2011)")
    plt.plot(temperatures, C_ps_tro, "g--", label="high temperature (HP2011)")
    plt.plot(temperatures, C_ps_troilite, "b-", label="troilite")
    plt.xlabel("T (K)")
    plt.ylabel("$C_p$ (J/K/mol)")
    plt.legend(loc="lower right")
    plt.show()

    # Spinel is a mineral with a Bragg-Williams type model
    sp = minerals.HP_2011_ds62.sp()

    sp2 = minerals.HP_2011_ds62.sp()
    sp2.property_modifiers = []

    C_ps_sp = sp.evaluate(["C_p"], pressures, temperatures)[0]
    C_ps_sp2 = sp2.evaluate(["C_p"], pressures, temperatures)[0]

    plt.plot(
        temperatures, C_ps_sp2, linestyle="--", label="spinel without B-W transition"
    )
    plt.plot(temperatures, C_ps_sp, label="spinel with B-W transition")
    plt.xlabel("T (K)")
    plt.ylabel("$C_p$ (J/K/mol)")
    plt.legend(loc="lower right")
    plt.show()

    # The Holland and Powell nepheline model
    # has a Landau-type transition at low temperature,
    # but it is also possible to simplify the model by just having an
    # excess energy, entropy and volume to estimate the
    # thermal properties at high temperature.

    ne_HP = burnman.minerals.HP_2011_ds62.ne()
    ne_HP2 = burnman.minerals.HP_2011_ds62.ne()
    ne_HP2.property_modifiers = []

    ne_HP.set_state(1.0e5, 1000.0)
    ne_HP2.set_state(1.0e5, 1000.0)

    ne_HP2.property_modifiers = [
        [
            "linear",
            {
                "delta_E": (ne_HP.molar_internal_energy - ne_HP2.molar_internal_energy),
                "delta_S": ne_HP.S - ne_HP2.S,
                "delta_V": ne_HP.V - ne_HP2.V,
            },
        ]
    ]

    temperatures = np.linspace(300.0, 800.0, 101)
    pressures = temperatures * 0.0 + 1.0e5

    Ss_HP = ne_HP.evaluate(["molar_entropy"], pressures, temperatures)[0]
    Ss_HP2 = ne_HP2.evaluate(["molar_entropy"], pressures, temperatures)[0]

    plt.plot(temperatures, Ss_HP, label="nepheline (HP2011)")
    plt.plot(temperatures, Ss_HP2, label="nepheline (disordered)")
    plt.xlabel("T (K)")
    plt.ylabel("S (J/K/mol)")
    plt.legend(loc="lower right")
    plt.show()
