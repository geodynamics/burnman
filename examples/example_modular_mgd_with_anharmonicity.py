# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_modular_mgd_with_anharmonicity
--------------------------------------

This example shows how to create a mineral object using
the Modular Mie-Grüneisen-Debye With Anharmonicity equation of state.

This equation of state is derived from the Modular Mie-Grüneisen-Debye equation of state,
whose use was demonstrated in example_modular_mgd.py. Like that equation of state,
the modified version demonstrated here allows users to pick and choose
different options for two key aspects of the mineral's behavior:

1. The equation of state (EOS) used to describe the mineral's thermodynamic properties at the reference temperature.
2. The model used to describe the mineral's Debye temperature as a function of volume.

In the Modular Mie-Grüneisen-Debye equation of state, anharmonicity is incorporated
via extra terms in the Helmholtz energy expression; in other words, the model is
somewhat removed from the underlying quasiharmonic form
(although it shares the same volume dependence on the reference Debye temperature).
In the modified equation of state described in this example, anharmonic effects are
more tightly coupled to the basic quasiharmonic framework, and use the thermodynamic
properties of that framework to describe the mineral's anharmonic behavior. The
model used here is based on the empirical ansatz introduced by Wu and Wentzcovitch
(2009; https://doi.org/10.1103/PhysRevB.79.104304), simplified somewhat to avoid
in internal V(P(V, T), T=0) inversion to calculate the anharmonic contributions.
This has the dual benefits of reducing computational complexity and improving
the shape of the thermodynamic property functions as a function of temperature,
especially at high volumes. The simplification still requires only one parameter,
`c_anh`, which controls the strength of the anharmonic effects, but the value of that
parameter may need to be adjusted relative to the model proposed by Wu and Wentzcovitch.

*Uses:*

* :class:`burnman.Mineral`


*Demonstrates:*

* How to create a mineral object using the Modular Mie-Grüneisen-Debye equation of state.

"""
# First, a couple of standard imports
import numpy as np
import matplotlib.pyplot as plt

# Now some specific BurnMan imports
# The Mineral class is the standard class for creating mineral objects
from burnman import Mineral

# For this example, we will use the SLB_2011 periclase mineral
# as a starting point
from burnman.minerals.SLB_2011 import periclase as mineral

# This helper function creates an Equation of State object
from burnman.eos.helper import create

# These imports contain the different models available for the EOS
from burnman.eos import debye_temperature_models

# Create the SLB mineral object
m_SLB = mineral()

# Now create a new mineral object using the modular MGD with
# anharmonicity equation of state.
params = m_SLB.params.copy()
params["equation_of_state"] = "modular_mgd_with_anharmonicity"

# In these two lines, we set the reference EOS and Debye temperature model
# to be an exact clone of the original SLB object
params["reference_eos"] = create("bm3")
params["debye_temperature_model"] = debye_temperature_models.SLB()

# Now we just need to set the anharmonicity parameter
# Positive values reduce the isochoric heat capacity, as expected from
# classic anharmonicity theory.
# Wu and Wentzcovitch (2009) proposed a value of c=0.1 for periclase,
# but as mentioned above, the model implemented here is somewhat simplified
# and c_anh needs to be adjusted downwards to avoid excessive depression of
# the isochoric heat capacity, Grueneisen parameter, and bulk modulus.
# Note also that the quasiharmonic reference state is somewhat different to
# the one used by Wu and Wentzcovitch (2009); in particular, the thermal
# expansion in the SLB model is significantly larger at high temperatures:
# at 0 GPa and 3000 K alpha ~ 8e-5 / K in Wu and Wentzcovitch (2009),
# but in the SLB model, alpha ~ 21.6e-5 / K at the same conditions.
params["c_anh"] = 0.04

# Create the new mineral object with the property modifiers from the SLB model
m = Mineral(params, property_modifiers=m_SLB.property_modifiers)

# And we're done! Now we can find the properties of the new mineral object
# using set_state or evaluate.
m.set_state(5.0e9, 4000.0)

# Print some properties to standard output
print(
    f"Properties of the new {m.name} object at "
    f"{m.pressure/1.e9:.1f} GPa and {m.temperature:.1f} K:"
)
print(f"Volume: {m.V * 1.e6:.1f} cm^3/mol")
print(f"Entropy: {m.S:.1f} J/K/mol")
print(f"Isochoric heat capacity: {m.molar_heat_capacity_v:.1f} J/K/mol")
print(f"Isobaric heat capacity: {m.molar_heat_capacity_p:.1f} J/K/mol")
print(f"Thermal expansion coefficient: {m.alpha*1.e6:.1f} * 1/MK")
print(f"Grueneisen parameter: {m.gr:.2f}")
print(f"Isothermal bulk modulus: {m.K_T / 1.e9:.1f} GPa")
print(f"Isentropic bulk modulus: {m.K_S / 1.e9:.1f} GPa")

# Create the figure and axes for the plots
fig = plt.figure(figsize=(12, 6))
ax = [fig.add_subplot(2, 4, i) for i in range(1, 8)]

# Create plots for the mineral properties as a function of temperature
# at different pressures.
temperatures = np.linspace(10.0, 3000.0, 101)
for P_GPa in [0, 10, 30, 60, 100]:
    P = P_GPa * 1.0e9
    pressures = temperatures * 0.0 + P

    volumes_SLB, Cv_SLB, Cp_SLB, alpha_SLB, gr_SLB, K_T_SLB = m_SLB.evaluate(
        [
            "V",
            "molar_heat_capacity_v",
            "molar_heat_capacity_p",
            "alpha",
            "gr",
            "K_T",
        ],
        pressures,
        temperatures,
    )
    ln = ax[0].plot(temperatures, volumes_SLB / m_SLB.params["V_0"], linestyle=":")
    ax[1].plot(temperatures, Cv_SLB, linestyle=":", color=ln[0].get_color())
    ax[2].plot(temperatures, Cp_SLB, linestyle=":", color=ln[0].get_color())
    ax[3].plot(temperatures, alpha_SLB * 1.0e5, linestyle=":", color=ln[0].get_color())
    ax[4].plot(temperatures, gr_SLB, linestyle=":", color=ln[0].get_color())
    ax[5].plot(temperatures, K_T_SLB / 1.0e9, linestyle=":", color=ln[0].get_color())
    ax[6].plot(
        temperatures,
        alpha_SLB * K_T_SLB / 1.0e6,
        linestyle=":",
        color=ln[0].get_color(),
    )

    volumes, Cv, Cp, alpha, gr, K_T = m.evaluate(
        [
            "V",
            "molar_heat_capacity_v",
            "molar_heat_capacity_p",
            "alpha",
            "gr",
            "K_T",
        ],
        pressures,
        temperatures,
    )
    ax[0].plot(temperatures, volumes / m.params["V_0"], c=ln[0].get_color())
    ax[1].plot(temperatures, Cv, c=ln[0].get_color(), label=f"{P/1.e9:.0f} GPa")
    ax[2].plot(temperatures, Cp, c=ln[0].get_color())
    ax[3].plot(temperatures, alpha * 1.0e5, c=ln[0].get_color())
    ax[4].plot(temperatures, gr, c=ln[0].get_color())
    ax[5].plot(temperatures, K_T / 1.0e9, c=ln[0].get_color())
    ax[6].plot(temperatures, alpha * K_T / 1.0e6, c=ln[0].get_color())

ax[1].legend()

for i in range(6):
    ax[i].set_xlabel("Temperature (K)")
    ax[i].set_xlim(0.0, 3000.0)
ax[1].set_ylim(0.0, 80.0)
ax[2].set_ylim(0.0, 80.0)
ax[3].set_ylim(0.0, 10.0)
ax[4].set_ylim(1.0, 2.0)
ax[5].set_ylim(100.0, 300.0)

ax[0].set_ylabel("$V/V_0$ (cm³/mol)")
ax[1].set_ylabel("$C_V$ (J/mol/K)")
ax[2].set_ylabel("$C_P$ (J/mol/K)")
ax[3].set_ylabel("$\\alpha$ (10$^{{-5}}$/K)")
ax[4].set_ylabel("$\\gamma$")
ax[5].set_ylabel("$K_T$ (GPa)")
ax[6].set_ylabel("$\\alpha K_T$ (MPa/K)")

fig.set_layout_engine("tight")
plt.show()
