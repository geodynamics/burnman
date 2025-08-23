# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_modular_mgd
-------------------

This example shows how to create a mineral object using
the Modular Mie-Grüneisen-Debye equation of state.

This equation of state allows users to pick and choose
different options for two key aspects of the mineral's behavior:

1. The equation of state (EOS) used to describe the mineral's thermodynamic properties at the reference temperature.
2. The model used to describe the mineral's Debye temperature as a function of volume.

and a further two options if the user wants to account for anharmonic effects:

1. The model used to describe the mineral's anharmonicity as a function of temperature.
2. The model used to describe the prefactor as a function of volume.

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

# For this example, we will use the SLB_2011 fayalite mineral
# as a starting point
from burnman.minerals.SLB_2011 import fayalite as mineral

# This helper function creates an Equation of State object
from burnman.eos.helper import create

# These imports contain the different models available for the EOS
from burnman.eos import debye_temperature_models
from burnman.eos import anharmonic_prefactor_models
from burnman.eos import anharmonic_thermal_models

# Create the SLB mineral object
m_SLB = mineral()

# Now create a new mineral object using the modular MGD EOS
params = m_SLB.params.copy()
params["equation_of_state"] = "modular_mgd"

# In these two lines, we set the reference EOS and Debye temperature model
# to be an exact clone of the original SLB object
params["reference_eos"] = create("bm3")
params["debye_temperature_model"] = debye_temperature_models.SLB()

# We are also able to add an anharmonicity model
# that is not a part of the SLB model. This model has two parts.
# The first is a volume-dependent prefactor model,
# which is here chosen to follow a power law with parameters
# a_anh and m_anh: i.e, A(V) = a_anh * (V/V_0)^m_anh
params["anharmonic_prefactor_model"] = anharmonic_prefactor_models.PowerLaw()
params["a_anh"] = 1.0
params["m_anh"] = 10.0

# The second part is a temperature-dependent anharmonic model,
# which is here chosen so that the modification to the isochoric heat capacity
# is based on the CDF log-normal distribution with parameters mu_anh and sigma_anh.
params["anharmonic_thermal_model"] = anharmonic_thermal_models.LogNormal()
params["mu_anh"] = 0.0
params["sigma_anh"] = 1.0

# Create the new mineral object with the property modifiers from the SLB model
m = Mineral(params, property_modifiers=m_SLB.property_modifiers)

# And we're done! Now we can find the properties of the new mineral object
# using set_state or evaluate.
m.set_state(5.0e9, 4000.0)

# Print some properties to standard output
print(
    "Properties of the new mineral object at "
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
temperatures = np.linspace(10.0, 5000.0, 101)
for P in np.linspace(0.0, 30.0e9, 7):
    pressures = temperatures * 0.0 + P
    try:
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
        ax[3].plot(temperatures, alpha_SLB, linestyle=":", color=ln[0].get_color())
        ax[4].plot(temperatures, gr_SLB, linestyle=":", color=ln[0].get_color())
        ax[5].plot(
            temperatures, K_T_SLB / 1.0e9, linestyle=":", color=ln[0].get_color()
        )
        ax[6].plot(
            temperatures,
            alpha_SLB * K_T_SLB / 1.0e6,
            linestyle=":",
            color=ln[0].get_color(),
        )
    except Exception:
        temperatures_lower = np.linspace(10.0, 2695.0 + P / 7.0e6, 101)
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
            temperatures_lower,
        )
        ln = ax[0].plot(
            temperatures_lower, volumes_SLB / m_SLB.params["V_0"], linestyle=":"
        )
        ax[1].plot(temperatures_lower, Cv_SLB, linestyle=":", color=ln[0].get_color())
        ax[2].plot(temperatures_lower, Cp_SLB, linestyle=":", color=ln[0].get_color())
        ax[3].plot(
            temperatures_lower, alpha_SLB, linestyle=":", color=ln[0].get_color()
        )
        ax[4].plot(temperatures_lower, gr_SLB, linestyle=":", color=ln[0].get_color())
        ax[5].plot(
            temperatures_lower, K_T_SLB / 1.0e9, linestyle=":", color=ln[0].get_color()
        )
        ax[6].plot(
            temperatures_lower,
            alpha_SLB * K_T_SLB / 1.0e6,
            linestyle=":",
            color=ln[0].get_color(),
        )
    try:
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
        ax[3].plot(temperatures, alpha, c=ln[0].get_color())
        ax[4].plot(temperatures, gr, c=ln[0].get_color())
        ax[5].plot(temperatures, K_T / 1.0e9, c=ln[0].get_color())
        ax[6].plot(temperatures, alpha * K_T / 1.0e6, c=ln[0].get_color())
    except Exception:
        temperatures_lower = np.linspace(10.0, 2695.0 + P / 1.0e7, 101)
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
            temperatures_lower,
        )
        ax[0].plot(temperatures_lower, volumes / m.params["V_0"], c=ln[0].get_color())
        ax[1].plot(
            temperatures_lower, Cv, c=ln[0].get_color(), label=f"{P/1.e9:.0f} GPa"
        )
        ax[2].plot(temperatures_lower, Cp, c=ln[0].get_color())
        ax[3].plot(temperatures_lower, alpha, c=ln[0].get_color())
        ax[4].plot(temperatures_lower, gr, c=ln[0].get_color())
        ax[5].plot(temperatures_lower, K_T / 1.0e9, c=ln[0].get_color())
        ax[6].plot(temperatures_lower, alpha * K_T / 1.0e6, c=ln[0].get_color())

ax[1].legend()

for i in range(6):
    ax[i].set_xlabel("Temperature (K)")
    ax[i].set_xlim(0.0, 5000.0)
ax[1].set_ylim(0.0, 180.0)
ax[2].set_ylim(0.0, 250.0)
ax[3].set_ylim(0.0, 3.0e-4)

ax[0].set_ylabel("$V/V_0$ (cm³/mol)")
ax[1].set_ylabel("$C_V$ (J/mol/K)")
ax[2].set_ylabel("$C_P$ (J/mol/K)")
ax[3].set_ylabel("$\\alpha$ (1/K)")
ax[4].set_ylabel("$\\gamma$")
ax[5].set_ylabel("$K_T$ (GPa)")
ax[6].set_ylabel("$\\alpha K_T$ (MPa/K)")

fig.set_layout_engine("tight")
plt.show()


# Now do exactly the same thing again, but using evaluate_with_volumes
fig = plt.figure(figsize=(12, 6))
ax = [fig.add_subplot(2, 4, i) for i in range(1, 8)]

temperatures = np.linspace(10.0, 5000.0, 101)
for Vrels in np.linspace(1.2, 0.8, 5):
    volumes = temperatures * 0.0 + Vrels * m_SLB.params["V_0"]

    try:
        pressures_SLB, Cv_SLB, Cp_SLB, alpha_SLB, gr_SLB, K_T_SLB = (
            m_SLB.evaluate_with_volumes(
                [
                    "pressure",
                    "molar_heat_capacity_v",
                    "molar_heat_capacity_p",
                    "alpha",
                    "gr",
                    "K_T",
                ],
                volumes,
                temperatures,
            )
        )
        ln = ax[0].plot(temperatures, pressures_SLB / 1.0e9, linestyle=":")
        ax[1].plot(temperatures, Cv_SLB, linestyle=":", color=ln[0].get_color())
        ax[2].plot(temperatures, Cp_SLB, linestyle=":", color=ln[0].get_color())
        ax[3].plot(temperatures, alpha_SLB, linestyle=":", color=ln[0].get_color())
        ax[4].plot(temperatures, gr_SLB, linestyle=":", color=ln[0].get_color())
        ax[5].plot(
            temperatures, K_T_SLB / 1.0e9, linestyle=":", color=ln[0].get_color()
        )
        ax[6].plot(
            temperatures,
            alpha_SLB * K_T_SLB / 1.0e6,
            linestyle=":",
            color=ln[0].get_color(),
        )
    except Exception:
        temperatures_lower = np.linspace(10.0, 2695.0 + (1.0 - Vrels) * 1000.0, 101)
        pressures_SLB, Cv_SLB, Cp_SLB, alpha_SLB, gr_SLB, K_T_SLB = (
            m_SLB.evaluate_with_volumes(
                [
                    "pressure",
                    "molar_heat_capacity_v",
                    "molar_heat_capacity_p",
                    "alpha",
                    "gr",
                    "K_T",
                ],
                volumes,
                temperatures_lower,
            )
        )
        ln = ax[0].plot(temperatures_lower, pressures_SLB / 1.0e9, linestyle=":")
        ax[1].plot(temperatures_lower, Cv_SLB, linestyle=":", color=ln[0].get_color())
        ax[2].plot(temperatures_lower, Cp_SLB, linestyle=":", color=ln[0].get_color())
        ax[3].plot(
            temperatures_lower, alpha_SLB, linestyle=":", color=ln[0].get_color()
        )
        ax[4].plot(temperatures_lower, gr_SLB, linestyle=":", color=ln[0].get_color())
        ax[5].plot(
            temperatures_lower, K_T_SLB / 1.0e9, linestyle=":", color=ln[0].get_color()
        )
        ax[6].plot(
            temperatures_lower,
            alpha_SLB * K_T_SLB / 1.0e6,
            linestyle=":",
            color=ln[0].get_color(),
        )
    try:
        pressures, Cv, Cp, alpha, gr, K_T = m.evaluate_with_volumes(
            [
                "pressure",
                "molar_heat_capacity_v",
                "molar_heat_capacity_p",
                "alpha",
                "gr",
                "K_T",
            ],
            volumes,
            temperatures,
        )
        ax[0].plot(temperatures, pressures / 1.0e9, c=ln[0].get_color(), linestyle="--")
        ax[1].plot(
            temperatures,
            Cv,
            c=ln[0].get_color(),
            linestyle="--",
            label=f"$V_{{rel}}$ = {Vrels:.1f}",
        )
        ax[2].plot(temperatures, Cp, c=ln[0].get_color(), linestyle="--")
        ax[3].plot(temperatures, alpha, c=ln[0].get_color(), linestyle="--")
        ax[4].plot(temperatures, gr, c=ln[0].get_color(), linestyle="--")
        ax[5].plot(temperatures, K_T / 1.0e9, c=ln[0].get_color(), linestyle="--")
        ax[6].plot(
            temperatures, alpha * K_T / 1.0e6, c=ln[0].get_color(), linestyle="--"
        )
    except Exception:
        temperatures_lower = np.linspace(10.0, 2695.0 + (1.0 - Vrels) * 1000.0, 101)
        pressures, Cv, Cp, alpha, gr, K_T = m.evaluate_with_volumes(
            [
                "pressure",
                "molar_heat_capacity_v",
                "molar_heat_capacity_p",
                "alpha",
                "gr",
                "K_T",
            ],
            volumes,
            temperatures_lower,
        )
        ax[0].plot(
            temperatures_lower, pressures / 1.0e9, c=ln[0].get_color(), linestyle="--"
        )
        ax[1].plot(
            temperatures_lower,
            Cv,
            c=ln[0].get_color(),
            linestyle="--",
            label=f"$V_{{rel}}$ = {Vrels:.1f}",
        )
        ax[2].plot(temperatures_lower, Cp, c=ln[0].get_color(), linestyle="--")
        ax[3].plot(temperatures_lower, alpha, c=ln[0].get_color(), linestyle="--")
        ax[4].plot(temperatures_lower, gr, c=ln[0].get_color(), linestyle="--")
        ax[5].plot(temperatures_lower, K_T / 1.0e9, c=ln[0].get_color(), linestyle="--")
        ax[6].plot(
            temperatures_lower, alpha * K_T / 1.0e6, c=ln[0].get_color(), linestyle="--"
        )

ax[1].legend()

for i in range(6):
    ax[i].set_xlabel("Temperature (K)")
    ax[i].set_xlim(0.0, 5000.0)
ax[2].set_ylim(0.0, 400.0)
ax[3].set_ylim(0.0, 3.0e-4)

ax[0].set_ylabel("$P$ (GPa)")
ax[1].set_ylabel("$C_V$ (J/mol/K)")
ax[2].set_ylabel("$C_P$ (J/mol/K)")
ax[3].set_ylabel("$\\alpha$ (1/K)")
ax[4].set_ylabel("$\\gamma$")
ax[5].set_ylabel("$K_T$ (GPa)")
ax[6].set_ylabel("$\\alpha K_T$ (MPa/K)")
fig.set_layout_engine("tight")
plt.show()
