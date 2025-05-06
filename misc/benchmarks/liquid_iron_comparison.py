# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.optimize import fsolve
import burnman
from burnman.tools.eos import check_eos_consistency


liq = burnman.minerals.other.liquid_iron()
check_eos_consistency(
    liq, P=10.0e9, T=7000.0, tol=1.0e-3, verbose=True, including_shear_properties=False
)

R = burnman.constants.gas_constant

# Find heat capacities
temperatures = np.linspace(1000.0, 15000.0, 101)
Cvs = np.empty_like(temperatures)
m = 0.055845
rhos = np.empty_like(temperatures)
densities = [5.0e3, 10.0e3, 15.0e3]
for rho in densities:
    V = m / rho
    for i, T in enumerate(temperatures):
        liq.set_state_with_volume(V, T)
        Cvs[i] = liq.molar_heat_capacity_v / R

    plt.plot(temperatures, Cvs)

fig1 = mpimg.imread(
    "../../burnman/data/input_figures/AA1994_liq_iron_TCv_different_densities.png"
)
plt.imshow(fig1, extent=[1000.0, 15000.0, 0.0, 6.0], aspect="auto")
plt.ylim(0.0, 6.0)
plt.title("AA1994, Figure 5")
plt.xlabel("Temperature (K)")
plt.ylabel("Cv/R/atom")
plt.show()


# Find volumes and temperatures up the reference isentrope
# Check standard state values
liq.set_state(1.0e5, 1811.0)
reference_entropy = liq.S


def isentrope(T, P, Sref, mineral):
    mineral.set_state(P, T[0])
    return Sref - mineral.S


pressures = np.linspace(0.0, 500.0e9, 21)
temperatures = np.empty_like(pressures)
rhos = np.empty_like(pressures)
Vps = np.empty_like(pressures)


for i, P in enumerate(pressures):
    temperatures[i] = fsolve(isentrope, 1811.0, args=(P, reference_entropy, liq))[0]
    rhos[i] = liq.density

fig1 = mpimg.imread(
    "../../burnman/data/input_figures/AA1994_liq_iron_PTrho_reference_isentrope.png"
)

plt.imshow(fig1, extent=[0.0, 500.0, 6.0, 15.0], aspect="auto")
plt.plot(pressures / 1.0e9, rhos / 1.0e3, marker="o", linestyle="None")
plt.title("1811 K isentrope; AA1994 Figure B1 (1/2)")
plt.xlabel("Pressure (GPa)")
plt.ylabel("Density (kg/m^3)")
plt.show()

plt.imshow(fig1, extent=[0.0, 500.0, 1500.0, 7000.0], aspect="auto")
plt.plot(pressures / 1.0e9, temperatures, marker="o", linestyle="None")
plt.title("1811 K isentrope; AA1994 Figure B1 (2/2)")
plt.xlabel("Pressure (GPa)")
plt.ylabel("Temperature (K)")
plt.show()


# Find densities at 1 bar
temperatures = np.linspace(1800.0, 2400.0, 100)
rhos = np.empty_like(temperatures)
rhos_mizuno = np.empty_like(temperatures)

P = 1.0e5
for i, T in enumerate(temperatures):
    liq.set_state(1.0e5, T)
    rhos[i] = liq.density / 1.0e3
    rhos_mizuno[i] = (7162.0 - 0.735 * (T - 1808.0)) / 1.0e3


fig1 = mpimg.imread("../../burnman/data/input_figures/AA1994_liq_iron_Trho_1bar.png")
plt.imshow(fig1, extent=[1800.0, 2400.0, 6.65, 7.1], aspect="auto")
plt.plot(temperatures, rhos, label="Model")
plt.plot(temperatures, rhos_mizuno, label="Mizuno et al.")
plt.ylim(6.65, 7.1)
plt.xlabel("Temperatures (K)")
plt.ylabel("Density (kg/m^3)")
plt.title("1 bar densities; AA1994 Figure 1")
plt.legend(loc="lower left")
plt.show()


def temperature(T, P, rho, mineral):
    mineral.set_state(P, T[0])
    return mineral.density - rho


# Check grueneisen values
Prhos = [
    [1.0e5, 7019.0],
    [0.2e9, 5500.0],
    [0.2e9, 6000.0],
    [0.2e9, 6500.0],
    [277.4e9, 12643.0],
    [331.5e9, 13015.0],
    [397.1e9, 13417.0],
]

print("Pressure (GPa), Temperature (K), Density (kg/m^3), Grueneisen parameter")
for Prho in Prhos:
    P, rho = Prho
    T = fsolve(temperature, 1811.0, args=(P, rho, liq))[0]
    print(
        "{0:.1f} {1:.0f} {2:.1f} {3:.4f}".format(
            P / 1.0e9, rho, T, liq.grueneisen_parameter
        )
    )
