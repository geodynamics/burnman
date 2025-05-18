# Benchmarks for the solid solution class

import numpy as np
import matplotlib.pyplot as plt

from burnman.tools.chemistry import equilibrium_temperature
from burnman.tools.eos import check_eos_consistency
from burnman.minerals import HGP_2018_ds633


# Here, we test that burnman outputs the correct properties for
# albite liquid from the Holland and Powell dataset (v. 6.33, 2018)
ab = HGP_2018_ds633.ab()
abL = HGP_2018_ds633.abL()

# Output from tc350beta1-linux using tc-ds633
P, T, abL_gibbs_thermocalc = np.array(
    [
        [1.0e5, 273.15 + 25.0, -3970.1027e3],
        [1.0e5, 273.15 + 1000.0, -4428.8506e3],
        [1.0e8, 273.15 + 25.0, -3959.2855e3],
        [1.0e8, 273.15 + 1000.0, -4417.6771e3],
        [1.0e10, 273.15 + 25.0, -2992.7298e3],
        [1.0e10, 273.15 + 1000.0, -3426.9106e3],
    ]
).T

abL_delta_gibbs_bm_tc = abL.evaluate(["gibbs"], P, T)[0] - abL_gibbs_thermocalc

print("abL Gibbs energy difference (burnman - thermocalc v3.50, J/mol):")
for i in range(len(P)):
    print(f"{P[i]/1.e9} GPa, {T[i]} K: {abL_delta_gibbs_bm_tc[i]:0.0f}")
print("")

# Check that the implemented equation of state has self-consistent properties.
check_eos_consistency(
    abL, P=1.0e10, T=2000.0, verbose=True, including_shear_properties=False
)

# Calculate S, V, G for solid and liquid albite,
# and plot the volumes.
temperatures = np.linspace(300.0, 2000.0, 101)
for pressure in [1.0e5, 1.0e9]:
    pressures = pressure + 0.0 * temperatures

    S_liq, V_liq, gibbs_liq = abL.evaluate(["S", "V", "gibbs"], pressures, temperatures)
    S_sol, V_sol, gibbs_sol = ab.evaluate(["S", "V", "gibbs"], pressures, temperatures)

    plt.plot(
        temperatures,
        V_liq * 1.0e6,
        linestyle=":",
        label=f"liquid albite ({pressure/1.e9} GPa)",
    )
    plt.plot(temperatures, V_sol * 1.0e6, label=f"solid albite ({pressure/1.e9} GPa)")

plt.xlabel("Temperature (K)")
plt.ylabel("Volume (cm^3/mol)")
plt.legend()

plt.show()

# Calculate the melting curve of albite to 1.8 GPa
# (for comparison with Figure 3 in Holland and Powell, 1998)
pressures = np.linspace(1.0e5, 1.8e9, 101)
temperatures = np.empty_like(pressures)

for i, P in enumerate(pressures):
    temperatures[i] = equilibrium_temperature([ab, abL], [1.0, -1.0], P, 1500.0)

plt.plot(pressures / 1.0e9, temperatures - 273.15)
plt.xlabel("Pressure (GPa)")
plt.ylabel("Temperature (C)")
plt.show()
