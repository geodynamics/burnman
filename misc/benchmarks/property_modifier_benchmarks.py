from __future__ import absolute_import
from __future__ import print_function

# Benchmarks for the solid solution class
from burnman import minerals
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


# Quartz is a mineral with a Landau type order-disorder model
qtz = minerals.HGP_2018_ds633.q()

P = 1.0e5
temperatures = np.linspace(300.0, 1000.0, 101)
C_ps = np.empty_like(temperatures)
Qs = np.empty_like(temperatures)
Hs = np.empty_like(temperatures)
for i, T in enumerate(temperatures):
    qtz.set_state(P, T)
    Hs[i] = qtz._property_modifiers["G"] - T * qtz._property_modifiers["dGdT"]
    C_ps[i] = -T * qtz._property_modifiers["d2GdT2"]
    Qs[i] = qtz.property_modifier_properties[0]["Q"]

print("Quartz:")
print("T: {0}, Q: {1:.3f}".format(temperatures[0], Qs[0]))
print("T: {0}, Q: {1:.3f}".format(temperatures[-1], Qs[-1]))
print()

fig = plt.figure(figsize=(12.5, 4))
ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

ax[0].plot(temperatures, C_ps)
ax[1].plot(temperatures, Qs)
ax[2].plot(temperatures, Hs / 1e3)

for i, ylabel in enumerate(["$C_p$ (J/K/mol)", "$Q$", "$H$ (kJ/mol)"]):
    ax[i].set_xlim(temperatures[0], temperatures[-1])
    ax[i].set_xlabel("T (K)")
    ax[i].set_ylabel(ylabel)

plt.show()


# Spinel is a mineral with a Bragg-Williams type order-disorder model
sp = minerals.HGP_2018_ds633.sp()

P = 1.0e5
temperatures = np.linspace(400.0, 1500.0, 101)
C_ps = np.empty_like(temperatures)
Qs = np.empty_like(temperatures)
Hs = np.empty_like(temperatures)
for i, T in enumerate(temperatures):
    sp.set_state(P, T)
    Hs[i] = sp._property_modifiers["G"] - T * sp._property_modifiers["dGdT"]
    C_ps[i] = -T * sp._property_modifiers["d2GdT2"]
    Qs[i] = sp.property_modifier_properties[0]["Q"]

print("Spinel:")
print("T: {0}, Q: {1:.3f}".format(temperatures[0], Qs[0]))
print("T: {0}, Q: {1:.3f}".format(temperatures[-1], Qs[-1]))
print()

fig = plt.figure(figsize=(12.5, 4))
ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

ax[0].plot(temperatures, C_ps)
ax[1].plot(temperatures, Qs)
ax[2].plot(temperatures, Hs / 1e3)

for i, ylabel in enumerate(["$C_p$ (J/K/mol)", "$Q$", "$H$ (kJ/mol)"]):
    ax[i].set_xlim(temperatures[0], temperatures[-1])
    ax[i].set_xlabel("T (K)")
    ax[i].set_ylabel(ylabel)

plt.show()


# Albite is also a mineral with a Bragg-Williams type model
# Here we reproduce Figure 11 in HP1996
ab = minerals.HP_2011_ds62.ab()
ab.property_modifiers = [
    [
        "bragg_williams",
        {
            "deltaH": 14000.0,
            "deltaV": 4.2e-07,
            "Wh": 13600.0,
            "Wv": 4.2e-07,
            "n": 3.0,
            "factor": 1.0,
        },
    ]
]

P = 1.0e5
temperatures = np.linspace(400.0, 1200.0, 101)
C_ps = np.empty_like(temperatures)
Qs = np.empty_like(temperatures)
Hs = np.empty_like(temperatures)
for i, T in enumerate(temperatures):
    ab.set_state(P, T)
    Hs[i] = ab._property_modifiers["G"] - T * ab._property_modifiers["dGdT"]
    C_ps[i] = -T * ab._property_modifiers["d2GdT2"]
    Qs[i] = ab.property_modifier_properties[0]["Q"]

print("Albite:")
print("T: {0}, Q: {1:.3f}".format(temperatures[0], Qs[0]))
print("T: {0}, Q: {1:.3f}".format(temperatures[-1], Qs[-1]))


fig = plt.figure(figsize=(12.5, 4))
ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

img_ab_Cp = mpimg.imread("figures/albite_disordering_Cp.png")
ax[0].imshow(img_ab_Cp, extent=[400.0, 1100.0, 0.0, 120.0], aspect="auto")  # paper * 2
img_ab_enthalpy_Q = mpimg.imread("figures/albite_disordering_enthalpy_Q.png")
ax[1].imshow(img_ab_enthalpy_Q, extent=[450.0, 1250.0, 0.0, 1.0], aspect="auto")
ax[2].imshow(img_ab_enthalpy_Q, extent=[450.0, 1250.0, 0.0, 15.0], aspect="auto")

ax[0].plot(temperatures, C_ps)
ax[1].plot(temperatures, Qs)
ax[2].plot(temperatures, Hs / 1e3)

for i, ylabel in enumerate(["$C_p$ (J/K/mol)", "$Q$", "$H$ (kJ/mol)"]):
    ax[i].set_xlim(temperatures[0], temperatures[-1])
    ax[i].set_xlabel("T (K)")
    ax[i].set_ylabel(ylabel)

plt.show()
