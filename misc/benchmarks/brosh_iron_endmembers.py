# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

# Benchmarks for the solid solution class


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import burnman
from burnman import minerals
from burnman import equilibrate
from burnman.tools.eos import check_eos_consistency


fcc = minerals.SE_2015.fcc_iron()
bcc = minerals.SE_2015.bcc_iron()
hcp = minerals.SE_2015.hcp_iron()
liq = minerals.SE_2015.liquid_iron()


calib = [
    [bcc, 1.0e5, 1.0, 69538.0, 7.0496],
    [bcc, 1.0e5, 300.0, -8183.0, 7.0947],
    [bcc, 1.0e9, 1.0, 76567.0, 7.0094],
    [bcc, 1.0e5, 1000.0, -42271.0, 7.3119],
    [bcc, 1.0e9, 1000.0, -34987.0, 7.2583],
    [fcc, 3.0e9, 1000.0, -20604.0, 7.0250],
    [fcc, 10.0e9, 1000.0, 27419.0, 6.7139],
    [hcp, 40.0e9, 1000.0, 214601.0, 5.8583],
]

for (m, P, T, G, V) in calib:
    m.set_state(P, T)
    # assert that the Gibbs energies are within 20 J/mol of each other
    assert m.gibbs - G < 20.0
    # assert that the volumes are within 0.01 % of each other
    assert (m.V * 1.0e6 - V) / V < 1.0e-4
    print(f"{P/1.e9}, {T}, {m.gibbs - G:.1f}, {m.V*1.e6 - V:.1e}")

check_eos_consistency(
    liq, P=1.0e10, T=2000.0, tol=1.0e-4, verbose=True, including_shear_properties=False
)

pressures = np.linspace(100.0e9, 400.0e9, 101)
temperatures = [300.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0]
for T in temperatures:
    Ts = T * np.ones_like(pressures)

    plt.plot(
        pressures / 1.0e9, liq.evaluate(["rho"], pressures, Ts)[0], label=f"{int(T)} K"
    )

plt.legend()
plt.xlabel("Pressure (GPa)")
plt.ylabel("Density (kg/m$^3$)")
plt.show()

# from SE15ver.dat bundled with PerpleX 6.8.3 (September 2018)
Fe_diag_img = mpimg.imread("figures/fe_perplex.png")
plt.imshow(Fe_diag_img, extent=[0.0, 350.0, 300, 8000], alpha=0.3, aspect="auto")

# alternative, from Brosh's 2007 paper.
# Slightly cropped to be the same scale as Saxena and Eriksson (2015)
# Fe_diag_img_Brosh = mpimg.imread('figures/fe_brosh.png')
# plt.imshow(Fe_diag_img_Brosh, extent=[0.0, 350.0, 300, 8000],
#            alpha=0.3, aspect='auto')


def invariant(m1, m2, m3, P=5.0e9, T=2000.0):
    composition = m1.formula
    assemblage = burnman.Composite([m1, m2, m3])
    assemblage.set_state(P, T)
    equality_constraints = [
        ("phase_fraction", (m1, 0.0)),
        ("phase_fraction", (m2, 0.0)),
    ]
    sol, prm = equilibrate(
        composition, assemblage, equality_constraints, store_iterates=False
    )
    return sol.x[0:2]


def univariant(m1, m2, condition_constraints, P=5.0e9, T=2000.0):
    composition = m1.formula
    assemblage = burnman.Composite([m1, m2])
    assemblage.set_state(P, T)
    equality_constraints = [condition_constraints, ("phase_fraction", (m1, 0.0))]
    sols, prm = equilibrate(
        composition, assemblage, equality_constraints, store_iterates=False
    )

    pressures = np.array([s.x[0] for s in sols])
    temperatures = np.array([s.x[1] for s in sols])
    return pressures, temperatures


PT = univariant(fcc, liq, ("P", np.array([21.0e9, 21.0e9])), P=21.0e9, T=2500.0)

print(f"FCC-liq equilibrium at {PT[0][0]/1.e9} GPa: {int(PT[1][0])} K")

Tmin = 1.0
Pmin = 1.0e5
Pmax = 350.0e9

# BCC-FCC-LIQ invariant
Pinv, Tinv = invariant(bcc, fcc, liq, P=5.0e9, T=2000.0)

pressures, temperatures = univariant(
    bcc, liq, ("P", np.linspace(Pmin, Pinv, 11)), P=Pmin, T=1800.0
)
plt.plot(pressures / 1.0e9, temperatures)

pressures, temperatures = univariant(
    bcc, fcc, ("P", np.linspace(Pmin, Pinv, 11)), P=Pmin, T=1800.0
)
plt.plot(pressures / 1.0e9, temperatures, label="bcc-fcc")


# FCC-HCP-LIQ invariance
Pinv2, Tinv2 = invariant(fcc, hcp, liq, P=90.0e9, T=3000.0)

pressures, temperatures = univariant(
    fcc, liq, ("P", np.linspace(Pinv, Pinv2, 11)), P=Pinv, T=Tinv
)
plt.plot(pressures / 1.0e9, temperatures, label="fcc-liq")

pressures, temperatures = univariant(
    hcp, liq, ("P", np.linspace(Pinv2, Pmax, 11)), P=Pinv2, T=Tinv2
)
plt.plot(pressures / 1.0e9, temperatures, label="hcp-liq")


# BCC-FCC-HCP invariance
Pinv3, Tinv3 = invariant(bcc, fcc, hcp, P=10.0e9, T=1000.0)

pressures, temperatures = univariant(
    bcc, hcp, ("T", np.linspace(Tmin, Tinv3, 11)), P=10.0e9, T=Tmin
)
plt.plot(pressures / 1.0e9, temperatures, label="bcc-hcp")

pressures, temperatures = univariant(
    bcc, fcc, ("P", np.linspace(Pmin, Pinv3, 11)), P=Pmin, T=Tinv3
)
plt.plot(pressures / 1.0e9, temperatures, label="bcc-fcc")

pressures, temperatures = univariant(
    fcc, hcp, ("P", np.linspace(Pinv3, Pinv2, 11)), P=Pinv3, T=Tinv3
)
plt.plot(pressures / 1.0e9, temperatures, label="fcc-hcp")


plt.scatter([21.0], [2516.0], label="Saxena and Eriksson, 2015, Fig2")

plt.xlim(0.0, 350.0)
plt.ylim(0.0, 8000.0)
plt.xlabel("Pressure (GPa)")
plt.ylabel("Temperature (K)")
plt.legend()
plt.show()
