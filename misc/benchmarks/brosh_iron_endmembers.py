# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

# Benchmarks for the Brosh-Calphad Equation of State
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

import burnman
from burnman import minerals
from burnman import equilibrate
from burnman.tools.eos import check_eos_consistency

fcc = minerals.SE_2015.fcc_iron()
bcc = minerals.SE_2015.bcc_iron()
hcp = minerals.SE_2015.hcp_iron()
liq = minerals.SE_2015.liquid_iron()


# PerpleX comparisons within the SGTE temperature range. The former 1 K
# Gibbs values included an unphysical extrapolation of the SGTE polynomial.
calib = [
    [bcc, 1.0e5, 300.0, -8183.0, 7.0947],
    [bcc, 1.0e5, 1000.0, -42271.0, 7.3119],
    [bcc, 1.0e9, 1000.0, -34987.0, 7.2583],
    [fcc, 3.0e9, 1000.0, -20604.0, 7.0250],
    [fcc, 10.0e9, 1000.0, 27419.0, 6.7139],
    [hcp, 40.0e9, 1000.0, 214601.0, 5.8583],
]

for m, P, T, G, V in calib:
    m.set_state(P, T)
    # assert that the Gibbs energies are within 20 J/mol of each other
    assert abs(m.gibbs - G) < 20.0
    # Volumes agree within 0.01%, allowing for reference rounding to 4 decimals.
    assert abs(m.V * 1.0e6 - V) < 1.0e-4 * V + 0.5e-4
    print(f"{P/1.e9}, {T}, {m.gibbs - G:.1f}, {m.V*1.e6 - V:.1e}")

check_eos_consistency(
    liq, P=1.0e10, T=2000.0, tol=1.0e-4, verbose=True, including_shear_properties=False
)


# Saxena and Eriksson (2015), Section 4, versus the Table 1 hcp variant.
# These published quantities are not mutually reproduced by one parameter set.
def melting_temperature(solid, pressure):
    def delta_g(T):
        solid.set_state(pressure, T)
        liq.set_state(pressure, T)
        return solid.gibbs - liq.gibbs

    return brentq(delta_g, 1000.0, 9000.0)


paper_hcp = minerals.SE_2015.hcp_iron()
paper_hcp.params["b"] = [0.7, 5.5]
P_core, T_core = 363.85e9, 6453.0
print(
    "Core benchmark: published Tm=6453 K, rho=13.71 g/cm^3; "
    "inner-core averages Cp=29.45 J/mol/K, S=92 J/mol/K, gamma=1.81"
)
for label, solid in [("BurnMan", hcp), ("Table 1 hcp b1", paper_hcp)]:
    Tm = melting_temperature(solid, P_core)
    solid.set_state(P_core, T_core)
    print(
        f"{label}: Tm={Tm:.2f} K, rho={solid.rho/1000.:.5f} g/cm^3, "
        f"Cp={solid.C_p:.5f} J/mol/K, S={solid.S:.5f} J/mol/K, "
        f"gamma={solid.gr:.5f}"
    )
    if solid is hcp:
        assert abs(Tm - T_core) < 5.0
        assert abs(solid.rho / 1000.0 - 13.71) < 0.01

# Table 2 gives approximate ambient-pressure polynomials without a fit range.
# Print residuals rather than asserting exact agreement with these fits.
for solid, a, k in [
    (liq, [3.7270e-5, 5.8520e-9, -0.63660e-12], [1122415.4, -284.63, 0.0229]),
    (hcp, [5.8690e-5, 19.690e-9, -2.2180e-12], [1568000.0, -263.3, 0.0202]),
]:
    for T in [1000.0, 2000.0, 4000.0]:
        solid.set_state(1.0e5, T)
        alpha_fit = a[0] + a[1] * T + a[2] * T**2
        K_fit = (k[0] + k[1] * T + k[2] * T**2) * 1.0e5
        print(
            f"Table 2 {solid.name}, {T:.0f} K: "
            f"alpha={solid.alpha:.7g}, fit={alpha_fit:.7g}, "
            f"K={solid.K_T/1.e9:.4f} GPa, fit={K_fit/1.e9:.4f} GPa"
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

# Perple_X gmet2 uses only the n=4 cold term,
# whereas BurnMan combines n=2..5 (Brosh 2007, Appendix A).
# Source: https://github.com/jadconnolly/Perple_X/blob/main/src/tlib.f
# Values below are hcp-liquid coexistence roots calculated with FRENDLY
# 7.1.13 (26 March 2025) and its bundled SE15ver.dat, in GPa and K.
perplex_melting = np.array(
    [
        [100.0, 3579.16],
        [150.0, 4237.35],
        [200.0, 4808.54],
        [250.0, 5324.94],
        [300.0, 5801.92],
        [350.0, 6247.77],
    ]
)

# Isolate the effect of that approximation without changing the default
# minerals or introducing another EoS strategy.
c4_hcp = minerals.SE_2015.hcp_iron()
c4_liq = minerals.SE_2015.liquid_iron()
for phase in (c4_hcp, c4_liq):
    phase.params["c"] = [0.0, 0.0, 1.0, 0.0]

c4_temperatures = []
for P_GPa, reference_T in perplex_melting:
    P = P_GPa * 1.0e9

    def delta_g(T):
        c4_hcp.set_state(P, T)
        c4_liq.set_state(P, T)
        return c4_hcp.gibbs - c4_liq.gibbs

    T_c4 = brentq(delta_g, 3000.0, 7500.0)
    T_full = melting_temperature(hcp, P)
    c4_temperatures.append(T_c4)
    print(
        f"Melting at {P_GPa:.0f} GPa: full Brosh={T_full:.2f} K, "
        f"n=4 cold curve={T_c4:.2f} K, Perple_X 7.1.13={reference_T:.2f} K"
    )
    # The remaining 2-3 K difference comes from Perple_X's reference
    # derivatives: its sqrt(T) terms omit their coefficient in Sr and Cp_r.
    assert abs(T_c4 - reference_T) < 5.0

plt.plot(
    perplex_melting[:, 0],
    c4_temperatures,
    "k--",
    label="n=4 cold curve, same thermal terms",
)
plt.plot(
    perplex_melting[:, 0],
    perplex_melting[:, 1],
    "ko",
    fillstyle="none",
    label="Perple_X 7.1.13 numerical reference",
)


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
plt.plot(pressures / 1.0e9, temperatures, color="C0", label="bcc-liq")

pressures, temperatures = univariant(
    bcc, fcc, ("P", np.linspace(Pmin, Pinv, 11)), P=Pmin, T=1800.0
)
plt.plot(pressures / 1.0e9, temperatures, color="C1", label="bcc-fcc")


# FCC-HCP-LIQ invariance
Pinv2, Tinv2 = invariant(fcc, hcp, liq, P=90.0e9, T=3000.0)

pressures, temperatures = univariant(
    fcc, liq, ("P", np.linspace(Pinv, Pinv2, 11)), P=Pinv, T=Tinv
)
plt.plot(pressures / 1.0e9, temperatures, color="C2", label="fcc-liq")

pressures, temperatures = univariant(
    hcp, liq, ("P", np.linspace(Pinv2, Pmax, 11)), P=Pinv2, T=Tinv2
)
plt.plot(pressures / 1.0e9, temperatures, color="C3", label="hcp-liq")


# BCC-FCC-HCP invariance
Pinv3, Tinv3 = invariant(bcc, fcc, hcp, P=10.0e9, T=1000.0)

pressures, temperatures = univariant(
    bcc, hcp, ("T", np.linspace(Tmin, Tinv3, 11)), P=10.0e9, T=Tmin
)
plt.plot(pressures / 1.0e9, temperatures, color="C4", label="bcc-hcp")

pressures, temperatures = univariant(
    bcc, fcc, ("P", np.linspace(Pmin, Pinv3, 11)), P=Pmin, T=Tinv3
)
plt.plot(pressures / 1.0e9, temperatures, color="C1")

pressures, temperatures = univariant(
    fcc, hcp, ("P", np.linspace(Pinv3, Pinv2, 11)), P=Pinv3, T=Tinv3
)
plt.plot(pressures / 1.0e9, temperatures, color="C5", label="fcc-hcp")


plt.xlim(0.0, 350.0)
plt.ylim(0.0, 8000.0)
plt.xlabel("Pressure (GPa)")
plt.ylabel("Temperature (K)")
plt.title("Melting: full Brosh versus the Perple_X approximation")
plt.legend(loc="upper left", fontsize="small")
plt.show()
