# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences, released under the GNU GPL v2 or later.

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

import burnman
from util import BurnManTest


class BroshCalphad(BurnManTest):
    def setUp(self):
        self.phases = [
            getattr(burnman.minerals.SE_2015, name + "_iron")()
            for name in ("bcc", "fcc", "hcp", "liquid")
        ]

    def test_low_temperature_continuation(self):
        # Integrate Appendix C's heat capacity from the SGTE reference state.
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            P, Tr = 300.0e9, p["T_0"]
            theta = eos._theta(P, p)
            I = eos._interpolating_function(P, p)

            def cp(t):
                x = theta / t
                return (
                    3.0
                    * p["n"]
                    * burnman.constants.gas_constant
                    * x**2
                    * np.exp(-x)
                    / (1.0 - np.exp(-x)) ** 2
                    + I * t / Tr * p["delta_Cpr"]
                )

            Gr = eos.gibbs_energy(P, Tr, 0.0, p)
            Sr = eos.entropy(P, Tr, 0.0, p)
            for T in (1.0, 100.0):
                with self.subTest(phase=p["name"], T=T):
                    expected_s = Sr + quad(lambda t: cp(t) / t, Tr, T)[0]
                    expected_g = (
                        Gr
                        - Sr * (T - Tr)
                        - quad(lambda t: (T - t) * cp(t) / t, Tr, T)[0]
                    )
                    self.assertFloatEqual(eos.entropy(P, T, 0.0, p), expected_s)
                    self.assertFloatEqual(
                        eos.gibbs_energy(P, T, 0.0, p), expected_g, tol=1.0e-12
                    )
                    self.assertFloatEqual(
                        eos.molar_heat_capacity_p(P, T, 0.0, p), cp(T), tol=1.0e-13
                    )
                    self.assertGreater(eos.molar_heat_capacity_p(P, T, 0.0, p), 0.0)

    def test_zero_temperature_limit(self):
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            P = 300.0e9
            with self.subTest(phase=p["name"]), np.errstate(
                over="raise", invalid="raise", divide="raise"
            ):
                V = eos.volume(P, 0.0, p)
                self.assertTrue(np.isfinite(eos.gibbs_energy(P, 0.0, V, p)))
                self.assertTrue(np.isfinite(eos.entropy(P, 0.0, V, p)))
                self.assertGreater(eos.isothermal_bulk_modulus_reuss(P, 0.0, V, p), 0.0)
                self.assertEqual(eos.molar_heat_capacity_p(P, 0.0, V, p), 0.0)
                self.assertEqual(eos.thermal_expansivity(P, 0.0, V, p), 0.0)
                self.assertTrue(np.isfinite(eos.entropy(P, 0.01, V, p)))

    def test_last_sgte_segment_is_extrapolated(self):
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            last = dict(p, gibbs_coefficients=[p["gibbs_coefficients"][-1]])
            for T in (12000.0, 12000.01, 15000.0):
                with self.subTest(phase=p["name"], T=T):
                    for method in (eos._gibbs_1bar, eos._S_1bar, eos._Cp_1bar):
                        self.assertEqual(method(T, p), method(T, last))

    def test_reference_calorics(self):
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            for T in (p["T_0"], 1000.0, 3000.0, 8000.0):
                with self.subTest(phase=p["name"], T=T):
                    self.assertFloatEqual(
                        eos.gibbs_energy(0.0, T, 0.0, p),
                        eos._gibbs_1bar(T, p),
                        tol=1.0e-13,
                        tol_zero=1.0e-9,
                    )
                    self.assertFloatEqual(
                        eos.entropy(0.0, T, 0.0, p), eos._S_1bar(T, p)
                    )
                    self.assertFloatEqual(
                        eos.molar_heat_capacity_p(0.0, T, 0.0, p), eos._Cp_1bar(T, p)
                    )

    def test_thermodynamic_derivatives(self):
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            for T in (100.0, 300.0, 2000.0, 8000.0):
                with self.subTest(phase=p["name"], T=T):
                    P, dP, dT = 1.0e10, 1.0e7, 0.01
                    V = eos.volume(P, T, p)
                    g = lambda p_, t: eos.gibbs_energy(p_, t, 0.0, p)
                    s = lambda p_, t: eos.entropy(p_, t, 0.0, p)
                    v = lambda p_, t: eos.volume(p_, t, p)
                    self.assertFloatEqual(
                        -(g(P, T + dT) - g(P, T - dT)) / (2 * dT), s(P, T), tol=1.0e-7
                    )
                    self.assertFloatEqual(
                        (g(P + dP, T) - g(P - dP, T)) / (2 * dP), V, tol=1.0e-6
                    )
                    self.assertFloatEqual(
                        T * (s(P, T + dT) - s(P, T - dT)) / (2 * dT),
                        eos.molar_heat_capacity_p(P, T, V, p),
                        tol=1.0e-7,
                    )
                    self.assertFloatEqual(
                        (v(P, T + dT) - v(P, T - dT)) / (2 * dT) / V,
                        eos.thermal_expansivity(P, T, V, p),
                        tol=1.0e-7,
                    )
                    self.assertFloatEqual(
                        -(s(P + dP, T) - s(P - dP, T)) / (2 * dP),
                        eos.thermal_expansivity(P, T, V, p) * V,
                        tol=1.0e-6,
                    )
                    self.assertFloatEqual(
                        (v(P + dP, T) - v(P - dP, T)) / (2 * dP),
                        eos.dvolume_dP(P, T, p),
                        tol=1.0e-6,
                        tol_zero=1.0e-24,
                    )

    def test_derivatives_at_sgte_joins(self):
        # Use the selected segment's right derivative, without crossing a join.
        for mineral in self.phases:
            eos, p = mineral.method, mineral.params
            for T in (1811.0, 6000.0, 12000.0):
                with self.subTest(phase=p["name"], T=T):
                    P, dT = 300.0e9, 0.001
                    V = eos.volume(P, T, p)
                    dS = eos.entropy(P, T + dT, V, p) - eos.entropy(P, T, V, p)
                    dV = eos.volume(P, T + dT, p) - V
                    self.assertFloatEqual(
                        T * dS / dT, eos.molar_heat_capacity_p(P, T, V, p), tol=1.0e-6
                    )
                    self.assertFloatEqual(
                        dV / dT / V, eos.thermal_expansivity(P, T, V, p), tol=1.0e-6
                    )

    def test_se2015_core_melting(self):
        # Section 4: 6453 K and hcp density 13.71 g/cm^3 at 363.85 GPa.
        hcp, liquid = self.phases[2:]
        P = 363.85e9

        def delta_g(T):
            hcp.set_state(P, T)
            liquid.set_state(P, T)
            return hcp.gibbs - liquid.gibbs

        T = brentq(delta_g, 6000.0, 7000.0)
        self.assertLess(abs(T - 6453.0), 5.0)
        hcp.set_state(P, T)
        self.assertLess(abs(hcp.density / 1000.0 - 13.71), 0.01)
