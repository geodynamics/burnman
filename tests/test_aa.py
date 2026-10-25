# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences, released under the GNU GPL v2 or later.

import copy
from util import BurnManTest

import numpy as np
from scipy.integrate import quad

import burnman
from burnman.tools.eos import check_eos_consistency


class AndersonAhrens(BurnManTest):
    def setUp(self):
        self.iron = burnman.minerals.other.liquid_iron()
        self.params = self.iron.params
        self.eos = self.iron.method

    def test_reference_state(self):
        p, eos = self.params, self.eos
        T, V, P = p["T_0"], p["V_0"], p["P_0"]
        self.assertEqual(eos.pressure(T, V, p), P)
        self.assertEqual(eos.entropy(P, T, V, p), p["S_0"])
        self.assertEqual(eos.molar_internal_energy(P, T, V, p), p["E_0"])
        self.assertEqual(eos.helmholtz_energy(P, T, V, p), p["E_0"] - T * p["S_0"])
        self.assertFloatEqual(
            eos.grueneisen_parameter(P, T, V, p), p["grueneisen_0"], tol=2e-15
        )
        self.assertFloatEqual(
            eos.isentropic_bulk_modulus_reuss(P, T, V, p), p["K_S"], tol=2e-15
        )
        calculated_volume = eos.volume(P, T, p)
        self.assertFloatEqual(calculated_volume, V, tol=2e-14)
        self.assertLess(abs(eos.pressure(T, calculated_volume, p) - P), 0.01)

    def test_reference_pressure_work_and_temperature(self):
        p = dict(self.params, P_0=2.5e9)
        zero = dict(p, P_0=0.0)
        for ratio in [0.55, 0.8, 1.0, 1.3]:
            with self.subTest(ratio=ratio):
                V = ratio * p["V_0"]
                work = self.eos._isentropic_energy_change(
                    V, p
                ) - self.eos._isentropic_energy_change(V, zero)
                self.assertFloatEqual(
                    work, p["P_0"] * (p["V_0"] - V), tol_zero=1e-9, tol=1e-13
                )
                Ts = self.eos._reference_temperature(V, p)
                self.assertFloatEqual(
                    self.eos.pressure(Ts, V, p),
                    self.eos._reference_pressure(V, p),
                    tol_zero=1e-4,
                    tol=1e-13,
                )
                self.assertFloatEqual(
                    self.eos.isentropic_bulk_modulus_reuss(0.0, Ts, V, p),
                    self.eos._reference_bulk_modulus(V, p),
                    tol=1e-13,
                )
                self.assertEqual(self.eos.entropy(0.0, Ts, V, p), p["S_0"])
                dl = 1e-5
                slope = (
                    np.log(self.eos._reference_temperature(V * np.exp(-dl), p))
                    - np.log(self.eos._reference_temperature(V * np.exp(dl), p))
                ) / (2 * dl)
                gamma = p["grueneisen_0"] + p["grueneisen_prime"] * ratio ** (
                    -p["grueneisen_n"]
                ) * self.eos._isentropic_energy_change(V, p)
                self.assertFloatEqual(slope, gamma, tol=1e-9)

    def test_reference_temperature_singular_exponents(self):
        for exponent in [0.0, -2.0 / 3.0, -4.0 / 3.0, -2.0, -8.0 / 3.0, 1.0]:
            for ratio in [0.6, 1.0, 1.3]:
                with self.subTest(exponent=exponent, ratio=ratio):
                    p = dict(self.params, grueneisen_n=exponent, P_0=2.5e9)
                    ell = -np.log(ratio)
                    integral = quad(
                        lambda z: p["grueneisen_0"]
                        + p["grueneisen_prime"]
                        * np.exp(exponent * z)
                        * self.eos._isentropic_energy_change(p["V_0"] * np.exp(-z), p),
                        0.0,
                        ell,
                        epsabs=1e-12,
                    )[0]
                    self.assertFloatEqual(
                        self.eos._reference_temperature(ratio * p["V_0"], p),
                        p["T_0"] * np.exp(integral),
                        tol=1e-13,
                    )

    def test_small_strain_energy(self):
        p = self.params
        ell = 1e-10
        V = p["V_0"] * np.exp(-ell)
        expected = p["V_0"] * (p["P_0"] * ell + 0.5 * p["K_S"] * ell**2)
        self.assertFloatEqual(
            self.eos._isentropic_energy_change(V, p), expected, tol=2e-6, tol_zero=0.0
        )

    def test_heat_capacity_and_figure5_expression_are_unchanged(self):
        p, eos = self.params, self.eos
        for rho in [5000.0, 7019.0, 10000.0, 15000.0]:
            for T in [1000.0, 1811.0, 5000.0, 15000.0]:
                with self.subTest(rho=rho, T=T):
                    V = p["molar_mass"] / rho
                    r = p["V_0"] / V
                    A = p["a"][0] + p["a"][1] / r
                    B = p["b"][0] + p["b"][1] / r**2
                    Theta = p["Theta"][0] * r ** p["Theta"][1]
                    xi = p["xi_0"] * r**-0.6
                    logistic = 1.0 / (1.0 + np.exp((r - p["F"][0]) / p["F"][1]))
                    lmda = logistic * (
                        p["lmda"][0] + p["lmda"][1] * r + p["lmda"][2] * r**0.4
                    )
                    expected = (
                        1.5 * p["n"] * burnman.constants.gas_constant
                        + A * T * T / (Theta * Theta + T * T)
                        + B * T**0.6
                        + (lmda * T + xi * p["theta"]) / (p["theta"] + T)
                    )
                    self.assertFloatEqual(
                        eos.molar_heat_capacity_v(0.0, T, V, p), expected, tol=2e-15
                    )

    def test_zero_reference_pressure_retains_previous_free_energy(self):
        # Recorded from the original AA implementation. Setting P_0=0 removes
        # the intended reference-pressure correction from the comparison.
        p = dict(self.params, P_0=0.0)
        cases = [
            (
                5000.0,
                5000.0,
                -630205.9208554473,
                166.79387460706806,
                37.874064801504474,
            ),
            (7019.0, 1811.0, -181100.00000000003, 100.0, 35.9580348274425),
            (
                10000.0,
                7000.0,
                -678896.7891424899,
                126.93412115230728,
                37.84536603348719,
            ),
            (
                13000.0,
                7000.0,
                -411867.4803151762,
                110.90353154347144,
                33.90557395400822,
            ),
            (
                15000.0,
                15000.0,
                -1118472.526226906,
                129.97493939581932,
                36.582954003664966,
            ),
        ]
        for rho, T, F, S, Cv in cases:
            with self.subTest(rho=rho, T=T):
                V = p["molar_mass"] / rho
                # SciPy releases expose slightly different precision for R.
                # Normalize only that known kinetic contribution in the fixture.
                dCkin = (
                    1.5 * p["n"] * (burnman.constants.gas_constant - 8.31446261815324)
                )
                Ts = self.eos._reference_temperature(V, p)
                log_ratio = np.log(T / Ts)
                F += dCkin * (T - Ts - T * log_ratio)
                S += dCkin * log_ratio
                Cv += dCkin
                self.assertFloatEqual(
                    self.eos._helmholtz_energy(T, V, p), F, tol=5e-13, tol_zero=1e-7
                )
                self.assertFloatEqual(self.eos.entropy(0.0, T, V, p), S, tol=5e-13)
                self.assertFloatEqual(
                    self.eos.molar_heat_capacity_v(0.0, T, V, p), Cv, tol=2e-15
                )

    def test_caloric_integrals(self):
        p, eos = self.params, self.eos
        for rho, T in [(5000.0, 7000.0), (10000.0, 5000.0), (15000.0, 4000.0)]:
            with self.subTest(rho=rho, T=T):
                V = p["molar_mass"] / rho
                Ts = eos._reference_temperature(V, p)
                expected_e = (
                    p["E_0"]
                    + eos._isentropic_energy_change(V, p)
                    + quad(lambda t: eos.molar_heat_capacity_v(0.0, t, V, p), Ts, T)[0]
                )
                expected_s = (
                    p["S_0"]
                    + quad(
                        lambda t: eos.molar_heat_capacity_v(0.0, t, V, p) / t, Ts, T
                    )[0]
                )
                self.assertFloatEqual(
                    eos.molar_internal_energy(0.0, T, V, p),
                    expected_e,
                    tol=2e-13,
                    tol_zero=1e-8,
                )
                self.assertFloatEqual(eos.entropy(0.0, T, V, p), expected_s, tol=2e-13)

    def test_analytic_derivatives_obey_first_law_and_maxwell_relation(self):
        p, eos = self.params, self.eos
        for rho, T in [
            (7019.0, 1811.0),
            (5000.0, 7000.0),
            (10000.0, 5000.0),
            (13000.0, 7000.0),
        ]:
            with self.subTest(rho=rho, T=T):
                V = p["molar_mass"] / rho
                dV, dT = 1e-5 * V, 0.05
                P, KT, aKT = eos._thermal_derivatives(T, V, p)
                f = lambda v, t: eos.helmholtz_energy(0.0, t, v, p)
                s = lambda v, t: eos.entropy(0.0, t, v, p)
                energy = lambda v, t: eos.molar_internal_energy(0.0, t, v, p)
                self.assertFloatEqual(
                    -(f(V + dV, T) - f(V - dV, T)) / (2 * dV),
                    P,
                    tol=5e-8,
                    tol_zero=100.0,
                )
                self.assertFloatEqual(
                    -(f(V, T + dT) - f(V, T - dT)) / (2 * dT), s(V, T), tol=1e-8
                )
                self.assertFloatEqual(
                    (energy(V, T + dT) - energy(V, T - dT)) / (2 * dT),
                    eos.molar_heat_capacity_v(0.0, T, V, p),
                    tol=1e-8,
                )
                self.assertFloatEqual(
                    (eos.pressure(T + dT, V, p) - eos.pressure(T - dT, V, p))
                    / (2 * dT),
                    aKT,
                    tol=1e-8,
                )
                self.assertFloatEqual(
                    (s(V + dV, T) - s(V - dV, T)) / (2 * dV), aKT, tol=1e-8
                )
                self.assertFloatEqual(
                    -V
                    * (eos.pressure(T, V + dV, p) - eos.pressure(T, V - dV, p))
                    / (2 * dV),
                    KT,
                    tol=1e-8,
                )
                self.assertFloatEqual(
                    (energy(V + dV, T) - energy(V - dV, T)) / (2 * dV),
                    T * aKT - P,
                    tol=1e-8,
                    tol_zero=100.0,
                )

    def test_derivatives_with_perturbed_coefficients(self):
        p = copy.deepcopy(self.params)
        p.update(
            a=[10.0, 5.0],
            b=[0.03, -0.02],
            Theta=[2300.0, 0.7],
            lmda=[-2.0, 4.0, 3.0],
            xi_0=20.0,
            F=[1.2, 0.3],
            theta=4000.0,
            grueneisen_n=-0.8,
        )
        T, V = 6500.0, 0.7 * p["V_0"]
        dV = 1e-5 * V
        eos = self.eos
        derivative = -(
            eos._helmholtz_energy(T, V + dV, p) - eos._helmholtz_energy(T, V - dV, p)
        ) / (2 * dV)
        self.assertFloatEqual(eos.pressure(T, V, p), derivative, tol=1e-8)
        derivative = (
            -V * (eos.pressure(T, V + dV, p) - eos.pressure(T, V - dV, p)) / (2 * dV)
        )
        self.assertFloatEqual(
            eos.isothermal_bulk_modulus_reuss(0.0, T, V, p), derivative, tol=1e-8
        )

    def test_volume_accuracy_is_independent_of_formula_unit_size(self):
        for scale in [1e-3, 1.0, 1e3]:
            p = copy.deepcopy(self.params)
            for key in ["V_0", "E_0", "S_0", "molar_mass", "n", "xi_0"]:
                p[key] *= scale
            for key in ["a", "b", "lmda"]:
                p[key] = [value * scale for value in p[key]]
            p["grueneisen_prime"] /= scale
            for rho, T in [(7019.0, 1811.0), (10000.0, 5000.0), (13000.0, 7000.0)]:
                with self.subTest(scale=scale, rho=rho, T=T):
                    V = p["molar_mass"] / rho
                    P = self.eos.pressure(T, V, p)
                    calculated = self.eos.volume(P, T, p)
                    self.assertFloatEqual(calculated, V, tol=2e-14)
                    self.assertLess(abs(self.eos.pressure(T, calculated, p) - P), 0.02)

    def test_derivatives_do_not_depend_on_energy_or_entropy_gauge(self):
        p = dict(self.params, E_0=1e12, S_0=1e9)
        T, V = 7000.0, 0.7 * p["V_0"]
        self.assertEqual(
            self.eos._thermal_derivatives(T, V, p),
            self.eos._thermal_derivatives(T, V, self.params),
        )
        expected = (p["E_0"] - self.params["E_0"]) - T * (p["S_0"] - self.params["S_0"])
        self.assertFloatEqual(
            self.eos._helmholtz_energy(T, V, p)
            - self.eos._helmholtz_energy(T, V, self.params),
            expected,
            tol=1e-15,
        )

    def test_gibbs_uses_supplied_volume(self):
        # Deliberately supply V independently of P: re-solving V would change PV.
        p = self.params
        P, T, V = 130e9, 5000.0, 0.8 * p["V_0"]
        self.assertFloatEqual(
            self.eos.gibbs_energy(P, T, V, p),
            self.eos._helmholtz_energy(T, V, p) + P * V,
            tol=1e-14,
        )

    def test_all_public_properties_agree_with_mineral(self):
        eos, p = self.eos, self.params
        P, T = 130e9, 5000.0
        self.iron.set_state(P, T)
        V = self.iron.V
        properties = {
            "helmholtz_energy": "molar_helmholtz",
            "gibbs_energy": "molar_gibbs",
            "molar_internal_energy": "molar_internal_energy",
            "enthalpy": "molar_enthalpy",
            "entropy": "molar_entropy",
            "molar_heat_capacity_v": "molar_heat_capacity_v",
            "molar_heat_capacity_p": "molar_heat_capacity_p",
            "grueneisen_parameter": "grueneisen_parameter",
            "isothermal_bulk_modulus_reuss": "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss": "isentropic_bulk_modulus_reuss",
            "thermal_expansivity": "thermal_expansivity",
            "shear_modulus": "shear_modulus",
        }
        for method, attribute in properties.items():
            with self.subTest(method=method):
                self.assertFloatEqual(
                    getattr(eos, method)(P, T, V, p),
                    getattr(self.iron, attribute),
                    tol=2e-14,
                    tol_zero=1e-8,
                )

    def test_preserves_heat_capacity_based_grueneisen_trend(self):
        p = self.params
        T, V, dt = p["T_0"], p["V_0"], 0.1
        slope = (
            self.eos.grueneisen_parameter(0.0, T + dt, V, p)
            - self.eos.grueneisen_parameter(0.0, T - dt, V, p)
        ) / (2 * dt)
        self.assertFloatEqual(slope, 1.0773006668962637e-4, tol=1e-8)
        self.assertGreater(slope, 0.0)

    def test_vectorized_properties(self):
        p = self.params
        T = np.array([2000.0, 5000.0, 7000.0])
        V = p["molar_mass"] / np.array([5000.0, 10000.0, 13000.0])
        for method in ["pressure", "_helmholtz_energy", "_molar_heat_capacity_v"]:
            f = getattr(self.eos, method)
            self.assertArraysAlmostEqual(
                f(T, V, p), [f(t, v, p) for t, v in zip(T, V)], tol=1e-14
            )
        for method in [
            "entropy",
            "isothermal_bulk_modulus_reuss",
            "molar_heat_capacity_p",
            "thermal_expansivity",
        ]:
            f = getattr(self.eos, method)
            self.assertArraysAlmostEqual(
                f(0.0, T, V, p), [f(0.0, t, v, p) for t, v in zip(T, V)], tol=1e-14
            )

    def test_burnman_consistency_checker(self):
        for P, T in [(1e5, 1811.0), (130e9, 5000.0), (320e9, 7000.0)]:
            with self.subTest(P=P, T=T):
                self.assertTrue(
                    check_eos_consistency(
                        self.iron, P=P, T=T, tol=1e-5, including_shear_properties=False
                    )
                )

    def test_parameter_validation(self):
        for key in [
            "P_0",
            "T_0",
            "S_0",
            "E_0",
            "V_0",
            "K_S",
            "Kprime_S",
            "Kprime_prime_S",
            "molar_mass",
            "n",
            "grueneisen_0",
            "grueneisen_prime",
            "grueneisen_n",
            "theta",
            "xi_0",
            "a",
            "b",
            "Theta",
            "lmda",
            "F",
        ]:
            with self.subTest(missing=key):
                p = dict(self.params)
                del p[key]
                with self.assertRaisesRegex(KeyError, key):
                    self.eos.validate_parameters(p)
        for key, value in [
            ("V_0", 0.0),
            ("T_0", -1.0),
            ("K_S", 0.0),
            ("molar_mass", 0.0),
            ("n", 0.0),
            ("theta", 0.0),
            ("grueneisen_prime", np.nan),
            ("Kprime_S", np.inf),
            ("a", [1.0]),
            ("b", [1.0, np.nan]),
            ("Theta", [0.0, 1.0]),
            ("F", [1.0, 0.0]),
            ("lmda", [1.0, 2.0]),
        ]:
            with self.subTest(key=key, value=value):
                with self.assertRaises(ValueError):
                    self.eos.validate_parameters(dict(self.params, **{key: value}))

    def test_invalid_states(self):
        p = self.params
        for T, V in [
            (0.0, p["V_0"]),
            (-1.0, p["V_0"]),
            (np.nan, p["V_0"]),
            (1811.0, 0.0),
            (1811.0, np.inf),
        ]:
            for method in ["pressure", "_helmholtz_energy", "_molar_heat_capacity_v"]:
                with self.subTest(T=T, V=V, method=method):
                    with self.assertRaises(ValueError):
                        getattr(self.eos, method)(T, V, p)
        with self.assertRaises(ValueError):
            self.eos.volume(np.inf, 1811.0, p)

    def test_logistic_does_not_overflow(self):
        p = self.params
        with np.errstate(over="raise"):
            lmda, xi = self.eos._lambdaxi(p["V_0"] / 1000.0, p)
        self.assertEqual(lmda, 0.0)
        self.assertGreater(xi, 0.0)
