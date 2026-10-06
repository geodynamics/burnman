import unittest
from unittest.mock import patch

import numpy as np

from burnman import constants
from burnman.minerals import HP_2011_fluids, Pitzer_Sterner_1994
from burnman.tools.eos import check_eos_consistency

# From burnman_cpp/tests/include/ps1994_reference.hpp, independently generated
# from Python BurnMan water_filter commit ffcb3688fdb4a28af9c76d5896d6df3abd55c7ba.
# The gas constant was set to PS1994's 8.31451 J/(mol K).
# Columns: P, T, V, K_T, alpha, delta G, delta S, delta Cp, stable roots.
# Deltas are relative to 1 bar at the same temperature; thermal references cancel.
REFERENCE_STATES = [
    (
        100000.0,
        500.0,
        0.04140096296511455,
        99585.33243479904,
        0.0020332222349105154,
        0.0,
        0.0,
        0.0,
        2.0,
    ),
    (
        2000000.0,
        500.0,
        0.0018893281618511445,
        1798947.6273069158,
        0.002879940900242607,
        12112.039244308515,
        -27.08617442352022,
        10.251662464684216,
        2.0,
    ),
    (
        3000000.0,
        500.0,
        2.1742755529782457e-05,
        812531476.3116649,
        0.0017149449851881025,
        13186.192436054698,
        -95.93106648597656,
        47.96702224554629,
        2.0,
    ),
    (
        5000000000.0,
        500.0,
        1.1970125398261703e-05,
        24135115061.75995,
        0.00022081822690930573,
        84510.1164014709,
        -125.20314274633537,
        17.556157233067793,
        1.0,
    ),
    (
        100000.0,
        1700.0,
        0.14135051001225282,
        100002.71441385848,
        0.0005883302708288478,
        0.0,
        0.0,
        0.0,
        1.0,
    ),
    (
        5000000000.0,
        1700.0,
        1.4955096679980533e-05,
        16153098004.150387,
        0.0001660301601020368,
        215068.692748687,
        -102.76236636600245,
        6.189852007898054,
        1.0,
    ),
    (
        1000.0,
        1200.0,
        9.977402634991858,
        999.9990614005583,
        0.0008333384995181995,
        -45946.75358937733,
        38.294062140618394,
        -0.01322468873671312,
        1.0,
    ),
    (
        1000000.0,
        600.0,
        0.004885973526491183,
        978914.7225092164,
        0.0018059726113527539,
        11395.461862664204,
        -19.590932861610156,
        1.8210710688759022,
        1.0,
    ),
    (
        4000000.0,
        550.0,
        0.0009931578143748743,
        3396069.0931882705,
        0.0030859387974489965,
        16327.04083421285,
        -33.860192358208394,
        15.39097841500638,
        2.0,
    ),
    (
        8000000.0,
        550.0,
        2.387014756137476e-05,
        461934004.6224184,
        0.0024365591737808734,
        17982.33911512245,
        -91.36344131590509,
        54.30493822868274,
        2.0,
    ),
    (
        50000000.0,
        550.0,
        2.233097918292227e-05,
        824267310.1939261,
        0.0016430201294631481,
        18949.11532622506,
        -93.2645639103325,
        44.7715640198037,
        1.0,
    ),
    (
        10000000.0,
        625.0,
        0.0004075318673552587,
        7389027.086086571,
        0.00394186871785531,
        22944.466329606134,
        -43.818569658303375,
        30.415693859053185,
        1.0,
    ),
    (
        20000000.0,
        625.0,
        3.0395255741777164e-05,
        89173580.81527053,
        0.007517011743528093,
        24880.782834310667,
        -83.55383255722433,
        112.15394215542966,
        1.0,
    ),
    (
        25000000.0,
        650.0,
        3.7124687926647476e-05,
        21158409.737423062,
        0.022194578942614885,
        27086.006809274026,
        -79.14575068338266,
        268.8648277125828,
        1.0,
    ),
    (
        100000000.0,
        700.0,
        2.7464741201411696e-05,
        373339944.53102344,
        0.002356205471032853,
        33468.531562076765,
        -84.2688254973965,
        51.050336359326955,
        1.0,
    ),
    (
        1000000000.0,
        1000.0,
        1.9852460985689748e-05,
        3675339487.373829,
        0.00045425154269857594,
        79964.40365391981,
        -91.4289127944962,
        17.430416910377673,
        1.0,
    ),
    (
        5000000000.0,
        1000.0,
        1.3252626645551434e-05,
        19857865943.43185,
        0.00018327196074269833,
        141811.987789983,
        -107.71104629511461,
        14.795657374228398,
        1.0,
    ),
    (
        5000000000.0,
        501.0,
        1.1972770490380761e-05,
        24122818895.626064,
        0.00022108002009724672,
        84635.3018598769,
        -125.16764679374523,
        17.973852160815838,
        1.0,
    ),
    (
        5000000000.0,
        1699.0,
        1.4952613794033352e-05,
        16157552668.094635,
        0.00016604286446667975,
        214965.92856115522,
        -102.76600931919154,
        6.192547370100115,
        1.0,
    ),
]


class PitzerSternerTests(unittest.TestCase):
    def test_reference_properties_and_stable_roots(self):
        with patch.object(constants, "gas_constant", 8.31451):
            for constructor in (
                HP_2011_fluids.H2O,
                Pitzer_Sterner_1994.H2O_Pitzer_Sterner,
            ):
                water, reference = constructor(), constructor()
                for p, t, v, bulk, alpha, dg, ds, dcp, roots in REFERENCE_STATES:
                    with self.subTest(reference=constructor.__name__, P=p, T=t):
                        water.set_state(p, t)
                        reference.set_state(1.0e5, t)
                        np.testing.assert_allclose(water.molar_volume, v, rtol=2.0e-12)
                        np.testing.assert_allclose(
                            water.isothermal_bulk_modulus_reuss, bulk, rtol=2.0e-8
                        )
                        np.testing.assert_allclose(
                            water.thermal_expansivity, alpha, rtol=2.0e-6
                        )
                        np.testing.assert_allclose(
                            water.molar_gibbs - reference.molar_gibbs,
                            dg,
                            rtol=2.0e-11,
                            atol=2.0e-8,
                        )
                        np.testing.assert_allclose(
                            water.molar_entropy - reference.molar_entropy,
                            ds,
                            rtol=2.0e-7,
                            atol=2.0e-7,
                        )
                        np.testing.assert_allclose(
                            water.molar_heat_capacity_p
                            - reference.molar_heat_capacity_p,
                            dcp,
                            rtol=2.0e-5,
                            atol=2.0e-3,
                        )
                        self.assertEqual(
                            water.method.shear_modulus(p, t, v, water.params), 0.0
                        )

    def test_thermocalc_absolute_water_energies(self):
        # Pure H2O states from tc350beta1 / tc-ds62.txt, extracted in
        # burnman_cpp/tests/include/thermocalc_reference.hpp. Expected values
        # come from the HPx-eos metapelite benchmark archive (2020-09-10).
        # Precision of the HP2011 reference enthalpy warrants 10 J/mol.
        water = HP_2011_fluids.H2O()
        for p, t, gibbs, enthalpy in (
            (2.5e8, 873.15, -367307.85, -240824.1),
            (1.0e9, 873.15, -351096.4, -234778.7),
            (1.0e8, 873.15, -372392.59, -236340.4),
            (1.0e9, 773.15, -338135.81, -240869.3),
        ):
            with self.subTest(P=p, T=t):
                water.set_state(p, t)
                self.assertAlmostEqual(water.molar_gibbs, gibbs, delta=10.0)
                self.assertAlmostEqual(water.molar_enthalpy, enthalpy, delta=10.0)

    def test_holland_powell_ideal_gas_reference(self):
        water = HP_2011_fluids.H2O()
        h, s = water.method._ideal_gas_reference(298.15, water.params)
        self.assertEqual(h, -241810.0)
        self.assertEqual(s, 188.80)
        # The reference component is evaluated at 298.15 K independently
        # of the fluid EOS's calibrated temperature range.
        t = 1000.0
        dt = 0.01
        hp, sp = water.method._ideal_gas_reference(t + dt, water.params)
        hm, sm = water.method._ideal_gas_reference(t - dt, water.params)
        expected_cp = 40.1 + 0.008656 * t + 487500.0 / t**2 - 251.2 / np.sqrt(t)
        np.testing.assert_allclose((hp - hm) / (2 * dt), expected_cp, rtol=1.0e-9)
        np.testing.assert_allclose(t * (sp - sm) / (2 * dt), expected_cp, rtol=1.0e-9)

    def test_thermodynamic_consistency(self):
        for constructor in (HP_2011_fluids.H2O, Pitzer_Sterner_1994.H2O_Pitzer_Sterner):
            with self.subTest(reference=constructor.__name__):
                self.assertTrue(
                    check_eos_consistency(
                        constructor(),
                        P=1.0e9,
                        T=1000.0,
                        tol=1.0e-5,
                        verbose=False,
                        including_shear_properties=False,
                    )
                )

    def test_water_domain(self):
        water = HP_2011_fluids.H2O()
        for pressure in (0.0, -1.0, np.nextafter(5.0e9, np.inf), np.inf, np.nan):
            with self.subTest(P=pressure):
                water.set_state(pressure, 1000.0)
                with self.assertRaises(ValueError):
                    _ = water.molar_volume
        for temperature in (
            np.nextafter(500.0, 0.0),
            np.nextafter(1700.0, np.inf),
            -1.0,
            np.inf,
            np.nan,
        ):
            with self.subTest(T=temperature):
                water.set_state(1.0e8, temperature)
                with self.assertRaises(ValueError):
                    _ = water.molar_gibbs
        for temperature in (500.0, 1700.0):
            water.set_state(5.0e9, temperature)
            self.assertTrue(np.isfinite(water.molar_gibbs))
            self.assertGreater(water.isothermal_bulk_modulus_reuss, 0.0)


if __name__ == "__main__":
    unittest.main()
