import unittest
from util import BurnManTest
import warnings
import numpy as np

import burnman
from burnman import minerals
import burnman.classes.averaging_schemes as avg
from burnman.tools.eos import check_eos_consistency


class mypericlase(burnman.Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """

    def __init__(self):
        self.params = {
            "equation_of_state": "slb3",
            "V_0": 11.24e-6,
            "K_0": 161.0e9,
            "Kprime_0": 3.8,
            "G_0": 131.0e9,
            "Gprime_0": 2.1,
            "molar_mass": 0.0403,
            "n": 2,
            "Debye_0": 773.0,
            "grueneisen_0": 1.5,
            "q_0": 1.5,
            "eta_s_0": 2.8,
        }
        burnman.Mineral.__init__(self)


class my_nonrigid_mineral(burnman.Mineral):
    def __init__(self):
        self.params = {
            "equation_of_state": "slb3",
            "V_0": 11.24e-6,
            "K_0": 161.0e9,
            "Kprime_0": 3.8,
            "G_0": 0.0e9,
            "Gprime_0": 0.0,
            "molar_mass": 0.0403,
            "n": 2,
            "Debye_0": 773.0,
            "grueneisen_0": 1.5,
            "q_0": 1.5,
            "eta_s_0": 0.0,
        }
        burnman.Mineral.__init__(self)


class VRH_average(BurnManTest):
    def test_one_object(self):
        v = avg.voigt_reuss_hill_function([2.0], [0.123])
        self.assertFloatEqual(0.123, v)

    def test_two_same(self):
        v = avg.voigt_reuss_hill_function([2.0, 2.0], [0.456, 0.456])
        self.assertFloatEqual(0.456, v)

    def test_one_no_volume(self):
        v = avg.voigt_reuss_hill_function([0.0, 2.0], [0.123, 0.456])
        self.assertFloatEqual(0.456, v)

    def test_mix(self):
        v = avg.voigt_reuss_hill_function([1.0, 2.0], [0.1, 0.2])
        self.assertFloatEqual(0.15833333333333, v)


class VRH(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_method("slb3")
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K_vrh[0] / 1.0e9)
        self.assertFloatEqual(150.901, G_vrh[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)

        rock.set_method("slb3")
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K_vrh[0] / 1.0e9)
        self.assertFloatEqual(150.901, G_vrh[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_method("slb3")
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(9956.708, v_p[0])
        self.assertFloatEqual(5606.378, v_s[0])
        self.assertFloatEqual(7564.880, v_phi[0])
        self.assertFloatEqual(279.408, K_vrh[0] / 1.0e9)
        self.assertFloatEqual(153.461, G_vrh[0] / 1.0e9)


class SelfConsistent(BurnManTest):
    def average(self, volumes, bulk, shear):
        scheme = avg.SelfConsistent()
        return np.array(
            [
                scheme.average_bulk_moduli(volumes, bulk, shear),
                scheme.average_shear_moduli(volumes, bulk, shear),
            ]
        )

    def test_single_and_identical_phases(self):
        for bulk, shear in [(160.0, 100.0), (10.0, 0.0), (0.0, 0.0)]:
            for volumes in [[2.0], [1.0, 2.0, 3.0]]:
                with self.subTest(bulk=bulk, shear=shear, volumes=volumes):
                    self.assertArraysAlmostEqual(
                        self.average(
                            volumes, [bulk] * len(volumes), [shear] * len(volumes)
                        ),
                        [bulk, shear],
                        tol=1.0e-12,
                    )

    def test_absent_phase(self):
        self.assertArraysAlmostEqual(
            self.average([0.0, 2.0, 0.0], [np.nan, 160.0, -1.0], [np.nan, 100.0, -1.0]),
            [160.0, 100.0],
        )

    def test_analytic_two_phase_solution(self):
        # At K_i = 4 G_i / 3, zeta = G and the equal-volume solution is
        # G = sqrt(G_1 G_2), K = 4 G / 3 (not the VRH mean).
        self.assertArraysAlmostEqual(
            self.average([1.0, 1.0], [4.0 / 3.0, 16.0 / 3.0], [1.0, 4.0]),
            [8.0 / 3.0, 2.0],
            tol=1.0e-12,
        )

    def test_common_shear_modulus(self):
        bulk = np.array([80.0, 240.0])
        shift = 4.0 * 60.0 / 3.0
        expected_bulk = 1.0 / np.sum(np.array([0.3, 0.7]) / (bulk + shift)) - shift
        self.assertArraysAlmostEqual(
            self.average([3.0, 7.0], bulk, [60.0, 60.0]),
            [expected_bulk, 60.0],
            tol=1.0e-12,
        )

    def test_multiphase_equations_and_bounds(self):
        volumes = np.array([2.0, 3.0, 5.0])
        fractions = volumes / sum(volumes)
        bulk = np.array([160.0, 80.0, 240.0])
        shear = np.array([100.0, 30.0, 60.0])
        k, g = self.average(volumes, bulk, shear)
        zeta = g * (9.0 * k + 8.0 * g) / (6.0 * (k + 2.0 * g))
        self.assertFloatEqual(
            np.sum(fractions * (bulk - k) / (bulk + 4.0 * g / 3.0)),
            0.0,
            tol_zero=1.0e-12,
        )
        self.assertFloatEqual(
            np.sum(fractions * (shear - g) / (shear + zeta)), 0.0, tol_zero=1.0e-12
        )
        lower = avg.HashinShtrikmanLower()
        upper = avg.HashinShtrikmanUpper()
        self.assertLessEqual(lower.average_bulk_moduli(volumes, bulk, shear), k)
        self.assertLessEqual(k, upper.average_bulk_moduli(volumes, bulk, shear))
        self.assertLessEqual(lower.average_shear_moduli(volumes, bulk, shear), g)
        self.assertLessEqual(g, upper.average_shear_moduli(volumes, bulk, shear))

    def test_phase_order_units_and_volume_scaling(self):
        volumes = np.array([2.0, 3.0, 5.0])
        bulk = np.array([160.0, 80.0, 240.0])
        shear = np.array([100.0, 30.0, 60.0])
        expected = self.average(volumes, bulk, shear)
        self.assertArraysAlmostEqual(
            self.average(volumes[::-1], bulk[::-1], shear[::-1]), expected, tol=1.0e-12
        )
        for volume_scale in [1.0e-100, 1.0e100]:
            for modulus_scale in [1.0e-9, 1.0e9]:
                with self.subTest(
                    volume_scale=volume_scale, modulus_scale=modulus_scale
                ):
                    actual = self.average(
                        volumes * volume_scale,
                        bulk * modulus_scale,
                        shear * modulus_scale,
                    )
                    self.assertArraysAlmostEqual(
                        actual / modulus_scale, expected, tol=1.0e-12
                    )

    def test_published_echoes_values(self):
        # Independent printed SC outputs in Echoes chapter 12, for spherical
        # pores with K_p = G_p = 1e-6 GPa in a K_s=72, G_s=32 GPa solid.
        # https://jfbarthelemy.github.io/echoes/4_schemes/schemes.html
        for porosity, reference, tolerance in [
            (0.2, [37.156, 19.629], 5.1e-4),  # printed to three decimals
            (0.3, [22.689245, 13.264364], 2.0e-5),  # iterative reference
        ]:
            with self.subTest(porosity=porosity):
                actual = self.average(
                    [1.0 - porosity, porosity], [72.0, 1.0e-6], [32.0, 1.0e-6]
                )
                for value, target in zip(actual, reference):
                    self.assertAlmostEqual(value, target, delta=tolerance)

    def test_all_fluids(self):
        self.assertArraysAlmostEqual(
            self.average([1.0, 3.0], [10.0, 20.0], [0.0, 0.0]),
            [16.0, 0.0],
            tol=1.0e-12,
        )

    def test_void_rigidity_threshold(self):
        for porosity in [0.49, 0.4999, 0.5, 0.7, 1.0]:
            with self.subTest(porosity=porosity):
                actual = self.average(
                    [1.0 - porosity, porosity], [72.0, 0.0], [32.0, 0.0]
                )
                if porosity < 0.5:
                    self.assertTrue(np.all(actual > 0.0))
                else:
                    self.assertArraysAlmostEqual(actual, [0.0, 0.0])

    def test_fluid_rigidity_threshold(self):
        for fraction in [0.59, 0.6, 0.8, 1.0]:
            with self.subTest(fraction=fraction):
                k, g = self.average(
                    [1.0 - fraction, fraction], [72.0, 2.0], [32.0, 0.0]
                )
                if fraction < 0.6:
                    self.assertGreater(g, 0.0)
                else:
                    self.assertFloatEqual(g, 0.0)
                    self.assertFloatEqual(
                        k, 1.0 / ((1.0 - fraction) / 72.0 + fraction / 2.0)
                    )

    def test_soft_pores_retain_positive_stiffness(self):
        k, g = self.average([0.4, 0.6], [72.0, 1.0e-6], [32.0, 1.0e-6])
        self.assertGreater(k, 1.0e-6)
        self.assertGreater(g, 1.0e-6)

    def test_invalid_inputs(self):
        for volumes, bulk, shear in [
            ([], [], []),
            ([1.0], [1.0, 2.0], [1.0]),
            ([[1.0]], [[1.0]], [[1.0]]),
            ([0.0, 0.0], [1.0, 2.0], [1.0, 2.0]),
            ([-1.0, 2.0], [1.0, 2.0], [1.0, 2.0]),
            ([np.nan], [1.0], [1.0]),
            ([np.inf], [1.0], [1.0]),
            ([1.0], [-1.0], [1.0]),
            ([1.0], [1.0], [-1.0]),
            ([1.0], [np.nan], [1.0]),
            ([1.0], [1.0], [np.inf]),
        ]:
            with self.subTest(volumes=volumes, bulk=bulk, shear=shear):
                with self.assertRaises(ValueError):
                    self.average(volumes, bulk, shear)

    def test_single_mineral_composite(self):
        mineral = minerals.SLB_2024.periclase()
        rock = burnman.Composite([mineral], [1.0])
        rock.set_averaging_scheme(avg.SelfConsistent())
        rock.set_state(10.0e9, 1500.0)
        self.assertArraysAlmostEqual(
            [rock.K_eff, rock.G_eff, rock.v_p, rock.v_s, rock.v_phi],
            [mineral.K_S, mineral.G, mineral.v_p, mineral.v_s, mineral.v_phi],
            tol=1.0e-12,
        )
        self.assertTrue(check_eos_consistency(mineral, P=10.0e9, T=1500.0))
        self.assertTrue(
            check_eos_consistency(
                rock, P=10.0e9, T=1500.0, including_shear_properties=False
            )
        )

    def test_slb_solution_composite(self):
        phases = [
            minerals.SLB_2024.bridgmanite([0.8, 0.1, 0.1, 0.0, 0.0, 0.0, 0.0]),
            minerals.SLB_2024.periclase(),
        ]
        rock = burnman.Composite(phases, [0.7, 0.3])
        rock.set_averaging_scheme("SelfConsistent")
        for pressure, temperature in [(10.0e9, 1500.0), (120.0e9, 2500.0)]:
            with self.subTest(P=pressure, T=temperature):
                rho, vp, vs, k, g = rock.evaluate(
                    ["rho", "v_p", "v_s", "K_eff", "G_eff"], [pressure], [temperature]
                )[:, 0]
                rock.set_state(pressure, temperature)
                expected = self.average(
                    rock.phase_volumes, [p.K_S for p in phases], [p.G for p in phases]
                )
                self.assertArraysAlmostEqual([k, g], expected, tol=1.0e-12)
                self.assertFloatEqual(
                    vp, np.sqrt((k + 4.0 * g / 3.0) / rho), tol=1.0e-12
                )
                self.assertFloatEqual(vs, np.sqrt(g / rho), tol=1.0e-12)
                # The EOS tester checks thermodynamic Reuss derivatives;
                # the seismic identities using effective moduli are above.
                self.assertTrue(
                    check_eos_consistency(
                        rock,
                        P=pressure,
                        T=temperature,
                        including_shear_properties=False,
                    )
                )
        self.assertMaterialAmountScaling(
            rock, extra_properties=("K_eff", "G_eff", "v_p", "v_s")
        )


class Reuss(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(9886.675, v_p[0])
        self.assertFloatEqual(5606.207, v_s[0])
        self.assertFloatEqual(7472.635, v_phi[0])
        self.assertFloatEqual(272.635, K[0] / 1.0e9)
        self.assertFloatEqual(153.452, G[0] / 1.0e9)

    def test_non_present_non_rigid_phase(self):
        rock = burnman.Composite([mypericlase(), my_nonrigid_mineral()], [1.0, 0.0])
        rock.set_averaging_scheme(avg.Reuss())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rho, v_p, v_s, v_phi, K, G = rock.evaluate(
                [
                    "rho",
                    "v_p",
                    "v_s",
                    "v_phi",
                    "K_eff",
                    "G_eff",
                ],
                10.0e9,
                300.0,
            )
            assert len(w) == 2  # we expect two errors to be thrown when p_nonrigid = 0
        self.assertFloatEqual(150.901, G / 1.0e9)

    def test_present_non_rigid_phase(self):
        rock = burnman.Composite([mypericlase(), my_nonrigid_mineral()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.Reuss())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rho, v_p, v_s, v_phi, K, G = rock.evaluate(
                [
                    "rho",
                    "v_p",
                    "v_s",
                    "v_phi",
                    "K_eff",
                    "G_eff",
                ],
                10.0e9,
                300.0,
            )
            assert (
                len(w) == 4
            )  # we expect four errors to be thrown when p_nonrigid != 0
        self.assertFloatEqual(0.0, G / 1.0e9)


class Voigt(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(10026.252, v_p[0])
        self.assertFloatEqual(5606.549, v_s[0])
        self.assertFloatEqual(7656.014, v_phi[0])
        self.assertFloatEqual(286.180, K[0] / 1.0e9)
        self.assertFloatEqual(153.471, G[0] / 1.0e9)


class HSLower(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(9951.001, v_p[0])
        self.assertFloatEqual(5606.377, v_s[0])
        self.assertFloatEqual(7557.367, v_phi[0])
        self.assertFloatEqual(278.853, K[0] / 1.0e9)
        self.assertFloatEqual(153.461, G[0] / 1.0e9)


class HSUpper(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(9951.842, v_p[0])
        self.assertFloatEqual(5606.387, v_s[0])
        self.assertFloatEqual(7558.466, v_phi[0])
        self.assertFloatEqual(278.934, K[0] / 1.0e9)
        self.assertFloatEqual(153.462, G[0] / 1.0e9)


class HSAverage(BurnManTest):
    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def same(self, number):
        rock = burnman.Composite([mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.0e9)
        self.assertFloatEqual(150.901, G[0] / 1.0e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()],
            [0.5, 0.5],
        )
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = rock.evaluate(
            ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"],
            [10.0e9],
            [300.0],
        )
        self.assertFloatEqual(4882.407, rho[0])
        self.assertFloatEqual(9951.421, v_p[0])
        self.assertFloatEqual(5606.382, v_s[0])
        self.assertFloatEqual(7557.917, v_phi[0])
        self.assertFloatEqual(278.893, K[0] / 1.0e9)
        self.assertFloatEqual(153.461, G[0] / 1.0e9)


if __name__ == "__main__":
    unittest.main()
