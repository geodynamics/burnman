"""Analytical radial structures, boundary conditions and query ordering."""

import unittest
import warnings
from util import BurnManTest

import numpy as np

import burnman
from burnman.classes.material import Material, material_property
from burnman.classes.layer import _radial_spline
from burnman import constants


class UniformMaterial(Material):
    def __init__(self, density):
        super().__init__()
        self.value = density

    @material_property
    def density(self):
        return self.value


class PressureStepMaterial(Material):
    @material_property
    def density(self):
        return 8000.0 if self.pressure >= 1e9 else 4000.0


def uniform_layer(name, inner, outer, density, n=65):
    layer = burnman.Layer(name, np.linspace(inner, outer, n))
    layer.set_material(UniformMaterial(density))
    layer.set_temperature_mode("user-defined", temperatures=np.full(n, 300.0))
    return layer


class test_radial_structure(BurnManTest):
    def test_radial_spline_bounds_density_jumps(self):
        for x in (
            np.linspace(0.0, 1.0, 17),
            np.array([0.0, 0.05, 0.2, 0.3, 0.45, 0.47, 0.6, 0.85, 1.0]),
        ):
            for interface in (0.0, 0.45, 0.99):
                for sign in (-1.0, 1.0):
                    for radius_scale, density_scale in ((1.0, 1e-6), (1e6, 1.0)):
                        with self.subTest(
                            x=x,
                            interface=interface,
                            sign=sign,
                            radius_scale=radius_scale,
                            density_scale=density_scale,
                        ):
                            r = x * radius_scale
                            rho = density_scale * (
                                6000.0
                                + sign * np.where(x <= interface, 2000.0, -2000.0)
                            )
                            spline = _radial_spline(r, rho)
                            self.assertArraysAlmostEqual(spline(r), rho, tol=1e-14)
                            points = r[:-1, None] + np.diff(r)[:, None] * np.linspace(
                                0.0, 1.0, 257
                            )
                            interpolated = spline(points)
                            tolerance = 128.0 * np.finfo(float).eps * max(rho)
                            self.assertTrue(
                                np.all(
                                    interpolated
                                    >= np.minimum(rho[:-1], rho[1:])[:, None]
                                    - tolerance
                                )
                            )
                            self.assertTrue(
                                np.all(
                                    interpolated
                                    <= np.maximum(rho[:-1], rho[1:])[:, None]
                                    + tolerance
                                )
                            )

    def test_radial_spline_preserves_smooth_profiles_and_extrema(self):
        for x in (np.linspace(0.0, 1.0, 9), np.array([0.0, 0.1, 0.25, 0.4, 0.8, 1.0])):
            for profile in (
                lambda r: 4000.0 + 0.0 * r,
                lambda r: 8000.0 - 4000.0 * r**2,
                lambda r: 8000.0 - 4000.0 * r**3,
                lambda r: 8000.0 - 4000.0 * (r - 0.53) ** 2,
                lambda r: 4000.0 + 4000.0 * (r - 0.53) ** 2,
            ):
                with self.subTest(x=x, profile=profile):
                    spline = _radial_spline(x * 1e6, profile(x))
                    points = np.linspace(0.0, 1.0, 1001)
                    self.assertArraysAlmostEqual(
                        spline(points * 1e6), profile(points), tol=1e-14
                    )

    def test_density_jump_uses_consistent_radial_integrals(self):
        from scipy.integrate import quad

        layer = burnman.Layer("density jump", np.linspace(1e6, 2e6, 65))
        layer.set_material(PressureStepMaterial())
        layer.set_temperature_mode("user-defined", temperatures=np.full(65, 300.0))
        with self.assertWarns(UserWarning):
            layer.set_pressure_mode(
                "user-defined", pressures=np.linspace(2e9, 0.0, 65), gravity_bottom=5.0
            )
        layer.make()
        rho = _radial_spline(layer.radii, layer.density)
        gravity = layer.gravity
        g = _radial_spline(layer.radii, gravity)

        def integral(function):
            return quad(
                function,
                layer.inner_radius,
                layer.outer_radius,
                points=layer.radii,
                epsabs=0.0,
                epsrel=1e-11,
                limit=200,
            )[0]

        self.assertFloatEqual(
            layer.mass, 4 * np.pi * integral(lambda r: rho(r) * r**2), tol=1e-13
        )
        self.assertFloatEqual(
            layer.moment_of_inertia,
            8 * np.pi / 3 * integral(lambda r: rho(r) * r**4),
            tol=1e-13,
        )
        self.assertFloatEqual(
            gravity[-1] * layer.outer_radius**2,
            5.0 * layer.inner_radius**2 + constants.G * layer.mass,
            tol=1e-14,
        )
        pressure = layer._compute_pressure(layer.density, gravity, 1e5)
        self.assertFloatEqual(
            pressure[0], 1e5 + integral(lambda r: rho(r) * g(r)), tol=1e-13
        )
        self.assertTrue(np.all(np.diff(pressure) < 0.0))

    def test_reference_core_eos_volume_precision(self):
        from burnman.eos.birch_murnaghan import (
            pressure_fourth_order,
            volume_fourth_order,
        )

        original = burnman.minerals.other.Liquid_Fe_Anderson().params
        for scale in (1e-3, 1.0, 1e3):
            params = dict(original, V_0=original["V_0"] * scale)
            for pressure in (0.0, 1e9, 100e9, 500e9):
                volume = volume_fourth_order(pressure, params)
                result = pressure_fourth_order(params["V_0"] / volume, params)
                self.assertFloatEqual(result, pressure, tol=2e-12, tol_zero=1e-3)

    def test_uniform_sphere(self):
        R, rho, surface = 1e6, 4000.0, 1e5
        for radii in ([0.0, R], [0.0, 0.2 * R, R], np.linspace(0, R, 9)):
            with self.subTest(radii=radii):
                layer = burnman.Layer("sphere", radii)
                layer.set_material(UniformMaterial(rho))
                layer.set_temperature_mode(
                    "user-defined", temperatures=np.full(len(radii), 300.0)
                )
                layer.set_pressure_mode(gravity_bottom=0.0, pressure_top=surface)
                layer.make()
                r = layer.radii
                self.assertArraysAlmostEqual(
                    layer.gravity,
                    4 * np.pi * constants.G * rho * r / 3,
                    tol_zero=1e-13,
                    tol=1e-7,
                )
                self.assertArraysAlmostEqual(
                    layer.pressures,
                    surface + 2 * np.pi * constants.G * rho**2 * (R**2 - r**2) / 3,
                    tol=1e-13,
                )
                mass = 4 * np.pi * rho * R**3 / 3
                self.assertFloatEqual(layer.mass, mass, tol=1e-14)
                self.assertFloatEqual(
                    layer.moment_of_inertia, 0.4 * mass * R**2, tol=1e-14
                )
                self.assertTrue(layer.converged)

    def test_hollow_shell_keeps_inner_gravity(self):
        inner, outer, rho, g0, p0 = 1e6, 2e6, 4000.0, 5.0, 1e5
        layer = uniform_layer("shell", inner, outer, rho, n=257)
        planet = burnman.Planet("shell", [layer], max_delta=1e-12)
        planet.set_pressure_mode(gravity_bottom=g0, pressure_top=p0, max_delta=1e-12)
        planet.make()
        a = 4 * np.pi * constants.G * rho / 3
        b = (g0 - a * inner) * inner**2
        r = layer.radii
        self.assertArraysAlmostEqual(layer.gravity, a * r + b / r**2, tol=1e-14)
        expected = p0 + rho * (a * (outer**2 - r**2) / 2 + b * (1 / r - 1 / outer))
        self.assertArraysAlmostEqual(layer.pressures, expected, tol=1e-8)
        self.assertEqual(layer.gravity_bottom, g0)
        self.assertEqual(planet.gravity[0], g0)

    def test_layer_interfaces_and_unordered_queries(self):
        core = uniform_layer("core", 0.0, 1e6, 8000.0)
        mantle = uniform_layer("mantle", 1e6, 2e6, 4000.0, n=257)
        planet = burnman.Planet("test", [mantle, core], max_delta=1e-12)
        planet.make()
        self.assertEqual(core.gravity_bottom, 0.0)
        self.assertEqual(mantle.gravity_bottom, core.gravity[-1])
        self.assertEqual(core.pressures[-1], mantle.pressures[0])
        self.assertArraysAlmostEqual(
            planet.evaluate(["density"], [1.8e6, 0.2e6, 1e6, 1.8e6]),
            [4000.0, 8000.0, 8000.0, 4000.0],
            tol=1e-7,
        )
        self.assertArraysAlmostEqual(planet.evaluate(["density"], []), [], tol=1e-7)
        for radius in (-1.0, 2e6 + 1.0, np.nan):
            with self.assertRaises(LookupError):
                planet.get_layer_by_radius(radius)
        mass = 4 * np.pi / 3 * (8000.0 * 1e18 + 4000.0 * (8e18 - 1e18))
        self.assertFloatEqual(planet.mass, mass, tol=1e-14)
        self.assertFloatEqual(
            planet.gravity[-1], constants.G * mass / (2e6) ** 2, tol=1e-14
        )

    def test_prescribed_planet_pressures_are_used_on_each_make(self):
        layer = uniform_layer("sphere", 0.0, 1e6, 4000.0)
        planet = burnman.Planet("test", [layer])
        for maximum in (1e9, 2e9):
            pressures = np.linspace(maximum, 0.0, len(planet.radii))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                planet.set_pressure_mode("user-defined", pressures=pressures)
            planet.make()
            self.assertArraysAlmostEqual(planet.pressures, pressures, tol=1e-7)
            self.assertArraysAlmostEqual(layer.pressures, pressures, tol=1e-7)
            self.assertArraysAlmostEqual(
                planet.gravity,
                4 * np.pi * constants.G * 4000.0 * planet.radii / 3,
                tol=1e-7,
            )
            self.assertTrue(planet.converged)

    def test_nonconvergence_is_reported(self):
        layer = burnman.Layer("MgO", np.linspace(660e3, 2890e3, 17))
        layer.set_material(burnman.minerals.SLB_2011.periclase())
        layer.set_temperature_mode("user-defined", temperatures=np.full(17, 2000.0))
        layer.set_pressure_mode(
            gravity_bottom=10.0,
            pressure_top=22.5e9,
            n_max_iterations=1,
            max_delta=1e-12,
        )
        with self.assertWarnsRegex(RuntimeWarning, "did not converge"):
            layer.make()
        self.assertFalse(layer.converged)
        self.assertEqual(layer.n_iterations, 1)

    def test_repeated_planet_make_refreshes_material_properties(self):
        layer = burnman.Layer("MgO", np.linspace(0.0, 1e6, 9))
        layer.set_material(burnman.minerals.SLB_2011.periclase())
        layer.set_temperature_mode("user-defined", temperatures=np.full(9, 1000.0))
        planet = burnman.Planet("test", [layer])
        values = []
        for maximum in (1e9, 2e9):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                planet.set_pressure_mode(
                    "user-defined", pressures=np.linspace(maximum, 0.0, 9)
                )
            planet.make()
            values.append(planet.density[0])
            check = burnman.minerals.SLB_2011.periclase()
            check.set_state(maximum, 1000.0)
            self.assertFloatEqual(planet.density[0], check.density, tol=1e-14)
        self.assertNotEqual(*values)

    def test_layer_profile_against_independent_hydrostatic_solution(self):
        from scipy.integrate import solve_bvp

        n = 513
        layer = burnman.Layer("MgO", np.linspace(660e3, 2890e3, n))
        layer.set_material(burnman.minerals.SLB_2011.periclase())
        layer.set_temperature_mode("user-defined", temperatures=np.full(n, 2000.0))
        layer.set_pressure_mode(
            gravity_bottom=10.0,
            pressure_top=22.5e9,
            max_delta=1e-12,
            n_max_iterations=100,
        )
        layer.make()
        material = layer.material.copy()
        x = layer.radii / 1e6

        # Independently integrate pressure (GPa) and q=(r/1e6)^2*g.
        def derivatives(radius, state):
            density = material.evaluate(
                ["density"], state[0] * 1e9, np.full(len(radius), 2000.0)
            )[0]
            return np.array(
                [
                    -density * state[1] / (1000.0 * radius**2),
                    4.0 * np.pi * constants.G * 1e6 * density * radius**2,
                ]
            )

        def boundaries(inner, outer):
            return np.array([inner[1] - x[0] ** 2 * 10.0, outer[0] - 22.5])

        initial = np.array([layer.pressures / 1e9, layer.gravity * x**2])
        reference = solve_bvp(
            derivatives, boundaries, x, initial, tol=1e-9, max_nodes=10000
        )
        self.assertTrue(reference.success, reference.message)
        pressure, q = reference.sol(x)
        self.assertArraysAlmostEqual(
            layer.pressures, pressure * 1e9, tol=0.0, tol_zero=10.0
        )
        self.assertArraysAlmostEqual(layer.gravity, q / x**2, tol=0.0, tol_zero=1e-8)

    def test_planet_profile_against_independent_hydrostatic_solution(self):
        from scipy.integrate import solve_bvp

        n = 257
        core = burnman.Layer("core", np.linspace(0.0, 3480e3, n))
        core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
        core.set_temperature_mode("user-defined", temperatures=np.full(n, 300.0))
        mantle = burnman.Layer("mantle", np.linspace(3480e3, 6371e3, n))
        mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        mantle.set_temperature_mode("adiabatic", temperature_top=1200.0)
        planet = burnman.Planet(
            "earth_like", [core, mantle], max_delta=1e-12, n_max_iterations=100
        )
        planet.make()
        core_material, mantle_material = core.material.copy(), mantle.material.copy()
        core_radius, mantle_width = 3.480, 2.891  # Million metres.
        x = np.linspace(0.0, 1.0, n)

        # Separate coordinates retain the density jump at the interface.
        # State: core P (GPa), core g, mantle P (GPa), mantle g, mantle T.
        def derivatives(position, state):
            core_density = core_material.evaluate(
                ["density"], state[0] * 1e9, np.full(len(position), 300.0)
            )[0]
            density, gamma, bulk = mantle_material.evaluate(
                ["density", "grueneisen_parameter", "isentropic_bulk_modulus_reuss"],
                state[2] * 1e9,
                state[4],
            )
            pressure_gradient = -mantle_width / 1000.0 * density * state[3]
            return np.array(
                [
                    -core_radius / 1000.0 * core_density * state[1],
                    4.0 * np.pi * constants.G * core_radius * 1e6 * core_density,
                    pressure_gradient,
                    4.0 * np.pi * constants.G * mantle_width * 1e6 * density
                    - 2.0
                    * mantle_width
                    * state[3]
                    / (core_radius + mantle_width * position),
                    gamma * state[4] / bulk * pressure_gradient * 1e9,
                ]
            )

        def boundaries(inner, outer):
            return np.array(
                [
                    inner[1],
                    outer[2],
                    outer[4] - 1200.0,
                    outer[0] - inner[2],
                    outer[1] - inner[3],
                ]
            )

        initial = np.array(
            [
                core.pressures / 1e9,
                core.gravity,
                mantle.pressures / 1e9,
                mantle.gravity,
                mantle.temperatures,
            ]
        )
        singular = np.zeros((5, 5))
        singular[1, 1] = -2.0  # Regular central gravity: dg/dx = 4*pi*G*rho*Rc - 2g/x.
        reference = solve_bvp(
            derivatives, boundaries, x, initial, S=singular, tol=1e-8, max_nodes=10000
        )
        self.assertTrue(reference.success, reference.message)
        pressure_core, gravity_core, pressure_mantle, gravity_mantle, temperature = (
            reference.sol(x)
        )
        self.assertArraysAlmostEqual(
            core.pressures, pressure_core * 1e9, tol=0.0, tol_zero=10.0
        )
        self.assertArraysAlmostEqual(
            mantle.pressures, pressure_mantle * 1e9, tol=0.0, tol_zero=10.0
        )
        self.assertArraysAlmostEqual(core.gravity, gravity_core, tol=0.0, tol_zero=1e-8)
        self.assertArraysAlmostEqual(
            mantle.gravity, gravity_mantle, tol=0.0, tol_zero=1e-8
        )
        self.assertArraysAlmostEqual(
            mantle.temperatures, temperature, tol=0.0, tol_zero=1e-6
        )


if __name__ == "__main__":
    unittest.main()
