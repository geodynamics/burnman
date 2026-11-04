import unittest

from util import BurnManTest
import numpy as np

import burnman
from burnman import Planet
from burnman import Layer


def make_simple_planet():
    core = Layer("core", np.linspace(0.0e3, 3480.0e3, 10))
    core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
    core.set_temperature_mode(
        "user-defined", temperatures=300.0 * np.ones_like(core.radii)
    )

    mantle = Layer("mantle", np.linspace(3480.0e3, 6371.0e3, 10))
    mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
    mantle.set_temperature_mode("adiabatic", temperature_top=1200)
    myplanet = Planet("earth_like", [core, mantle], max_delta=1e-11)
    myplanet.make()
    return myplanet, core, mantle


class test_planet(BurnManTest):
    def test_extensive_profiles_follow_each_layer_material_amount(self):
        core = Layer("core", np.linspace(0.0, 5e5, 4))
        mantle = Layer("mantle", np.linspace(5e5, 1e6, 4))
        for layer, material, amount in (
            (core, burnman.minerals.SLB_2011.forsterite(), 0.4),
            (mantle, burnman.minerals.SLB_2011.periclase(), 2.5),
        ):
            material.number_of_moles = amount
            layer.set_material(material)
            layer.set_temperature_mode("user-defined", temperatures=np.full(4, 1500.0))
        planet = Planet("scaled samples", [core, mantle])
        planet.make()
        physical_mass, physical_volume = planet.mass, planet.volume
        self.assertEqual(np.ndim(physical_mass), 0)
        self.assertEqual(np.ndim(physical_volume), 0)
        self.assertFloatEqual(physical_mass, sum(layer.mass for layer in planet.layers))
        self.assertFloatEqual(
            physical_volume, sum(layer.volume for layer in planet.layers)
        )
        amounts = np.r_[np.full(4, 0.4), np.full(4, 2.5)]
        totals = {
            "energy": "molar_internal_energy",
            "internal_energy": "molar_internal_energy",
            "gibbs": "molar_gibbs",
            "helmholtz": "molar_helmholtz",
            "V": "molar_volume",
            "material_volume": "molar_volume",
            "material_mass": "molar_mass",
            "S": "molar_entropy",
            "H": "molar_enthalpy",
            "C_p": "molar_heat_capacity_p",
            "C_v": "molar_heat_capacity_v",
        }
        for name, molar in totals.items():
            with self.subTest(property=name):
                expected = amounts * getattr(planet, molar)
                self.assertArraysAlmostEqual(getattr(planet, name), expected)
                self.assertArraysAlmostEqual(planet.evaluate([name]), expected)
                radii = np.array([2e5, 8e5])
                self.assertArraysAlmostEqual(
                    planet.evaluate([name], radii),
                    np.array([0.4, 2.5]) * planet.evaluate([molar], radii),
                )
        for name, value in (("mass", physical_mass), ("volume", physical_volume)):
            self.assertArraysAlmostEqual(
                planet.evaluate([name]), np.full(len(planet.radii), value)
            )
            self.assertArraysAlmostEqual(
                planet.evaluate([name], radii), np.full(len(radii), value)
            )
        self.assertFloatEqual(
            planet.volume, 4.0 * np.pi * planet.radius_planet**3 / 3.0
        )
        for layer in planet.layers:
            layer.material.number_of_moles = 1.0
        planet.make()
        self.assertFloatEqual(planet.mass, physical_mass)
        self.assertEqual(planet.volume, physical_volume)

    def test_planet_1(self):
        myplanet, core, mantle = make_simple_planet()
        assert myplanet.get_layer("core") == core
        assert myplanet.get_layer_by_radius(2000.0e3) == core
        assert myplanet.get_layer_by_radius(4000.0e3) == mantle
        assert myplanet.get_layer_by_radius(5000.0e3) == mantle
        self.assertFloatEqual(
            myplanet.pressure[0],
            445231466723.2055,
            tol=1.0e-4,
        )

    def test_sort_layers(self):
        myplanet, core, mantle = make_simple_planet()
        assert myplanet.layers == [core, mantle]

    def test_temperature(self):
        myplanet, core, mantle = make_simple_planet()
        self.assertArraysAlmostEqual(mantle.S, [mantle.S[-1] for S in mantle.S])

    def test_evaluate(self):
        myplanet, core, mantle = make_simple_planet()

        # point-by-point properties
        property_names = [
            "depth",
            "gravity",
            "temperature",
            "molar_internal_energy",
            "molar_gibbs",
            "molar_helmholtz",
            "molar_mass",
            "molar_volume",
            "molar_entropy",
            "molar_enthalpy",
            "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss",
            "isothermal_compressibility_reuss",
            "isentropic_compressibility_reuss",
            "bulk_sound_velocity",
            "grueneisen_parameter",
            "thermal_expansivity",
            "molar_heat_capacity_v",
            "molar_heat_capacity_p",
            "P",
            "T",
            "energy",
            "helmholtz",
            "gibbs",
            "V",
            "rho",
            "S",
            "H",
            "K_T",
            "K_S",
            "beta_T",
            "beta_S",
            "v_phi",
            "gr",
            "alpha",
            "C_v",
            "C_p",
        ]
        properties = myplanet.evaluate(property_names, myplanet.radii)

        for i, name in enumerate(property_names):
            self.assertFloatEqual(properties[i][-2], getattr(myplanet, name)[-2])

        # layer-by-layer properties
        property_names = [
            "average_density",
            "mass",
            "moment_of_inertia",
            "moment_of_inertia_factor",
        ]
        for name in property_names:
            self.assertTrue(isinstance(getattr(myplanet, name), float))

        # solid properties
        property_names = [
            "bullen",
            "brunt_vasala",
            "shear_modulus",
            "p_wave_velocity",
            "shear_wave_velocity",
            "G",
            "v_p",
            "v_s",
        ]
        properties = myplanet.evaluate(property_names, myplanet.layers[1].radii[1:-1])

        for i, name in enumerate(property_names):
            self.assertFloatEqual(
                properties[i][-1], getattr(myplanet.layers[1], name)[-2]
            )

    def test_evaluate_solid_planet(self):
        myplanet, core, mantle = make_simple_planet()
        core.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        myplanet.make()

        # point-by-point properties
        property_names = [
            "depth",
            "gravity",
            "temperature",
            "molar_internal_energy",
            "molar_gibbs",
            "molar_helmholtz",
            "molar_mass",
            "molar_volume",
            "molar_entropy",
            "molar_enthalpy",
            "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss",
            "isothermal_compressibility_reuss",
            "isentropic_compressibility_reuss",
            "bulk_sound_velocity",
            "grueneisen_parameter",
            "thermal_expansivity",
            "molar_heat_capacity_v",
            "molar_heat_capacity_p",
            "P",
            "T",
            "energy",
            "helmholtz",
            "gibbs",
            "V",
            "rho",
            "S",
            "H",
            "K_T",
            "K_S",
            "beta_T",
            "beta_S",
            "v_phi",
            "gr",
            "alpha",
            "C_v",
            "C_p",
            "bullen",
            "brunt_vasala",
            "shear_modulus",
            "p_wave_velocity",
            "shear_wave_velocity",
            "G",
            "v_p",
            "v_s",
        ]
        properties = myplanet.evaluate(property_names, myplanet.radii)

        for i, name in enumerate(property_names):
            self.assertFloatEqual(properties[i][-2], getattr(myplanet, name)[-2])

        # layer-by-layer properties
        property_names = [
            "average_density",
            "mass",
            "moment_of_inertia",
            "moment_of_inertia_factor",
        ]
        for name in property_names:
            self.assertTrue(isinstance(getattr(myplanet, name), float))


if __name__ == "__main__":
    unittest.main()
