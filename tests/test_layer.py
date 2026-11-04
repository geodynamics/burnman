import unittest
from util import BurnManTest

import numpy as np
import burnman
from burnman import Layer, BoundaryLayerPerturbation


class min1(burnman.Mineral):
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


class min2(min1):
    def __init__(self):
        min1.__init__(self)
        self.params["V_0"] = 10e-6


# Reference values from independent hydrostatic boundary-value integration.


class test_layer(BurnManTest):
    def test_extensive_profiles_follow_material_amount(self):
        material = burnman.minerals.SLB_2011.forsterite()
        material.number_of_moles = 2.5
        layer = Layer("sample", np.linspace(1e6, 2e6, 5))
        layer.set_material(material)
        layer.set_temperature_mode(
            "user-defined", temperatures=np.linspace(1000.0, 1500.0, 5)
        )
        layer.set_pressure_mode(
            "user-defined", pressures=np.linspace(6e9, 4e9, 5), gravity_bottom=10.0
        )
        layer.make()
        physical_mass, physical_volume = layer.mass, layer.volume
        self.assertEqual(np.ndim(physical_mass), 0)
        self.assertEqual(np.ndim(physical_volume), 0)
        self.assertFloatEqual(physical_volume, 4.0 * np.pi / 3.0 * (2e6**3 - 1e6**3))
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
                expected = 2.5 * getattr(layer, molar)
                self.assertArraysAlmostEqual(getattr(layer, name), expected)
                self.assertArraysAlmostEqual(layer.evaluate([name]), expected)
                radii = np.array([1.2e6, 1.8e6])
                self.assertArraysAlmostEqual(
                    layer.evaluate([name], radii),
                    material.evaluate(
                        [molar],
                        np.interp(radii, layer.radii, layer.pressures),
                        np.interp(radii, layer.radii, layer.temperatures),
                    )[0]
                    * 2.5,
                )
        for name, value in (("mass", physical_mass), ("volume", physical_volume)):
            self.assertArraysAlmostEqual(
                layer.evaluate([name]), np.full(len(layer.radii), value)
            )
            self.assertArraysAlmostEqual(
                layer.evaluate([name], radii), np.full(len(radii), value)
            )
        material.number_of_moles = 1.0
        layer.make()
        self.assertFloatEqual(layer.mass, physical_mass)
        self.assertEqual(layer.volume, physical_volume)

    def layer1(self):
        rock = min1()
        rock.set_method("slb3")

        radii = np.linspace(660.0e3, 2890.0e3, 257)
        temperatures = np.full(len(radii), 2000.0)

        layer = Layer("m1", radii)
        layer.set_material(rock)
        layer.set_temperature_mode("user-defined", temperatures=temperatures)
        layer.set_pressure_mode(
            "self-consistent", pressure_top=22.5e9, gravity_bottom=10.0, max_delta=1e-11
        )
        layer.make()
        return layer

    def test_pressures(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.pressures[[0, -1]] / 1.0e9, [60.95639688227379, 22.5]
        )

    def test_vs1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.shear_wave_velocity[[0, -1]] / 1.0e3,
            [6.711860911802672, 5.94945976586669],
        )

    def test_vp1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.p_wave_velocity[[0, -1]] / 1.0e3,
            [11.801888100498513, 10.215384021484955],
        )

    def test_vphi1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.bulk_sound_velocity[[0, -1]] / 1.0e3,
            [8.900512730783674, 7.560377550020552],
        )

    def test_evaluate(self):
        m = self.layer1()

        d = m.evaluate(
            ["thermal_expansivity", "molar_heat_capacity_p"],
            [700.0e3, 800.0e3, 750.0e3],
        )

        d0 = [1.9266858705237977e-05, 1.9961283585897523e-05, 1.9626669161472518e-05]
        d1 = [51.49279132548997, 51.613458411319186, 51.555292976809795]

        self.assertArraysAlmostEqual(d[0], d0)
        self.assertArraysAlmostEqual(d[1], d1)

    def test_properties(self):
        m = self.layer1()

        self.assertTrue(isinstance(m.mass, np.floating))
        self.assertTrue(isinstance(m.moment_of_inertia, np.floating))
        self.assertTrue(isinstance(m.gravity, np.ndarray))
        self.assertTrue(isinstance(m.bullen, np.ndarray))
        self.assertTrue(isinstance(m.brunt_vasala, np.ndarray))
        self.assertTrue(isinstance(m.P, np.ndarray))
        self.assertTrue(isinstance(m.T, np.ndarray))
        self.assertTrue(isinstance(m.energy, np.ndarray))
        self.assertTrue(isinstance(m.gibbs, np.ndarray))
        self.assertTrue(isinstance(m.helmholtz, np.ndarray))
        self.assertTrue(isinstance(m.H, np.ndarray))
        self.assertTrue(isinstance(m.V, np.ndarray))
        self.assertTrue(isinstance(m.molar_mass, np.ndarray))
        self.assertTrue(isinstance(m.rho, np.ndarray))
        self.assertTrue(isinstance(m.K_T, np.ndarray))
        self.assertTrue(isinstance(m.K_S, np.ndarray))
        self.assertTrue(isinstance(m.beta_T, np.ndarray))
        self.assertTrue(isinstance(m.beta_S, np.ndarray))
        self.assertTrue(isinstance(m.G, np.ndarray))
        self.assertTrue(isinstance(m.v_phi, np.ndarray))
        self.assertTrue(isinstance(m.gr, np.ndarray))
        self.assertTrue(isinstance(m.alpha, np.ndarray))
        self.assertTrue(isinstance(m.C_p, np.ndarray))
        self.assertTrue(isinstance(m.C_v, np.ndarray))

    def test_tbl(self):
        p = BoundaryLayerPerturbation(
            radius_bottom=3480.0e3,
            radius_top=5711.0e3,
            rayleigh_number=1.0e7,
            temperature_change=840.0,
            boundary_layer_ratio=0.0,
        )

        Ts = p.temperature(np.array([3480.0e3, 5711.0e3]))
        self.assertArraysAlmostEqual(Ts, [840.0, 0.0])

        p = BoundaryLayerPerturbation(
            radius_bottom=3480.0e3,
            radius_top=5711.0e3,
            rayleigh_number=1.0e7,
            temperature_change=1000.0,
            boundary_layer_ratio=0.25,
        )

        Ts = p.temperature(np.array([3480.0e3, 5711.0e3]))
        self.assertArraysAlmostEqual(Ts, [800.0, -200.0])

        p.set_model_thermal_gradients(-18.0 / 1.0e3, -0.6 / 1.0e3)

        dTdrs = p.dTdr(np.array([3480.0e3, 5711.0e3]))
        self.assertArraysAlmostEqual(dTdrs, [-18.0 / 1.0e3, -0.6 / 1.0e3])


if __name__ == "__main__":
    unittest.main()
