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


class test_layer(BurnManTest):
    def layer1(self):
        rock = min1()
        rock.set_method("slb3")

        depths = np.linspace(660.0e3, 2890.0e3, 10)
        temperatures = np.linspace(2000.0, 2000.0, 10)

        layer = Layer("m1", depths)
        layer.set_material(rock)
        layer.set_temperature_mode("user-defined", temperatures=temperatures)
        layer.set_pressure_mode(
            "self-consistent", pressure_top=22.5e9, gravity_bottom=10.0
        )
        layer.make()
        return layer

    def test_pressures(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(m.pressures[0::9] / 1.0e9, [61.62679132, 22.5])

    def test_vs1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.shear_wave_velocity[0::9] / 1.0e3, [6.72149374, 5.94945973]
        )

    def test_vp1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.p_wave_velocity[0::9] / 1.0e3, [11.82332527, 10.21538394]
        )

    def test_vphi1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(
            m.bulk_sound_velocity[0::9] / 1.0e3, [8.91925163, 7.56037749]
        )

    def test_evaluate(self):
        m = self.layer1()

        d = m.evaluate(
            ["thermal_expansivity", "molar_heat_capacity_p"],
            [700.0e3, 800.0e3, 750.0e3],
        )

        d0 = [1.90903185e-05, 1.97610550e-05, 1.94185684e-05]
        d1 = [5.14621379e01, 5.15786478e01, 5.15191396e01]

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
