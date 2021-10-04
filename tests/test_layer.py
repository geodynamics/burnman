from __future__ import absolute_import
import unittest
from util import BurnManTest

import numpy as np
import burnman_path
import burnman
from burnman import Layer

assert burnman_path  # silence pyflakes warning


class min1 (burnman.Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8}
        burnman.Mineral.__init__(self)


class min2 (min1):

    def __init__(self):
        min1.__init__(self)
        self.params['V_0'] = 10e-6


class test_layer(BurnManTest):

    def layer1(self):
        rock = min1()
        rock.set_method('slb3')

        depths = np.linspace(660.e3, 2890.e3, 10)
        temperatures = np.linspace(2000., 2000., 10)

        layer = Layer("m1", depths)
        layer.set_material(rock)
        layer.set_temperature_mode('user-defined', temperatures=temperatures)
        layer.set_pressure_mode('self-consistent', pressure_top=22.5e9,
                                gravity_bottom=10.)
        layer.make()
        return layer

    def test_pressures(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(m.pressures[0::9]/1.e9, [61.62679132, 22.5])

    def test_vs1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(m.shear_wave_velocity[0::9]/1.e3,
                                     [6.72149374, 5.94945973])

    def test_vp1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(m.p_wave_velocity[0::9]/1.e3,
                                     [11.82332527, 10.21538394])

    def test_vphi1(self):
        m = self.layer1()
        self.assertArraysAlmostEqual(m.bulk_sound_velocity[0::9]/1.e3,
                                     [8.91925163, 7.56037749])

    def test_evaluate(self):
        m = self.layer1()

        d = m.evaluate(['thermal_expansivity', 'molar_heat_capacity_p'],
                       [700.e3, 800.e3, 750.e3])

        d0 = [1.90903185e-05, 1.97610550e-05, 1.94185684e-05]
        d1 = [5.14621379e+01, 5.15786478e+01, 5.15191396e+01]

        self.assertArraysAlmostEqual(d[0], d0)
        self.assertArraysAlmostEqual(d[1], d1)


if __name__ == '__main__':
    unittest.main()
