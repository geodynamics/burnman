import unittest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman.mineral import Mineral
from util import BurnManTest
from burnman.tools import *
        
class test_tools(BurnManTest):

    def test_eqm_T(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        H_ex = 1.e3
        S_ex = 2.
        fo2.params['H_0'] = fo.params['H_0'] + H_ex
        fo2.params['S_0'] = fo.params['S_0'] + S_ex

        T = H_ex / S_ex
        T_calc = equilibrium_temperature(
            [fo, fo2], [1.0, -1.0], fo.params['P_0'], 1200.)

        self.assertArraysAlmostEqual([T], [T_calc])

    def test_eqm_P(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        fo.params['K_0'] = 1.e20
        fo2.params['K_0'] = 1.e20

        H_ex = 1
        V_ex = 1.e-8
        fo2.params['H_0'] = fo.params['H_0'] + H_ex
        fo2.params['V_0'] = fo.params['V_0'] - V_ex

        P = fo.params['P_0'] + H_ex / V_ex
        P_calc = equilibrium_pressure([fo, fo2], [1.0, -1.0], fo.params['T_0'])
        self.assertArraysAlmostEqual([P], [P_calc])

    def test_hugoniot(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        T_ref = 298.15
        P_ref = 1.e5
        pressures = np.array([P_ref])
        temperatures, volumes = hugoniot(fo, P_ref, T_ref, pressures)

        self.assertArraysAlmostEqual(temperatures, np.array([T_ref]))

    def test_fraction_conversion_mass(self):
        pv = burnman.minerals.HP_2011_ds62.mpv()
        en = burnman.minerals.HP_2011_ds62.en()
        c = burnman.Composite([pv, en], [0.5, 0.5])
        molar_fractions = convert_fractions(c, [0.5, 0.5], 'mass', 'molar')
        self.assertArraysAlmostEqual(molar_fractions, [2. / 3., 1. / 3.])

    def test_fraction_conversion_mass_2(self):
        per = burnman.minerals.SLB_2011.periclase()
        stv = burnman.minerals.SLB_2011.stishovite()
        fo = burnman.minerals.SLB_2011.forsterite()
        c = burnman.Composite([per, stv, fo], [0.5, 0.25, 0.25])
        mass_fractions = convert_fractions(
            c, [0.5, 0.25, 0.25], 'molar', 'mass')
        self.assertArraysAlmostEqual(
            [mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]])

    def test_fraction_conversion_volume(self):
        per = burnman.minerals.SLB_2011.periclase()
        c = burnman.Composite([per, per, per], [0.25, 0.25, 0.5])
        c.set_state(1.e5, 300.)
        mass_fractions = convert_fractions(
            c, [0.25, 0.25, 0.5], 'molar', 'volume')
        self.assertArraysAlmostEqual(
            [mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]])

    def test_bracket(self):
        def fn(x):
            return (x - 1.) * (x - 2.)

        # Test that it can find the root at one from the right
        sol = bracket(fn, 1.2, 1.e-2)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[1] < 1.0)
        self.assertTrue(sol[0] > 1.0)

        # Test that it can find the root at one from the left
        sol = bracket(fn, 0.6, 1.e-2)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[0] < 1.0)
        self.assertTrue(sol[1] > 1.0)

        # Test that it can find the root at two from the left
        sol = bracket(fn, 1.7, 1.e-2)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[0] < 2.0)
        self.assertTrue(sol[1] > 2.0)

        # Test that it can find the root at two from the right
        sol = bracket(fn, 2.5, 1.e-2)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[1] < 2.0)
        self.assertTrue(sol[0] > 2.0)

        # Test that it can find the root at one if we take too large of a dx
        sol = bracket(fn, 1.2, 2.0)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[1] < 1.0)
        self.assertTrue(sol[0] > 1.0)

        # Test that it can find the root at two if we take too large of a dx
        sol = bracket(fn, 1.8, 2.0)
        self.assertTrue(sol[2] * sol[3] < 0.)
        self.assertTrue(sol[0] < 2.0)
        self.assertTrue(sol[1] > 2.0)

    def test_bracket_failure(self):

        mineral = burnman.minerals.SLB_2011.fayalite()
        # This should be too high pressure for the EoS
        mineral.set_state(300.e9, 300.)

        def fn():
            return mineral.molar_volume
        with np.errstate(all='ignore'):
            self.assertRaises(Exception, fn)

    def test_padding_1D(self):
        array = np.array([1., 2., 3., 5.])
        padding = (1,)
        padded_array = burnman.tools._pad_ndarray_inverse_mirror(array, padding)
        self.assertArraysAlmostEqual(padded_array, np.array([0., 1., 2., 3., 5., 7.]))
        
    def test_padding_2D(self):
        array = np.array([[1., 2.], [3., 5.]])
        padding = (1,1)
        padded_array = burnman.tools._pad_ndarray_inverse_mirror(array, padding)
        self.assertArraysAlmostEqual(np.ndarray.flatten(padded_array),
                                     np.ndarray.flatten(np.array([[-3., -1., -1., 1.],
                                                                  [0., 1., 2., 3.],
                                                                  [1., 3., 5., 7.],
                                                                  [4., 5., 8., 9.]])))
        
    def test_smoothing(self):
        array = np.array([0., 1., 2., 3., 4.])
        smoothed_array = smooth_array(array, [1.], [1.])
        self.assertArraysAlmostEqual(array, smoothed_array)

    def test_interp_smoothing_ij(self):
        array = np.array([[0., 1., 2.], [0., 1., 2.], [0., 1., 2.]])
        axis_values = np.array([0., 1., 2.])
        f, dfdx, dfdy = interp_smoothed_array_and_derivatives(array,
                                                              axis_values, axis_values,
                                                              0, 0, indexing='ij')
        self.assertArraysAlmostEqual([f(0., 1.)[0], dfdx(0., 1.)[0], dfdy(0., 1.)[0]],
                                     [array[0][1], 0., 1.])
        
    def test_interp_smoothing_xy(self):
        array = np.array([[0., 1., 2.], [0., 1., 2.], [0., 1., 2.]])
        axis_values = np.array([0., 1., 2.])
        f, dfdx, dfdy = interp_smoothed_array_and_derivatives(array,
                                                              axis_values, axis_values,
                                                              0, 0, indexing='xy')
        self.assertArraysAlmostEqual([f(0., 1.)[0], dfdx(0., 1.)[0], dfdy(0., 1.)[0]],
                                     [array[1][0], 1., 0.])
        
        
if __name__ == '__main__':
    unittest.main()
