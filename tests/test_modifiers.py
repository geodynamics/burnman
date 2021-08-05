from __future__ import absolute_import
import unittest
from util import BurnManTest

import burnman_path
import burnman.eos.property_modifiers as pm

assert burnman_path  # silence pyflakes warning


class Modifiers(BurnManTest):

    def test_excess_functions(self):
        linear_params = {'delta_E': 1200., 'delta_S': 5., 'delta_V': 1.e-7}
        landau_params = {'Tc_0': 800., 'S_D': 5., 'V_D': 1.e-7}
        landau_params_2 = {'Tc_0': 1200., 'S_D': 5., 'V_D': 1.e-7}
        landau_hp_params = {
            'P_0': 1.e5, 'T_0': 298.15, 'Tc_0': 800., 'S_D': 5., 'V_D': 1.e-7}
        landau_hp_params_2 = {
            'P_0': 1.e5, 'T_0': 298.15, 'Tc_0': 1200., 'S_D': 5., 'V_D': 1.e-7}
        bragg_williams_params = {
            'n': 1., 'factor': 0.8, 'Wh': 1000., 'Wv': 1.e-7, 'deltaH': 1000., 'deltaV': 1.e-7}
        magnetic_params = {'structural_parameter': 0.4, 'curie_temperature': [
            800., 1.e-8], 'magnetic_moment': [2.2, 1.e-10]}
        magnetic_params_2 = {'structural_parameter': 0.4, 'curie_temperature': [
            1200., 1.e-8], 'magnetic_moment': [2.2, 1.e-10]}

        P = 1.e11
        T = 1000.

        linear_excesses = pm._linear_excesses(P, T, linear_params)
        landau_excesses = pm._landau_excesses(P, T, landau_params)
        landau_excesses_2 = pm._landau_excesses(P, T, landau_params_2)
        landau_hp_excesses = pm._landau_hp_excesses(P, T, landau_hp_params)
        landau_hp_excesses_2 = pm._landau_hp_excesses(P, T, landau_hp_params_2)
        bragg_williams_excesses = pm._bragg_williams_excesses(
            P, T, bragg_williams_params)
        magnetic_excesses = pm._magnetic_excesses_chs(P, T, magnetic_params)
        magnetic_excesses_2 = pm._magnetic_excesses_chs(
            P, T, magnetic_params_2)

        self.assertFloatEqual(linear_excesses[0]['G'],
                              1200. - 5.*T + 1.e-7*P)

        gibbs_excesses = [6200.0, 1137.8563, 1019.3716, -2534.185,
                          -1696.3556, -551.3463, -14069.1506, -21582.5564]
        gibbs_excesses_output = []
        for excesses in [linear_excesses,
                         landau_excesses, landau_excesses_2,
                         landau_hp_excesses, landau_hp_excesses_2,
                         bragg_williams_excesses,
                         magnetic_excesses, magnetic_excesses_2]:

            gibbs_excesses_output.append(round(excesses[0]['G'], 4))

        self.assertArraysAlmostEqual(gibbs_excesses_output, gibbs_excesses)


if __name__ == '__main__':
    unittest.main()
