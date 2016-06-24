from __future__ import absolute_import

import os
import sys
sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman.eos.property_modifiers as pm

import unittest
from util import BurnManTest


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
            'n': 1., 'factor': 0.8, 'Wh': 13000., 'Wv': 1.e-7, 'deltaH': 13000., 'deltaV': 1.e-7}
        magnetic_params = {'structural_parameter': 0.4, 'curie_temperature': [
            800., 1.e-8], 'magnetic_moment': [2.2, 1.e-10]}
        magnetic_params_2 = {'structural_parameter': 0.4, 'curie_temperature': [
            1200., 1.e-8], 'magnetic_moment': [2.2, 1.e-10]}

        linear_excesses = []
        landau_excesses = []
        landau_excesses_2 = []
        landau_hp_excesses = []
        landau_hp_excesses_2 = []
        bragg_williams_excesses = []
        magnetic_excesses = []
        magnetic_excesses_2 = []

        P = 1.e11
        T = 1000.

        dP = 1000.
        dT = 1.

        pressures = [P - dP, P, P + dP]
        temperatures = [T - dT, T, T + dT]
        for P in pressures:
            for T in temperatures:
                linear_excesses.append(
                    pm._linear_excesses(P, T, linear_params))
                landau_excesses.append(
                    pm._landau_excesses(P, T, landau_params))
                landau_excesses_2.append(
                    pm._landau_excesses(P, T, landau_params_2))
                landau_hp_excesses.append(
                    pm._landau_hp_excesses(P, T, landau_hp_params))
                landau_hp_excesses_2.append(
                    pm._landau_hp_excesses(P, T, landau_hp_params_2))
                bragg_williams_excesses.append(
                    pm._bragg_williams_excesses(P, T, bragg_williams_params))
                magnetic_excesses.append(
                    pm._magnetic_excesses_chs(P, T, magnetic_params))
                magnetic_excesses_2.append(
                    pm._magnetic_excesses_chs(P, T, magnetic_params_2))

        all_excesses = [['linear', linear_excesses],
                        ['landau below Tc', landau_excesses],
                        ['landau above Tc', landau_excesses_2],
                        ['landau_hp below Tc', landau_hp_excesses],
                        ['landau_hp above Tc', landau_hp_excesses_2],
                        ['bragg-williams', bragg_williams_excesses],
                        ['magnetic below Tc', magnetic_excesses],
                        ['magnetic above Tc', magnetic_excesses_2]]

        """
        0) P-dP, T-dT
        1) P-dP, T
        2) P-dP, T+dT
        3) P,    T-dT
        4) P,    T
        5) P,    T+dT
        6) P+dP, T-dT
        7) P+dP, T
        8) P+dP, T+dT
        """
        numerical = []
        analytical = []
        for excess_type, excesses in all_excesses:
            numerical.append((excesses[5]['G'] - excesses[3]['G']) / (2. * dT))
            analytical.append(excesses[4]['dGdT'])
            numerical.append((excesses[7]['G'] - excesses[1]['G']) / (2. * dP))
            analytical.append(excesses[4]['dGdP'])
            numerical.append(
                (excesses[5]['dGdT'] - excesses[3]['dGdT']) / (2. * dT))
            analytical.append(excesses[4]['d2GdT2'])
            numerical.append(
                (excesses[7]['dGdP'] - excesses[1]['dGdP']) / (2. * dP))
            analytical.append(excesses[4]['d2GdP2'])
            numerical.append(
                (excesses[5]['dGdP'] - excesses[3]['dGdP']) / (2. * dT))
            analytical.append(excesses[4]['d2GdPdT'])

        self.assertArraysAlmostEqual(numerical, analytical)

if __name__ == '__main__':
    unittest.main()
