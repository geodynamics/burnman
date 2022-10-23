from __future__ import absolute_import
import unittest
from util import BurnManTest

import burnman.eos.property_modifiers as pm
from burnman.minerals import SLB_2011, HP_2011_ds62
from burnman.tools.eos import check_eos_consistency

linear_params = {"delta_E": 1200.0, "delta_S": 5.0, "delta_V": 1.0e-7}
landau_params = {"Tc_0": 800.0, "S_D": 5.0, "V_D": 1.0e-7}
landau_params_2 = {"Tc_0": 1200.0, "S_D": 5.0, "V_D": 1.0e-7}
landau_hp_params = {
    "P_0": 1.0e5,
    "T_0": 298.15,
    "Tc_0": 800.0,
    "S_D": 5.0,
    "V_D": 1.0e-7,
}
landau_hp_params_2 = {
    "P_0": 1.0e5,
    "T_0": 298.15,
    "Tc_0": 1200.0,
    "S_D": 5.0,
    "V_D": 1.0e-7,
}
bragg_williams_params = {
    "n": 1.0,
    "factor": 0.8,
    "Wh": 1000.0,
    "Wv": 1.0e-7,
    "deltaH": 1000.0,
    "deltaV": 1.0e-7,
}
magnetic_params = {
    "structural_parameter": 0.4,
    "curie_temperature": [800.0, 1.0e-8],
    "magnetic_moment": [2.2, 1.0e-10],
}
magnetic_params_2 = {
    "structural_parameter": 0.4,
    "curie_temperature": [1200.0, 1.0e-8],
    "magnetic_moment": [2.2, 1.0e-10],
}
thermal_params = {"Theta_0": 1200.0, "Cv_inf": 1.0}
thermal_delta_params = {"Theta_0": 1200.0, "S_inf": 1.0}


fn_params = [
    [pm._linear_excesses, "linear", linear_params],
    [pm._landau_excesses, "landau", landau_params],
    [pm._landau_excesses, "landau", landau_params_2],
    [pm._landau_hp_excesses, "landau_hp", landau_hp_params],
    [pm._landau_hp_excesses, "landau_hp", landau_hp_params_2],
    [pm._bragg_williams_excesses, "bragg_williams", bragg_williams_params],
    [pm._magnetic_excesses_chs, "magnetic_chs", magnetic_params],
    [pm._magnetic_excesses_chs, "magnetic_chs", magnetic_params_2],
    [pm._debye_excesses, "debye", thermal_params],
    [pm._debye_delta_excesses, "debye_delta", thermal_delta_params],
    [pm._einstein_excesses, "einstein", thermal_params],
    [pm._einstein_delta_excesses, "einstein_delta", thermal_delta_params],
]


class Modifiers(BurnManTest):
    def test_excess_functions(self):
        P = 1.0e11
        T = 1000.0

        gibbs_excesses = [
            6200.0,
            1137.8563,
            1019.3716,
            -2534.185,
            -1696.3556,
            -551.3463,
            -14069.1506,
            -21582.5564,
            -565.3150,
            -620.7977,
            -358.3824,
            -517.2153,
        ]
        gibbs_excesses_output = []
        for (fn, name, params) in fn_params:
            excesses = fn(P, T, params)
            gibbs_excesses_output.append(round(excesses[0]["G"], 4))

        self.assertArraysAlmostEqual(gibbs_excesses_output, gibbs_excesses)

    def test_modifier_with_pvt_eos(self):
        qtz = SLB_2011.quartz()
        assert check_eos_consistency(qtz, 1.0e5, 600.0)
        assert check_eos_consistency(qtz, 1.0e5, 1100.0)
        assert check_eos_consistency(qtz, 1.0e9, 1100.0)

    def test_modifier_with_vpt_eos(self):
        qtz = HP_2011_ds62.q()
        assert check_eos_consistency(
            qtz, 1.0e5, 600.0, including_shear_properties=False
        )
        assert check_eos_consistency(
            qtz, 1.0e5, 1100.0, including_shear_properties=False
        )
        assert check_eos_consistency(
            qtz, 1.0e9, 1100.0, including_shear_properties=False
        )

    def test_consistency(self):
        per = SLB_2011.periclase()
        for (fn, name, params) in fn_params:
            per.property_modifiers = [[name, params]]
            # Bragg-Williams currently calculates derivatives
            # numerically, so the modifier derivatives aren't
            # sufficiently accurate to test here.
            # TODO: Fix this in property_modifiers.py
            if name != "bragg_williams":
                assert check_eos_consistency(
                    per, 1.0e9, 200.0, including_shear_properties=False
                )
                assert check_eos_consistency(
                    per, 1.0e9, 1100.0, including_shear_properties=False
                )


if __name__ == "__main__":
    unittest.main()
