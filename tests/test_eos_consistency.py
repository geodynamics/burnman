import unittest
from util import BurnManTest

import burnman
from burnman.tools.eos import check_eos_consistency


class EosConsistency(BurnManTest):
    def test_HP(self):
        P = 10.0e9
        T = 3000.0
        per = burnman.minerals.HP_2011_ds62.per()
        self.assertEqual(
            check_eos_consistency(
                per,
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_HP98(self):
        P = 10.0e9
        T = 3000.0

        params = burnman.minerals.HP_2011_ds62.per().params
        params["equation_of_state"] = "hp98"
        params["dKdT_0"] = -1.0e7
        per = burnman.Mineral(params)
        self.assertEqual(
            check_eos_consistency(
                per,
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_HP_liquid(self):
        P = 10.0e9
        T = 3000.0
        self.assertEqual(
            check_eos_consistency(
                burnman.minerals.HP_2011_ds62.foL(),
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_CORK(self):
        P = 10.0e9
        T = 3000.0
        self.assertEqual(
            check_eos_consistency(
                burnman.minerals.HP_2011_fluids.CO2(),
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_SLB(self):
        P = 10.0e9
        T = 3000.0
        self.assertEqual(
            check_eos_consistency(burnman.minerals.SLB_2011.periclase(), P, T), True
        )

    def test_SLB_2022(self):
        P = 2.0e9
        T = 500.0
        self.assertEqual(
            check_eos_consistency(burnman.minerals.SLB_2022.almandine(), P, T), True
        )

    def test_DKS_liquid(self):
        P = 10.0e9
        T = 3000.0
        self.assertEqual(
            check_eos_consistency(
                burnman.minerals.DKS_2013_liquids.Mg3Si2O7_liquid(),
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_DKS_solid(self):
        P = 10.0e9
        T = 3000.0
        params = burnman.minerals.DKS_2013_solids.periclase().params
        params2 = burnman.minerals.SLB_2011.periclase().params
        for p in ["G_0", "Gprime_0", "eta_s_0"]:
            params[p] = params2[p]
        m = burnman.Mineral(params)
        self.assertEqual(
            check_eos_consistency(m, P, T, including_shear_properties=True),
            True,
        )

    def test_modifier(self):
        P = 10.0e9
        T = 3000.0
        self.assertEqual(
            check_eos_consistency(
                burnman.minerals.Sundman_1991.bcc_iron(),
                P,
                T,
                including_shear_properties=False,
            ),
            True,
        )

    def test_solution(self):
        P = 10.0e9
        T = 3000.0
        m = burnman.minerals.SLB_2011.garnet(molar_fractions=[0.2, 0.2, 0.2, 0.2, 0.2])
        self.assertEqual(check_eos_consistency(m, P, T), True)


if __name__ == "__main__":
    unittest.main()
