import unittest
from util import BurnManTest

from burnman.eos.anharmonic_debye_pade import AnharmonicDebyePade
from burnman.eos.debye_temperature_models import SLB as theta_SLB


params = {
    "a_anh": 1500.0,
    "V_0": 2.0e-6,
    "m_anh": 3.0,
    "Debye_0": 1000.0,
    "grueneisen_0": 1.2,
    "q_0": 1.1,
    "debye_temperature_model": theta_SLB(),
    "T_0": 300.0,
}


class Anharmonic(BurnManTest):
    def test_derivatives(self):
        V_0 = params["V_0"]
        V = V_0 * 0.9
        T = 2000.0

        S = AnharmonicDebyePade.entropy(T, V, params)
        C_v = AnharmonicDebyePade.heat_capacity_v(T, V, params)
        P = AnharmonicDebyePade.pressure(T, V, params)
        K_T = AnharmonicDebyePade.isothermal_bulk_modulus(T, V, params)
        alphaK_T = AnharmonicDebyePade.dSdV(T, V, params)

        # Test partial derivatives
        dV = 1.0e-11
        dT = 1.0e-2

        dFdV = (
            AnharmonicDebyePade.helmholtz_energy(T, V + dV, params)
            - AnharmonicDebyePade.helmholtz_energy(T, V - dV, params)
        ) / (2.0 * dV)
        dFdT = (
            AnharmonicDebyePade.helmholtz_energy(T + dT, V, params)
            - AnharmonicDebyePade.helmholtz_energy(T - dT, V, params)
        ) / (2.0 * dT)

        dSdV = (
            AnharmonicDebyePade.entropy(T, V + dV, params)
            - AnharmonicDebyePade.entropy(T, V - dV, params)
        ) / (2.0 * dV)
        dSdT = (
            AnharmonicDebyePade.entropy(T + dT, V, params)
            - AnharmonicDebyePade.entropy(T - dT, V, params)
        ) / (2.0 * dT)

        dPdV = (
            AnharmonicDebyePade.pressure(T, V + dV, params)
            - AnharmonicDebyePade.pressure(T, V - dV, params)
        ) / (2.0 * dV)
        dPdT = (
            AnharmonicDebyePade.pressure(T + dT, V, params)
            - AnharmonicDebyePade.pressure(T - dT, V, params)
        ) / (2.0 * dT)

        self.assertAlmostEqual(S, -dFdT, places=6)
        self.assertAlmostEqual(C_v, T * dSdT, places=6)
        self.assertAlmostEqual(P / 1.0e8, -dFdV / 1.0e8, places=6)
        self.assertAlmostEqual(K_T / 1.0e8, -V * dPdV / 1.0e8, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dSdV / 1.0e3, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dPdT / 1.0e3, places=6)


if __name__ == "__main__":
    unittest.main()
