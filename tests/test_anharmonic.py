import unittest
from util import BurnManTest

from burnman.eos.anharmonic_debye import AnharmonicDebye
from burnman.eos import anharmonic_thermal_models
from burnman.eos import debye_temperature_models
from burnman.eos import anharmonic_prefactor_models
from burnman.eos import debye
from burnman.constants import gas_constant

params = {
    "a_anh": 1.5,
    "V_0": 2.0e-6,
    "m_anh": 3.0,
    "Debye_0": 1000.0,
    "grueneisen_0": 1.2,
    "q_0": 1.1,
    "debye_temperature_model": debye_temperature_models.SLB(),
    "T_0": 300.0,
    "mu_anh": 0.2,
    "sigma_anh": 1.2,
}


class Anharmonic(BurnManTest):
    def test_pade_derivatives(self):
        pade = anharmonic_thermal_models.Pade()
        T = 2000.0
        S = pade.nondimensional_entropy(T, params["Debye_0"])
        C_v = pade.nondimensional_heat_capacity(T, params["Debye_0"])
        dFdTheta_ana = pade.nondimensional_dhelmholtz_dTheta(T, params["Debye_0"])
        dSdTheta_ana = pade.nondimensional_dentropy_dTheta(T, params["Debye_0"])
        d2FdTheta2_ana = pade.nondimensional_d2helmholtz_dTheta2(T, params["Debye_0"])

        # Test partial derivatives
        dT = 1.0e-2
        dTheta = 1.0e-2
        S_num = -(
            pade.nondimensional_helmholtz_energy(T + dT, params["Debye_0"])
            - pade.nondimensional_helmholtz_energy(T - dT, params["Debye_0"])
        ) / (2.0 * dT)
        C_v_num = (
            T
            * (
                pade.nondimensional_entropy(T + dT, params["Debye_0"])
                - pade.nondimensional_entropy(T - dT, params["Debye_0"])
            )
            / (2.0 * dT)
        )
        dFdTheta_num = (
            pade.nondimensional_helmholtz_energy(T, params["Debye_0"] + dTheta)
            - pade.nondimensional_helmholtz_energy(T, params["Debye_0"] - dTheta)
        ) / (2.0 * dTheta)
        dSdTheta_num = (
            pade.nondimensional_entropy(T, params["Debye_0"] + dTheta)
            - pade.nondimensional_entropy(T, params["Debye_0"] - dTheta)
        ) / (2.0 * dTheta)
        d2FdTheta2_num = (
            pade.nondimensional_helmholtz_energy(T, params["Debye_0"] + dTheta)
            - 2 * pade.nondimensional_helmholtz_energy(T, params["Debye_0"])
            + pade.nondimensional_helmholtz_energy(T, params["Debye_0"] - dTheta)
        ) / (dTheta**2)

        self.assertAlmostEqual(S, S_num, places=6)
        self.assertAlmostEqual(C_v, C_v_num, places=6)
        self.assertAlmostEqual(dFdTheta_ana, dFdTheta_num, places=6)
        self.assertAlmostEqual(dSdTheta_ana, dSdTheta_num, places=6)
        self.assertAlmostEqual(d2FdTheta2_ana, d2FdTheta2_num, places=6)

    def test_lognormal_derivatives(self):
        lognormal = anharmonic_thermal_models.LogNormal()
        T = 1500.0
        S = lognormal.nondimensional_entropy(T, params["Debye_0"], params)
        C_v = lognormal.nondimensional_heat_capacity(T, params["Debye_0"], params)
        dFdTheta_ana = lognormal.nondimensional_dhelmholtz_dTheta(
            T, params["Debye_0"], params
        )
        dSdTheta_ana = lognormal.nondimensional_dentropy_dTheta(
            T, params["Debye_0"], params
        )
        d2FdTheta2_ana = lognormal.nondimensional_d2helmholtz_dTheta2(
            T, params["Debye_0"], params
        )

        # Test partial derivatives
        dT = 1.0e-2
        dTheta = 1.0e-2
        S_num = -(
            lognormal.nondimensional_helmholtz_energy(T + dT, params["Debye_0"], params)
            - lognormal.nondimensional_helmholtz_energy(
                T - dT, params["Debye_0"], params
            )
        ) / (2.0 * dT)
        C_v_num = (
            T
            * (
                lognormal.nondimensional_entropy(T + dT, params["Debye_0"], params)
                - lognormal.nondimensional_entropy(T - dT, params["Debye_0"], params)
            )
            / (2.0 * dT)
        )
        dFdTheta_num = (
            lognormal.nondimensional_helmholtz_energy(
                T, params["Debye_0"] + dTheta, params
            )
            - lognormal.nondimensional_helmholtz_energy(
                T, params["Debye_0"] - dTheta, params
            )
        ) / (2.0 * dTheta)
        dSdTheta_num = (
            lognormal.nondimensional_entropy(T, params["Debye_0"] + dTheta, params)
            - lognormal.nondimensional_entropy(T, params["Debye_0"] - dTheta, params)
        ) / (2.0 * dTheta)
        d2FdTheta2_num = (
            lognormal.nondimensional_helmholtz_energy(
                T, params["Debye_0"] + dTheta, params
            )
            - 2
            * lognormal.nondimensional_helmholtz_energy(T, params["Debye_0"], params)
            + lognormal.nondimensional_helmholtz_energy(
                T, params["Debye_0"] - dTheta, params
            )
        ) / (dTheta**2)

        self.assertAlmostEqual(S, S_num, places=6)
        self.assertAlmostEqual(C_v, C_v_num, places=6)
        self.assertAlmostEqual(dFdTheta_ana, dFdTheta_num, places=6)
        self.assertAlmostEqual(dSdTheta_ana, dSdTheta_num, places=6)
        self.assertAlmostEqual(d2FdTheta2_ana, d2FdTheta2_num, places=6)

    def test_derivatives_lognormal(self):
        thermal_model = anharmonic_thermal_models.LogNormal()
        prefactor_model = anharmonic_prefactor_models.PowerLaw()
        params["anharmonic_thermal_model"] = thermal_model
        params["anharmonic_prefactor_model"] = prefactor_model
        V_0 = params["V_0"]
        V = V_0 * 0.9
        T = 1500.0

        S = AnharmonicDebye.entropy(T, V, params)
        C_v = AnharmonicDebye.heat_capacity_v(T, V, params)
        P = AnharmonicDebye.pressure(T, V, params)
        K_T = AnharmonicDebye.isothermal_bulk_modulus(T, V, params)
        alphaK_T = AnharmonicDebye.dSdV(T, V, params)

        # Test partial derivatives
        dV = 1.0e-11
        dT = 1.0e-2

        dFdV = (
            AnharmonicDebye.helmholtz_energy(T, V + dV, params)
            - AnharmonicDebye.helmholtz_energy(T, V - dV, params)
        ) / (2.0 * dV)
        dFdT = (
            AnharmonicDebye.helmholtz_energy(T + dT, V, params)
            - AnharmonicDebye.helmholtz_energy(T - dT, V, params)
        ) / (2.0 * dT)

        dSdV = (
            AnharmonicDebye.entropy(T, V + dV, params)
            - AnharmonicDebye.entropy(T, V - dV, params)
        ) / (2.0 * dV)
        dSdT = (
            AnharmonicDebye.entropy(T + dT, V, params)
            - AnharmonicDebye.entropy(T - dT, V, params)
        ) / (2.0 * dT)

        dPdV = (
            AnharmonicDebye.pressure(T, V + dV, params)
            - AnharmonicDebye.pressure(T, V - dV, params)
        ) / (2.0 * dV)
        dPdT = (
            AnharmonicDebye.pressure(T + dT, V, params)
            - AnharmonicDebye.pressure(T - dT, V, params)
        ) / (2.0 * dT)

        self.assertAlmostEqual(S, -dFdT, places=6)
        self.assertAlmostEqual(C_v, T * dSdT, places=6)
        self.assertAlmostEqual(P / 1.0e8, -dFdV / 1.0e8, places=6)
        self.assertAlmostEqual(K_T / 1.0e8, -V * dPdV / 1.0e8, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dSdV / 1.0e3, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dPdT / 1.0e3, places=6)

    def test_derivatives_pade(self):
        thermal_model = anharmonic_thermal_models.Pade()
        prefactor_model = anharmonic_prefactor_models.PowerLaw()
        params["anharmonic_thermal_model"] = thermal_model
        params["anharmonic_prefactor_model"] = prefactor_model
        V_0 = params["V_0"]
        V = V_0 * 0.9
        T = 1500.0

        S = AnharmonicDebye.entropy(T, V, params)
        C_v = AnharmonicDebye.heat_capacity_v(T, V, params)
        P = AnharmonicDebye.pressure(T, V, params)
        K_T = AnharmonicDebye.isothermal_bulk_modulus(T, V, params)
        alphaK_T = AnharmonicDebye.dSdV(T, V, params)

        # Test partial derivatives
        dV = 1.0e-11
        dT = 1.0e-2

        dFdV = (
            AnharmonicDebye.helmholtz_energy(T, V + dV, params)
            - AnharmonicDebye.helmholtz_energy(T, V - dV, params)
        ) / (2.0 * dV)
        dFdT = (
            AnharmonicDebye.helmholtz_energy(T + dT, V, params)
            - AnharmonicDebye.helmholtz_energy(T - dT, V, params)
        ) / (2.0 * dT)

        dSdV = (
            AnharmonicDebye.entropy(T, V + dV, params)
            - AnharmonicDebye.entropy(T, V - dV, params)
        ) / (2.0 * dV)
        dSdT = (
            AnharmonicDebye.entropy(T + dT, V, params)
            - AnharmonicDebye.entropy(T - dT, V, params)
        ) / (2.0 * dT)

        dPdV = (
            AnharmonicDebye.pressure(T, V + dV, params)
            - AnharmonicDebye.pressure(T, V - dV, params)
        ) / (2.0 * dV)
        dPdT = (
            AnharmonicDebye.pressure(T + dT, V, params)
            - AnharmonicDebye.pressure(T - dT, V, params)
        ) / (2.0 * dT)

        self.assertAlmostEqual(S, -dFdT, places=6)
        self.assertAlmostEqual(C_v, T * dSdT, places=6)
        self.assertAlmostEqual(P / 1.0e8, -dFdV / 1.0e8, places=6)
        self.assertAlmostEqual(K_T / 1.0e8, -V * dPdV / 1.0e8, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dSdV / 1.0e3, places=6)
        self.assertAlmostEqual(alphaK_T / 1.0e3, dPdT / 1.0e3, places=6)

    def test_standard_entropy_pade(self):
        thermal_model = anharmonic_thermal_models.Pade()
        prefactor_model = anharmonic_prefactor_models.PowerLaw()
        params["anharmonic_thermal_model"] = thermal_model
        params["anharmonic_prefactor_model"] = prefactor_model
        V = params["V_0"]
        T = params["debye_temperature_model"].value(Vrel=1.0, params=params)
        self.assertAlmostEqual(T, params["Debye_0"], places=6)

        S = AnharmonicDebye.entropy(T, V, params) / params["a_anh"]
        S1 = -debye.thermal_energy(1.0, 1.0, 1.0) / gas_constant / 3.0
        self.assertAlmostEqual(S, S1, places=4)
        self.assertAlmostEqual(S, -0.6744, places=4)

    def test_prefactor_power_law(self):
        prefactor_model = anharmonic_prefactor_models.PowerLaw()
        model_params = {"a_anh": 1.1, "m_anh": 5.0}

        # Check first and second derivatives
        # functions are called value, dVrel and dVrel2
        Vrel = 0.8
        dVrel = 1.0e-6

        dAdVrel_num = (
            prefactor_model.value(Vrel + dVrel, model_params)
            - prefactor_model.value(Vrel - dVrel, model_params)
        ) / (2.0 * dVrel)
        dAdVrel_analytical = prefactor_model.dVrel(Vrel, model_params)
        self.assertAlmostEqual(dAdVrel_num, dAdVrel_analytical, places=6)

        d2AdVrel2_num = (
            prefactor_model.dVrel(Vrel + dVrel, model_params)
            - prefactor_model.dVrel(Vrel - dVrel, model_params)
        ) / (2.0 * dVrel)
        d2AdVrel2_analytical = prefactor_model.dVrel2(Vrel, model_params)
        self.assertAlmostEqual(d2AdVrel2_num, d2AdVrel2_analytical, places=6)

    def test_prefactor_sigmoid(self):
        prefactor_model = anharmonic_prefactor_models.Sigmoid()
        model_params = {"a_anh": 1.1, "b_anh": 1.2, "c_anh": 1.3}

        # Check first and second derivatives
        # functions are called value, dVrel and dVrel2
        Vrel = 0.8
        dVrel = 1.0e-6

        dAdVrel_num = (
            prefactor_model.value(Vrel + dVrel, model_params)
            - prefactor_model.value(Vrel - dVrel, model_params)
        ) / (2.0 * dVrel)
        dAdVrel_analytical = prefactor_model.dVrel(Vrel, model_params)
        self.assertAlmostEqual(dAdVrel_num, dAdVrel_analytical, places=6)

        d2AdVrel2_num = (
            prefactor_model.dVrel(Vrel + dVrel, model_params)
            - prefactor_model.dVrel(Vrel - dVrel, model_params)
        ) / (2.0 * dVrel)
        d2AdVrel2_analytical = prefactor_model.dVrel2(Vrel, model_params)
        self.assertAlmostEqual(d2AdVrel2_num, d2AdVrel2_analytical, places=6)


if __name__ == "__main__":
    unittest.main()
