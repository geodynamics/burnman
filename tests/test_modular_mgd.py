import unittest
from util import BurnManTest

from burnman import Mineral, Composite
from burnman.minerals.SLB_2011 import forsterite as forsterite_slb
from burnman.utils.chemistry import dictionarize_formula, formula_mass
from burnman.eos.helper import create
from burnman.tools.eos import check_eos_consistency
from burnman.eos.debye_temperature_models import SLB as theta_SLB
from burnman.eos.debye_temperature_models import PowerLawGamma as PLG
from burnman.eos.debye_temperature_models import PowerLawGammaSimple as PLGS


formula = "Mg2SiO4"
formula = dictionarize_formula(formula)
mineral_params = {
    "name": "Forsterite",
    "formula": formula,
    "equation_of_state": "modular_mgd",
    "F_0": -2055403.0,
    "V_0": 4.3603e-05,
    "K_0": 1.279555e11,
    "Kprime_0": 4.21796,
    "T_0": 298.15,
    "n": sum(formula.values()),
    "molar_mass": formula_mass(formula),
}


class ModularMGD(BurnManTest):
    def test_setup(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1

        _ = Mineral(params)

    def test_SLB_derivatives(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1

        t = theta_SLB()

        dVrel = 1.0e-6
        Vrel = 0.7
        dtdVrel_numerical = (t(Vrel + dVrel, params) - t(Vrel - dVrel, params)) / (
            2.0 * dVrel
        )
        dtdVrel_analytical = t.dVrel(Vrel, params)
        d2tdVrel2_numerical = (
            t.dVrel(Vrel + dVrel, params) - t.dVrel(Vrel - dVrel, params)
        ) / (2.0 * dVrel)
        d2tdVrel2_analytical = t.d2dVrel2(Vrel, params)
        self.assertAlmostEqual(dtdVrel_numerical, dtdVrel_analytical, places=6)
        self.assertAlmostEqual(d2tdVrel2_numerical, d2tdVrel2_analytical, places=6)

    def test_PLG_derivatives(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLG()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.0
        params["c_1"] = 0.2
        params["c_2"] = 0.1
        params["q_1"] = 1.0
        params["q_2"] = 2.0

        t = PLG()

        dVrel = 1.0e-6
        Vrel = 0.7
        dtdVrel_numerical = (t(Vrel + dVrel, params) - t(Vrel - dVrel, params)) / (
            2.0 * dVrel
        )
        dtdVrel_analytical = t.dVrel(Vrel, params)
        d2tdVrel2_numerical = (
            t.dVrel(Vrel + dVrel, params) - t.dVrel(Vrel - dVrel, params)
        ) / (2.0 * dVrel)
        d2tdVrel2_analytical = t.d2dVrel2(Vrel, params)
        self.assertAlmostEqual(dtdVrel_numerical, dtdVrel_analytical, places=6)
        self.assertAlmostEqual(d2tdVrel2_numerical, d2tdVrel2_analytical, places=6)

    def test_PLGS_derivatives(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLGS()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1

        t = PLGS()

        dVrel = 1.0e-6
        Vrel = 0.7
        dtdVrel_numerical = (t(Vrel + dVrel, params) - t(Vrel - dVrel, params)) / (
            2.0 * dVrel
        )
        dtdVrel_analytical = t.dVrel(Vrel, params)
        d2tdVrel2_numerical = (
            t.dVrel(Vrel + dVrel, params) - t.dVrel(Vrel - dVrel, params)
        ) / (2.0 * dVrel)
        d2tdVrel2_analytical = t.d2dVrel2(Vrel, params)
        self.assertAlmostEqual(dtdVrel_numerical, dtdVrel_analytical, places=6)
        self.assertAlmostEqual(d2tdVrel2_numerical, d2tdVrel2_analytical, places=6)

    def test_check_SLB_consistency(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1

        m = Mineral(params)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

    def test_check_PLG_consistency(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLG()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.0
        params["c_1"] = 0.2
        params["c_2"] = 0.1
        params["q_1"] = 1.0
        params["q_2"] = 2.0

        m = Mineral(params)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

    def test_check_PLGS_consistency(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLGS()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.0
        params["q_0"] = 1.0

        m = Mineral(params)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

    def test_SLB_grueneisen(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1
        params["P_0"] = 1.0e5

        m = Mineral(params)
        m.set_state(1.0e5, 298.15)
        self.assertAlmostEqual(m.gr, params["grueneisen_0"], places=6)

    def test_PLG_grueneisen(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLG()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.0
        params["c_1"] = 0.2
        params["c_2"] = 0.1
        params["q_1"] = 1.0
        params["q_2"] = 2.0
        params["P_0"] = 1.0e5

        m = Mineral(params)
        m.set_state(1.0e5, 298.15)
        self.assertAlmostEqual(m.gr, params["grueneisen_0"], places=6)

    def test_PLGS_grueneisen(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = PLGS()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1
        params["P_0"] = 1.0e5

        m = Mineral(params)
        m.set_state(1.0e5, 298.15)
        self.assertAlmostEqual(m.gr, params["grueneisen_0"], places=6)

    def test_SLB_electronic_contribution(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1
        params["P_0"] = 1.0e5
        m_without_component = Mineral(params)

        params2 = params.copy()
        params2["bel_0"] = 0.005
        params2["gel"] = 1.5

        m = Mineral(params2)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

        m.set_state(2.0e9, 2000.0)
        m_without_component.set_state(2.0e9, 2000.0)
        self.assertFalse(m_without_component.S == m.S)

    def test_SLB_anharmonic_contribution(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1
        params["P_0"] = 1.0e5
        m_without_component = Mineral(params)

        params2 = params.copy()
        params2["a_anh"] = 1000.0
        params2["m_anh"] = 3.0

        m = Mineral(params2)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

        m.set_state(2.0e9, 2000.0)
        m_without_component.set_state(2.0e9, 2000.0)
        self.assertFalse(m_without_component.S == m.S)

    def test_SPOCK_isothermal_contribution(self):
        params = mineral_params.copy()
        params["reference_eos"] = create("spock")
        params["debye_temperature_model"] = theta_SLB()
        params["Debye_0"] = 1000.0
        params["grueneisen_0"] = 1.2
        params["q_0"] = 1.1
        params["P_0"] = 1.0e5
        params["Kprime_inf"] = 2.0
        params["Kprime_0"] = 4.21796
        params["Kdprime_0"] = -4 / 1.0e11

        m = Mineral(params)
        consistent = check_eos_consistency(
            m, 2.0e9, 2000.0, including_shear_properties=False, tol=1.0e-4
        )
        self.assertTrue(consistent)

    def test_SLB_is_same_as_MGD(self):
        fo_slb = forsterite_slb()
        params = fo_slb.params.copy()
        params["equation_of_state"] = "modular_mgd"
        params["reference_eos"] = create("bm3")
        params["debye_temperature_model"] = theta_SLB()
        fo_mgd = Mineral(params)

        minerals = Composite([fo_slb, fo_mgd])
        minerals.set_state(1.0e5, 298.15)
        self.assertAlmostEqual(fo_slb.helmholtz, fo_mgd.helmholtz, places=3)

        minerals.set_state(1.0e9, 2000.0)
        self.assertAlmostEqual(fo_slb.helmholtz, fo_mgd.helmholtz, places=3)


if __name__ == "__main__":
    unittest.main()
