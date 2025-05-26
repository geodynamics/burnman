from util import BurnManTest
import unittest
import numpy as np

import burnman
from burnman.utils.chemistry import formula_to_string, sum_formulae
from burnman.minerals.SLB_2011 import pyrope, grossular, enstatite
from burnman.classes.elasticsolutionmodel import (
    ElasticIdealSolution,
    ElasticSymmetricRegularSolution,
    ElasticAsymmetricRegularSolution,
    ElasticSubregularSolution,
    ElasticFunctionSolution,
)


class pyrope_ss(burnman.ElasticSolution):
    # One-mineral solid solution
    def __init__(self, molar_fractions=None):
        self.name = "Dummy solid solution"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[[pyrope(), "[Mg]2SiO4"]], energy_interaction=[]
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class pyrope_pyrope_ss(burnman.ElasticSolution):
    # Two-mineral solid solution
    def __init__(self, molar_fractions=None):
        self.name = "Fo-Fo solid solution"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[[pyrope(), "[Mg]2SiO4"], [pyrope(), "[Mg]2SiO4"]],
            energy_interaction=[[0.0]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class two_identical_gt_ss(burnman.ElasticSolution):
    # Two-mineral solid solution
    def __init__(self, molar_fractions=None):
        self.name = "Fo-Fo solid solution"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[[pyrope(), "[Mg]2SiO4"], [pyrope(), "[Fe]2SiO4"]],
            energy_interaction=[[0.0]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class garnet_ideal_ss(burnman.ElasticSolution):
    # Ideal solid solution
    def __init__(self, molar_fractions=None):
        self.name = "Fo-Fo solid solution"
        self.solution_model = ElasticIdealSolution(
            endmembers=[[pyrope(), "[Mg]2SiO4"], [grossular(), "[Fe]2SiO4"]]
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class garnet_ss(burnman.ElasticSolution):
    # garnet solid solution
    def __init__(self, molar_fractions=None):
        self.name = "garnet"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[[pyrope(), "[Mg]2SiO4"], [grossular(), "[Fe]2SiO4"]],
            energy_interaction=[[8.4e3]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


def garnet_helmholtz_function(volume, temperature, molar_amounts):
    n_moles = sum(molar_amounts)
    molar_fractions = molar_amounts / n_moles
    return n_moles * (
        (8.4e3 + 10.0e9 * volume) * molar_fractions[0] * molar_fractions[1]
    )


class garnet_ss_function(burnman.ElasticSolution):
    # Three-endmember, two site symmetric solid solution
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss"
        self.solution_model = ElasticFunctionSolution(
            endmembers=[
                [pyrope(), "[Mg]3[Al]2Si3O12"],
                [grossular(), "[Ca]3[Al]2Si3O12"],
            ],
            excess_helmholtz_function=garnet_helmholtz_function,
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class orthopyroxene(burnman.ElasticSolution):
    # Orthopyroxene solid solution
    def __init__(self, molar_fractions=None):
        self.name = "orthopyroxene"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[
                [enstatite(), "[Mg][Mg]Si2O6"],
                [enstatite(), "[Mg1/2Al1/2][Mg1/2Al1/2]AlSiO6"],
            ],
            energy_interaction=[[burnman.constants.gas_constant * 1.0e3]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class two_site_ss(burnman.ElasticSolution):
    # Three-endmember, two site symmetric solid solution
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss"
        self.solution_model = ElasticSymmetricRegularSolution(
            endmembers=[
                [pyrope(), "[Mg]3[Al]2Si3O12"],
                [pyrope(), "[Fe]3[Al]2Si3O12"],
                [pyrope(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[[10.0e3, 5.0e3], [-10.0e3]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


def ss_helmholtz_function(volume, temperature, molar_amounts):
    n_moles = sum(molar_amounts)
    molar_fractions = molar_amounts / n_moles
    return n_moles * (
        10.0e3 * molar_fractions[0] * molar_fractions[1]
        + 5.0e3 * molar_fractions[0] * molar_fractions[2]
        + -10.0e3 * molar_fractions[1] * molar_fractions[2]
    )


class two_site_ss_function(burnman.ElasticSolution):
    # Three-endmember, two site symmetric solid solution
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss"
        self.solution_model = ElasticFunctionSolution(
            endmembers=[
                [pyrope(), "[Mg]3[Al]2Si3O12"],
                [pyrope(), "[Fe]3[Al]2Si3O12"],
                [pyrope(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            excess_helmholtz_function=ss_helmholtz_function,
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class two_site_ss_subregular(burnman.ElasticSolution):
    # Three-endmember, two site solid solution
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular symmetric)"
        self.solution_model = ElasticSubregularSolution(
            endmembers=[
                [pyrope(), "[Mg]3[Al]2Si3O12"],
                [pyrope(), "[Fe]3[Al]2Si3O12"],
                [pyrope(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[
                [[10.0e3, 10.0e3], [5.0e3, 5.0e3]],
                [[-10.0e3, -10.0e3]],
            ],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class two_site_ss_subregular_ternary(burnman.ElasticSolution):
    # Three-endmember, two site solid solution with ternary term
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular symmetric)"
        self.solution_model = ElasticSubregularSolution(
            endmembers=[
                [pyrope(), "[Mg]3[Al]2Si3O12"],
                [grossular(), "[Fe]3[Al]2Si3O12"],
                [pyrope(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[
                [[10.0e3, -10.0e3], [5.0e3, 3.0e3]],
                [[-10.0e3, -10.0e3]],
            ],
            pressure_interaction=[[[1.0e9, -1.0e9], [0.0, 1.0e9]], [[0.0, 0.0]]],
            entropy_interaction=[[[1.0, -2.0], [0.0, 1.0]], [[0.0, 0.0]]],
            energy_ternary_terms=[[0, 1, 2, 3.0e3]],
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class temkin_ss(burnman.ElasticSolution):
    # Temkin solution
    def __init__(self, molar_fractions=None):
        self.name = "ol-q-water melt (HGP 2018)"
        self.solution_model = ElasticIdealSolution(
            endmembers=[
                [pyrope(), "[Mg]4[Sitet]1[Vac]2"],
                [grossular(), "[Fe]4[Sitet]1[Vac]2"],
                [pyrope(), "[]0[Sinet]1[Vac]2"],
                [grossular(), "[]0[]0[H]2"],
            ]
        )

        burnman.ElasticSolution.__init__(self, molar_fractions)


class test_ElasticSolution(BurnManTest):
    def setup_1min_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = pyrope()
        fo.set_state(P, T)
        py_ss = pyrope_ss()
        py_ss.set_composition([1.0])
        py_ss.set_state(P, T)
        return fo, py_ss

    def setup_2min_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = pyrope()
        fo.set_state(P, T)
        py_py_ss = pyrope_pyrope_ss()
        py_py_ss.set_composition([0.3, 0.7])
        py_py_ss.set_state(P, T)
        return fo, py_py_ss

    def setup_gt_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = pyrope()
        fo.set_state(P, T)

        gt_ss = garnet_ss()
        gt_ss.set_composition([1.0, 0.0])
        gt_ss.set_state(P, T)
        return fo, gt_ss

    def test_1_gibbs(self):
        fo, py_ss = self.setup_1min_ss()
        endmember_properties = [
            fo.gibbs,
            fo.H,
            fo.S,
            fo.V,
            fo.C_p,
            fo.C_v,
            fo.alpha,
            fo.K_T,
            fo.K_S,
            fo.gr,
            fo.G,
        ]
        ss_properties = [
            py_ss.gibbs,
            py_ss.H,
            py_ss.S,
            py_ss.V,
            py_ss.C_p,
            py_ss.C_v,
            py_ss.alpha,
            py_ss.K_T,
            py_ss.K_S,
            py_ss.gr,
            py_ss.G,
        ]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_2_gibbs(self):
        fo, py_ss = self.setup_2min_ss()
        endmember_properties = [
            fo.gibbs,
            fo.H,
            fo.S,
            fo.V,
            fo.C_p,
            fo.C_v,
            fo.alpha,
            fo.K_T,
            fo.K_S,
            fo.gr,
            fo.G,
        ]
        ss_properties = [
            py_ss.gibbs,
            py_ss.H,
            py_ss.S,
            py_ss.V,
            py_ss.C_p,
            py_ss.C_v,
            py_ss.alpha,
            py_ss.K_T,
            py_ss.K_S,
            py_ss.gr,
            py_ss.G,
        ]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_gt_gibbs(self):
        fo, py_ss = self.setup_gt_ss()
        endmember_properties = [
            fo.gibbs,
            fo.H,
            fo.S,
            fo.V,
            fo.C_p,
            fo.C_v,
            fo.alpha,
            fo.K_T,
            fo.K_S,
            fo.gr,
        ]
        ss_properties = [
            py_ss.gibbs,
            py_ss.H,
            py_ss.S,
            py_ss.V,
            py_ss.C_p,
            py_ss.C_v,
            py_ss.alpha,
            py_ss.K_T,
            py_ss.K_S,
            py_ss.gr,
        ]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_set_state_with_volume(self):
        m = garnet_ss()
        m.set_composition([0.5, 0.5])

        P0 = 6.0e9
        T0 = 1000.0
        T1 = 298.15

        m.set_state(P0, T0)
        V = m.V
        m.set_state_with_volume(V, T1)
        P1 = m.pressure
        m.set_state(P0, T0)  # forget new state
        m.set_state(P1, T1)  # return to new state
        self.assertFloatEqual(V, m.V)

    def test_gt_Wh(self):
        gt_ss = garnet_ss()
        H_excess = gt_ss.solution_model.excess_enthalpy(
            1.0e5, 1000.0, [0.5, 0.5]
        )  # Hxs = Exs if Vxs=0
        We = gt_ss.solution_model.We[0][1]
        self.assertArraysAlmostEqual([We / 4.0], [H_excess])

    def test_order_disorder(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.0e5, 300.0)

        en = enstatite()
        en.set_state(1.0e5, 300.0)

        self.assertArraysAlmostEqual([opx.gibbs], [en.gibbs])

    def test_site_totals(self):
        ss = two_site_ss()
        ss.set_composition([0.3, 0.3, 0.4])
        ss.set_state(1.0e5, 300.0)

        site_fractions = np.dot(
            ss.molar_fractions, ss.solution_model.endmember_occupancies
        )
        i = 0
        site_fill = []
        ones = [1.0] * ss.solution_model.n_sites
        for site in ss.solution_model.sites:
            site_fill.append(sum(site_fractions[i : i + len(site)]))
            i += len(site)

        self.assertArraysAlmostEqual(site_fill, ones)

    def test_molar_mass(self):
        ss = garnet_ss()
        ss.set_composition(np.array([0.5, 0.5]))
        self.assertArraysAlmostEqual(
            [ss.molar_mass],
            [
                0.5 * pyrope().params["molar_mass"]
                + 0.5 * grossular().params["molar_mass"]
            ],
        )

    def test_subregular_model_ternary_volume(self):
        ss = two_site_ss_subregular_ternary()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_composition(f0)

        P = 1.0e9
        T = 1000.0
        dP = 1000.0
        ss.set_state(P, T)

        V = ss.V
        ss.set_state(P + dP / 2.0, T)
        G1 = ss.gibbs
        ss.set_state(P - dP / 2.0, T)
        G0 = ss.gibbs

        self.assertFloatEqual(V, (G1 - G0) / dP)

    def test_subregular_model_ternary_partial_gibbs_multicomponent_change(self):
        ss = two_site_ss_subregular_ternary()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_composition(f0)
        ss.set_state(1.0e9, 1000.0)

        dGdx = ss.partial_gibbs

        df = 0.0001
        dGdx2 = np.empty(3)
        for i, f_mod in enumerate(np.eye(3) * df):
            ss.set_composition((f0 - f_mod / 2) / (1.0 - df / 2.0))
            G0 = ss.gibbs * (1.0 - df / 2.0)
            ss.set_composition((f0 + f_mod / 2) / (1.0 + df / 2.0))
            G1 = ss.gibbs * (1.0 + df / 2.0)

            dGdx2[i] = (G1 - G0) / df

        self.assertArraysAlmostEqual(dGdx, dGdx2)

    def test_subregular_model_ternary_partial_entropy_multicomponent_change(self):
        ss = two_site_ss_subregular_ternary()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_composition(f0)
        ss.set_state(1.0e9, 1000.0)

        dSdx = ss.partial_entropies

        df = 0.0001
        dSdx2 = np.empty(3)
        for i, f_mod in enumerate(np.eye(3) * df):
            ss.set_composition((f0 - f_mod / 2) / (1.0 - df / 2.0))
            S0 = ss.S * (1.0 - df / 2.0)
            ss.set_composition((f0 + f_mod / 2) / (1.0 + df / 2.0))
            S1 = ss.S * (1.0 + df / 2.0)

            dSdx2[i] = (S1 - S0) / df

        self.assertArraysAlmostEqual(dSdx, dSdx2)

    def test_subregular_model_ternary_partial_volume_multicomponent_change(self):
        ss = two_site_ss_subregular_ternary()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_composition(f0)
        ss.set_state(1.0e9, 1000.0)

        dVdx = ss.partial_volumes

        df = 0.0001
        dVdx2 = np.empty(3)
        for i, f_mod in enumerate(np.eye(3) * df):
            ss.set_composition((f0 - f_mod / 2) / (1.0 - df / 2.0))
            V0 = ss.V * (1.0 - df / 2.0)
            ss.set_composition((f0 + f_mod / 2) / (1.0 + df / 2.0))
            V1 = ss.V * (1.0 + df / 2.0)

            dVdx2[i] = (V1 - V0) / df

        self.assertArraysAlmostEqual(dVdx, dVdx2)

    def test_subregular(self):
        ss0 = two_site_ss()
        ss1 = two_site_ss_subregular()

        ss0.set_composition([0.3, 0.3, 0.4])
        ss0.set_state(1.0e5, 300.0)

        ss1.set_composition([0.3, 0.3, 0.4])
        ss1.set_state(1.0e5, 300.0)

        self.assertArraysAlmostEqual(ss0.partial_gibbs, ss1.partial_gibbs)

    def test_activities_ideal(self):
        gt = two_identical_gt_ss()
        gt.set_composition(np.array([0.5, 0.5]))
        gt.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(gt.activities, [0.25, 0.25])

    def test_activity_coefficients_ideal(self):
        gt = two_identical_gt_ss()
        gt.set_composition(np.array([0.5, 0.5]))
        gt.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(gt.activity_coefficients, [1.0, 1.0])

    def test_activity_coefficients_non_ideal(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(opx.activity_coefficients, [np.exp(1.0), 1.0])

    def test_formula(self):
        gt = garnet_ideal_ss()
        gt.set_composition([0.5, 0.5])
        self.assertEqual(formula_to_string(gt.formula), "Ca3/2Mg3/2Al2Si3O12")
        self.assertEqual(
            formula_to_string(sum_formulae(gt.endmember_formulae, [0.5, 0.5])),
            "Ca3/2Mg3/2Al2Si3O12",
        )

    def test_stoichiometric_matrix_binary_solution(self):
        olivine = burnman.minerals.SLB_2011.mg_fe_olivine()

        self.assertTrue(len(olivine.endmember_names) == 2)
        self.assertTrue(len(olivine.elements) == 4)
        self.assertTrue(olivine.stoichiometric_matrix.shape == (2, 4))
        self.assertTrue(olivine.reaction_basis.shape[0] == 0)
        self.assertArraysAlmostEqual(
            olivine.compositional_null_basis[0], [-1.0 / 2.0, -1.0 / 2.0, 1, 0]
        )
        self.assertArraysAlmostEqual(
            olivine.compositional_null_basis[1], [-2.0, -2.0, 0, 1]
        )
        self.assertArraysAlmostEqual(olivine.independent_element_indices, [0, 1])
        self.assertArraysAlmostEqual(olivine.dependent_element_indices, [2, 3])
        self.assertTrue(olivine.n_reactions == 0)

    def test_stoichiometric_matrix_complex_solution(self):
        opx = burnman.minerals.JH_2015.orthopyroxene()

        self.assertTrue(len(opx.endmember_names) == 7)
        self.assertTrue(len(opx.elements) == 7)
        self.assertTrue(opx.stoichiometric_matrix.shape == (7, 7))
        self.assertTrue(opx.reaction_basis.shape[0] == 1)
        self.assertArraysAlmostEqual(
            opx.compositional_null_basis[0],
            [-3.0 / 2.0, -3.0 / 2.0, -3.0 / 2.0, -3.0 / 2.0, -3.0 / 2.0, -3.0 / 2.0, 1],
        )
        self.assertArraysAlmostEqual(
            opx.independent_element_indices, [0, 1, 2, 3, 4, 5]
        )
        self.assertArraysAlmostEqual(opx.dependent_element_indices, [6])
        self.assertTrue(opx.n_reactions == 1)
        self.assertArraysAlmostEqual(
            opx.reaction_basis[0], [-1.0 / 2.0, -1.0 / 2.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        )

    def test_temkin_entropies(self):
        ss = temkin_ss()
        # '[Mg]4[Sitet]1[Vac]2'],
        # '[Fe]4[Sitet]1[Vac]2'],
        # '[]0[Sinet]1[Vac]2'],
        # '[]0[]0[H]2']]
        R = burnman.constants.gas_constant
        ss.set_state(1.0e9, 600.0)

        S = np.empty((4, 2))

        ss.set_composition([0.5, 0.5, 0.0, 0.0])
        S[0] = [ss.excess_entropy, -4.0 * R * np.log(0.5)]

        ss.set_composition([0.5, 0.0, 0.5, 0.0])
        S[1] = [ss.excess_entropy, -R * np.log(0.5)]

        ss.set_composition([0.25, 0.25, 0.5, 0.0])
        S[2] = [ss.excess_entropy, -3.0 * R * np.log(0.5)]

        ss.set_composition([0.0, 0.0, 0.5, 0.5])
        S[3] = [ss.excess_entropy, -2.0 * R * np.log(0.5)]

        self.assertArraysAlmostEqual(S[:, 0], S[:, 1])

    def test_temkin_partial_gibbs(self):
        ss = temkin_ss()
        f0 = np.array([0.25, 0.35, 0.3, 0.1])
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)

        dGdx = ss.partial_gibbs

        df = 0.0001
        dGdx2 = np.empty(4)
        for i, f_mod in enumerate(np.eye(4) * df):
            ss.set_composition((f0 - f_mod / 2) / (1.0 - df / 2.0))
            G0 = ss.gibbs * (1.0 - df / 2.0)
            ss.set_composition((f0 + f_mod / 2) / (1.0 + df / 2.0))
            G1 = ss.gibbs * (1.0 + df / 2.0)

            dGdx2[i] = (G1 - G0) / df

        self.assertArraysAlmostEqual(dGdx, dGdx2)

    def test_function_solution(self):
        ss = [two_site_ss(), two_site_ss_function()]
        for s in ss:
            s.set_state(1.0e5, 300.0)
            s.set_composition([0.2, 0.3, 0.5])

        self.assertArraysAlmostEqual(ss[0].partial_gibbs, ss[1].partial_gibbs)
        self.assertArraysAlmostEqual(ss[0].partial_entropies, ss[1].partial_entropies)
        self.assertArraysAlmostEqual(ss[0].partial_volumes, ss[1].partial_volumes)
        self.assertArraysAlmostEqual(ss[0].activities, ss[1].activities)
        self.assertArraysAlmostEqual(
            ss[0].activity_coefficients, ss[1].activity_coefficients
        )

    def test_function_solution_low_proportions(self):
        ss = [two_site_ss(), two_site_ss_function()]
        for s in ss:
            s.set_state(1.0e5, 300.0)
            s.set_composition([0.5, 0.0, 0.5])

        self.assertArraysAlmostEqual(ss[0].partial_gibbs, ss[1].partial_gibbs)
        self.assertArraysAlmostEqual(ss[0].partial_entropies, ss[1].partial_entropies)
        self.assertArraysAlmostEqual(ss[0].partial_volumes, ss[1].partial_volumes)
        self.assertArraysAlmostEqual(ss[0].activities, ss[1].activities)
        self.assertArraysAlmostEqual(
            ss[0].activity_coefficients, ss[1].activity_coefficients
        )

    def test_function_solution_p_t_derivatives(self):
        ss = garnet_ss_function()
        f0 = np.array([0.6, 0.4])
        ss.set_composition(f0)

        P = 1.0e9
        T = 1000.0
        ss.set_state(P, T)

        V = ss.V
        S = ss.S

        dP = 1000.0
        ss.set_state(P + dP / 2.0, T)
        G1 = ss.gibbs
        ss.set_state(P - dP / 2.0, T)
        G0 = ss.gibbs
        self.assertFloatEqual(V, (G1 - G0) / dP)

        dT = 1.0
        ss.set_state(P, T + dT / 2.0)
        G1 = ss.gibbs
        ss.set_state(P, T - dT / 2.0)
        G0 = ss.gibbs
        self.assertFloatEqual(S, -(G1 - G0) / dT)


if __name__ == "__main__":
    unittest.main()
