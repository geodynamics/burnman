from __future__ import absolute_import

from util import BurnManTest
import unittest
import numpy as np
import warnings

import burnman
from burnman import Mineral
from burnman import CombinedMineral
from burnman.classes.solutionmodel import IdealSolution
from burnman.classes.solutionmodel import SymmetricRegularSolution
from burnman.classes.solutionmodel import SubregularSolution
from burnman.classes.solutionmodel import AsymmetricRegularSolution
from burnman.classes.solutionmodel import FunctionSolution
from burnman.classes.solutionmodel import PolynomialSolution
from burnman.utils.chemistry import dictionarize_formula, formula_mass
from burnman.utils.chemistry import formula_to_string, sum_formulae
from burnman.minerals import HGP_2018_ds633
from burnman.minerals.SLB_2011 import mg_post_perovskite
from burnman.minerals.SLB_2011 import fe_post_perovskite
from burnman.minerals.SLB_2011 import al_post_perovskite


class forsterite(Mineral):
    def __init__(self):
        formula = "Mg2.0Si1.0O4.0"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "fo",
            "formula": formula,
            "equation_of_state": "hp_tmt",
            "H_0": -2172590.0,
            "S_0": 95.1,
            "V_0": 4.366e-05,
            "Cp": [233.3, 0.001494, -603800.0, -1869.7],
            "a_0": 2.85e-05,
            "K_0": 1.285e11,
            "Kprime_0": 3.84,
            "Kdprime_0": -3e-11,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }
        Mineral.__init__(self)


class fayalite(Mineral):
    def __init__(self):
        formula = "Fe2.0Si1.0O4.0"
        formula = dictionarize_formula(formula)
        self.params = {
            "name": "fa",
            "formula": formula,
            "equation_of_state": "hp_tmt",
            "H_0": -1477720.0,
            "S_0": 151.0,
            "V_0": 4.631e-05,
            "Cp": [201.1, 0.01733, -1960600.0, -900.9],
            "a_0": 2.82e-05,
            "K_0": 1.256e11,
            "Kprime_0": 4.68,
            "Kdprime_0": -3.7e-11,
            "n": sum(formula.values()),
            "molar_mass": formula_mass(formula),
        }
        Mineral.__init__(self)


made_forsterite = CombinedMineral([forsterite(), forsterite()], [0.5, 0.5])


# One-mineral solid solution
class forsterite_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "Dummy solid solution"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[[forsterite(), "[Mg]2SiO4"]], energy_interaction=[]
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Two-mineral solid solution
class forsterite_forsterite_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "Fo-Fo solid solution"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[[forsterite(), "[Mg]2SiO4"], [forsterite(), "[Mg]2SiO4"]],
            energy_interaction=[[0.0]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Ideal solid solution
class olivine_ideal_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "Fo-Fo solid solution"
        self.solution_model = IdealSolution(
            endmembers=[[forsterite(), "[Mg]2SiO4"], [fayalite(), "[Fe]2SiO4"]]
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Olivine solid solution
class olivine_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "Olivine"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[[forsterite(), "[Mg]2SiO4"], [fayalite(), "[Fe]2SiO4"]],
            energy_interaction=[[8.4e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Olivine solid solution with combined endmember
class olivine_ss2(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "Olivine"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[[made_forsterite, "[Mg]2SiO4"], [fayalite(), "[Fe]2SiO4"]],
            energy_interaction=[[8.4e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Orthopyroxene solid solution
class orthopyroxene(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        # Name
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [forsterite(), "[Mg][Mg]Si2O6"],
                [forsterite(), "[Mg1/2Al1/2][Mg1/2Al1/2]AlSiO6"],
            ],
            energy_interaction=[[burnman.constants.gas_constant * 1.0e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Three-endmember, two site symmetric solid solution
class two_site_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[[10.0e3, 5.0e3], [-10.0e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Three-endmember, two site asymmetric solid solution
class two_site_ss_asymmetric(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (asymmetric)"
        self.solution_model = AsymmetricRegularSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            alphas=[1.0, 2.0, 2.0],
            energy_interaction=[[10.0e3, 5.0e3], [-10.0e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Three-endmember, two site solid solution
class two_site_ss_subregular(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular symmetric)"
        self.solution_model = SubregularSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[
                [[10.0e3, 10.0e3], [5.0e3, 5.0e3]],
                [[-10.0e3, -10.0e3]],
            ],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Three-endmember, two site solid solution
class two_site_ss_subregular_asymmetric(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular symmetric)"
        self.solution_model = SubregularSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[
                [[10.0e3, -10.0e3], [5.0e3, 3.0e3]],
                [[-10.0e3, -10.0e3]],
            ],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Three-endmember, two site solid solution with ternary term
class two_site_ss_subregular_ternary(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular)"
        self.solution_model = SubregularSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[
                [[10.0e3, -10.0e3], [5.0e3, 3.0e3]],
                [[-10.0e3, -10.0e3]],
            ],
            entropy_interaction=[[[1.0, -2.0], [0.0, 1.0]], [[0.0, 0.0]]],
            energy_ternary_terms=[[0, 1, 2, 3.0e3]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# As above, but using the PolynomialSolution model
class two_site_ss_polynomial_ternary(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular from polynomial)"
        self.solution_model = PolynomialSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            ESV_interactions=[
                # Wj and Wi terms in Helffrich and Wood
                [10.0e3, 1.0, 0.0, 0, 1, 1],
                [-10.0e3, -2.0, 0.0, 0, 0, 1],
                [5.0e3, 0.0, 0.0, 0, 2, 2],
                [3.0e3, 1.0, 0.0, 0, 0, 2],
                [-10.0e3, 0.0, 0.0, 1, 2, 2],
                [-10.0e3, 0.0, 0.0, 1, 1, 2],
                # 0.5*Wk terms in Helffrich and Wood
                # Because there are only three endmembers, there
                # is only one additional term for each term above:
                [5.0e3, 0.5, 0.0, 0, 1, 2],
                [-5.0e3, -1.0, 0.0, 0, 1, 2],
                [2.5e3, 0.0, 0.0, 0, 1, 2],
                [1.5e3, 0.5, 0.0, 0, 1, 2],
                [-5.0e3, 0.0, 0.0, 0, 1, 2],
                [-5.0e3, 0.0, 0.0, 0, 1, 2],
                # Ternary term
                [3.0e3, 0.0, 0.0, 0, 1, 2],
            ],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# As above, but using a transformation matrix to change the order of
# endmembers.
# The defined interactions are a function of the transformed endmembers,
# not the ones defined by the endmembers parameter
# So in this example, p[0] for the ESV_interactions corresponds to
# the proportion of the second endmember in the endmembers list,
# p[1] to the last endmember, and p[2] to the first.
class two_site_ss_polynomial_ternary_transformed(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular from polynomial)"
        self.solution_model = PolynomialSolution(
            endmembers=[
                [forsterite(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
                [forsterite(), "[Mg]3[Al]2Si3O12"],
                [forsterite(), "[Fe]3[Al]2Si3O12"],
            ],
            ESV_interactions=[
                # Wj and Wi terms in Helffrich and Wood
                [10.0e3, 1.0, 0.0, 0, 1, 1],
                [-10.0e3, -2.0, 0.0, 0, 0, 1],
                [5.0e3, 0.0, 0.0, 0, 2, 2],
                [3.0e3, 1.0, 0.0, 0, 0, 2],
                [-10.0e3, 0.0, 0.0, 1, 2, 2],
                [-10.0e3, 0.0, 0.0, 1, 1, 2],
                # 0.5*Wk terms in Helffrich and Wood
                # Because there are only three endmembers, there
                # is only one additional term for each term above:
                [5.0e3, 0.5, 0.0, 0, 1, 2],
                [-5.0e3, -1.0, 0.0, 0, 1, 2],
                [2.5e3, 0.0, 0.0, 0, 1, 2],
                [1.5e3, 0.5, 0.0, 0, 1, 2],
                [-5.0e3, 0.0, 0.0, 0, 1, 2],
                [-5.0e3, 0.0, 0.0, 0, 1, 2],
                # Ternary term
                [3.0e3, 0.0, 0.0, 0, 1, 2],
            ],
            transformation_matrix=np.array(
                [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]
            ),
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Polynomial solution using extra endmembers
# to define the mixing terms
fo1 = forsterite()
fo2 = forsterite()
fo2.property_modifiers = [
    ["linear", {"delta_E": 1200.0, "delta_S": 5.0, "delta_V": 1.0e-6}]
]


class polynomial_mbr_mixing(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "binary mbr mixing"
        self.solution_model = PolynomialSolution(
            endmembers=[
                [forsterite(), "[Mg]2"],
                [fayalite(), "[Mg]2"],
            ],
            interaction_endmembers=[fo2, fo1],
            endmember_coefficients_and_interactions=[[1.0, -1.0, 0, 1]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


fo3 = forsterite()
fo3.property_modifiers = [
    ["linear", {"delta_E": 600.0, "delta_S": 2.5, "delta_V": 0.5e-6}]
]


class polynomial_esv_and_mbr_mixing(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "two_site_ss (subregular symmetric)"
        self.solution_model = PolynomialSolution(
            endmembers=[
                [forsterite(), "[Mg]2"],
                [fayalite(), "[Mg]2"],
            ],
            ESV_interactions=[[600.0, 2.5, 0.5e-6, 0, 1]],
            interaction_endmembers=[fo3, fo1],
            endmember_coefficients_and_interactions=[[1.0, -1.0, 0, 1]],
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


# Temkin solution
class temkin_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name = "ol-q-water melt (HGP 2018)"
        self.solution_model = IdealSolution(
            endmembers=[
                [HGP_2018_ds633.foL(), "[Mg]4[Sitet]1[Vac]2"],
                [HGP_2018_ds633.faL(), "[Fe]4[Sitet]1[Vac]2"],
                [HGP_2018_ds633.qL(), "[]0[Sinet]1[Vac]2"],
                [HGP_2018_ds633.h2oL(), "[]0[]0[H]2"],
            ]
        )

        burnman.SolidSolution.__init__(self, molar_fractions)


class ppv_symmetric(burnman.Solution):
    def __init__(self, molar_fractions=None):
        self.name = "post-perovskite/bridgmanite"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [mg_post_perovskite(), "[Mg][Si]O3"],
                [fe_post_perovskite(), "[Fe][Si]O3"],
                [al_post_perovskite(), "[Al][Al]O3"],
            ],
            energy_interaction=[[0.0, 60.0e3], [0.0]],
        )

        burnman.Solution.__init__(self, molar_fractions=molar_fractions)


def excess_gibbs_function_ppv(pressure, temperature, molar_amounts):
    n_moles = sum(molar_amounts)
    molar_fractions = molar_amounts / n_moles
    return n_moles * 60.0e3 * molar_fractions[0] * molar_fractions[2]


class ppv_function(burnman.Solution):
    def __init__(self, molar_fractions=None):
        self.name = "post-perovskite/bridgmanite"
        self.solution_model = FunctionSolution(
            endmembers=[
                [mg_post_perovskite(), "[Mg][Si]O3"],
                [fe_post_perovskite(), "[Fe][Si]O3"],
                [al_post_perovskite(), "[Al][Al]O3"],
            ],
            excess_gibbs_function=excess_gibbs_function_ppv,
        )

        burnman.Solution.__init__(self, molar_fractions=molar_fractions)


class test_solidsolution(BurnManTest):
    def setup_1min_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = forsterite()
        fo.set_state(P, T)
        fo_ss = forsterite_ss()
        fo_ss.set_composition([1.0])
        fo_ss.set_state(P, T)
        return fo, fo_ss

    def setup_2min_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = forsterite()
        fo.set_state(P, T)
        fo_fo_ss = forsterite_forsterite_ss()
        fo_fo_ss.set_composition([0.3, 0.7])
        fo_fo_ss.set_state(P, T)
        return fo, fo_fo_ss

    def setup_ol_ss(self):
        P = 1.0e5
        T = 1000.0
        fo = forsterite()
        fo.set_state(P, T)

        ol_ss = olivine_ss()
        ol_ss.set_composition([1.0, 0.0])
        ol_ss.set_state(P, T)
        return fo, ol_ss

    def setup_ol_ss2(self):
        P = 1.0e5
        T = 1000.0
        fo = forsterite()
        fo.set_state(P, T)

        ol_ss = olivine_ss2()
        ol_ss.set_composition([1.0, 0.0])
        ol_ss.set_state(P, T)
        return fo, ol_ss

    def test_1_gibbs(self):
        fo, fo_ss = self.setup_1min_ss()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
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
                fo_ss.gibbs,
                fo_ss.H,
                fo_ss.S,
                fo_ss.V,
                fo_ss.C_p,
                fo_ss.C_v,
                fo_ss.alpha,
                fo_ss.K_T,
                fo_ss.K_S,
                fo_ss.gr,
                fo_ss.G,
            ]
            assert len(w) == 3  # we expect to trigger 3 shear modulus warnings
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_2_gibbs(self):
        fo, fo_ss = self.setup_2min_ss()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
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
                fo_ss.gibbs,
                fo_ss.H,
                fo_ss.S,
                fo_ss.V,
                fo_ss.C_p,
                fo_ss.C_v,
                fo_ss.alpha,
                fo_ss.K_T,
                fo_ss.K_S,
                fo_ss.gr,
                fo_ss.G,
            ]
            assert len(w) == 4  # we expect to trigger 4 shear modulus warnings
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_gibbs(self):
        fo, fo_ss = self.setup_ol_ss()
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
            fo_ss.gibbs,
            fo_ss.H,
            fo_ss.S,
            fo_ss.V,
            fo_ss.C_p,
            fo_ss.C_v,
            fo_ss.alpha,
            fo_ss.K_T,
            fo_ss.K_S,
            fo_ss.gr,
        ]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_gibbs2(self):
        fo, fo_ss = self.setup_ol_ss2()
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
            fo_ss.gibbs,
            fo_ss.H,
            fo_ss.S,
            fo_ss.V,
            fo_ss.C_p,
            fo_ss.C_v,
            fo_ss.alpha,
            fo_ss.K_T,
            fo_ss.K_S,
            fo_ss.gr,
        ]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_set_state_with_volume(self):
        m = olivine_ss()
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

    def test_ol_Wh(self):
        ol_ss = olivine_ss()
        H_excess = ol_ss.solution_model.excess_enthalpy(
            1.0e5, 1000.0, [0.5, 0.5]
        )  # Hxs = Exs if Vxs=0
        We = ol_ss.solution_model.We[0][1]
        self.assertArraysAlmostEqual([We / 4.0], [H_excess])

    def test_order_disorder(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.0e5, 300.0)
        self.assertArraysAlmostEqual([opx.excess_gibbs], [0.0])

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

    def test_set_method(self):
        ss = olivine_ss()
        ss.set_method("hp_tmt")

    def test_molar_mass(self):
        ss = olivine_ss()
        ss.set_composition(np.array([0.5, 0.5]))
        self.assertArraysAlmostEqual(
            [ss.molar_mass],
            [
                0.5 * forsterite().params["molar_mass"]
                + 0.5 * fayalite().params["molar_mass"]
            ],
        )

    def test_asymmetric_model_hessian_one_component_change(self):
        ss = two_site_ss_asymmetric()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_state(1.0e5, 300.0)
        dp = 1.0e-6
        H = np.zeros((3, 3))
        for i in range(3):
            f = np.array(f0)
            f[i] -= dp / 2.0
            f /= np.sum(f)
            ss.set_composition(f)
            dGdx1 = ss.partial_gibbs

            f = np.array(f0)
            f[i] += dp / 2.0
            f /= np.sum(f)
            ss.set_composition(f)
            dGdx2 = ss.partial_gibbs

            H[i, :] = (dGdx2 - dGdx1) / dp

        ss.set_composition(f0)
        for i in range(3):
            self.assertArraysAlmostEqual(ss.gibbs_hessian[i], H[i])

    def test_subregular_model_hessian_one_component_change(self):
        ss = two_site_ss_subregular()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_state(1.0e5, 300.0)
        dp = 1.0e-6
        H = np.zeros((3, 3))
        for i in range(3):
            f = np.array(f0)
            f[i] -= dp / 2.0
            f /= np.sum(f)
            ss.set_composition(f)
            dGdx1 = ss.partial_gibbs

            f = np.array(f0)
            f[i] += dp / 2.0
            f /= np.sum(f)
            ss.set_composition(f)
            dGdx2 = ss.partial_gibbs

            H[i, :] = (dGdx2 - dGdx1) / dp

        ss.set_composition(f0)
        for i in range(3):
            self.assertArraysAlmostEqual(ss.gibbs_hessian[i], H[i])

    def test_asymmetric_model_hessian_multicomponent_change(self):
        ss = two_site_ss_asymmetric()
        f0 = np.array([0.25, 0.35, 0.4])
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)
        H0 = ss.gibbs_hessian

        df = np.array([2.0e-6, -1.5e-6, -0.5e-6])
        ss.set_composition(f0 - df / 2.0)
        dGdx1 = ss.partial_gibbs
        ss.set_composition(f0 + df / 2.0)
        dGdx2 = ss.partial_gibbs

        self.assertArraysAlmostEqual(H0.dot(df), dGdx2 - dGdx1)

    def test_subregular_model_hessian_multicomponent_change(self):
        ss = two_site_ss_subregular_asymmetric()
        f0 = [0.25, 0.35, 0.4]
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)
        H0 = ss.gibbs_hessian

        df = np.array([2.0e-6, -1.5e-6, -0.5e-6])
        ss.set_composition(f0 - df / 2.0)
        dGdx1 = ss.partial_gibbs
        ss.set_composition(f0 + df / 2.0)
        dGdx2 = ss.partial_gibbs

        self.assertArraysAlmostEqual(H0.dot(df), dGdx2 - dGdx1)

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

    def test_subregular_model_ternary_hessian_multicomponent_change(self):
        ss = two_site_ss_subregular_ternary()
        f0 = [0.25, 0.35, 0.4]
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)
        H0 = ss.gibbs_hessian

        df = np.array([2.0e-6, -1.5e-6, -0.5e-6])
        ss.set_composition(f0 - df / 2.0)
        dGdx1 = ss.partial_gibbs
        ss.set_composition(f0 + df / 2.0)
        dGdx2 = ss.partial_gibbs

        self.assertArraysAlmostEqual(H0.dot(df), dGdx2 - dGdx1)

    def test_polynomial_using_subregular_scalars(self):

        ss1 = two_site_ss_subregular_ternary()
        ss2 = two_site_ss_polynomial_ternary()
        f0 = [0.25, 0.35, 0.4]
        for ss in [ss1, ss2]:
            ss.set_state(1.0e5, 300.0)
            ss.set_composition(f0)

        ss1_properties = [
            ss1.gibbs,
            ss1.H,
            ss1.S,
            ss1.V,
            ss1.C_p,
            ss1.C_v,
            ss1.alpha,
            ss1.K_T,
            ss1.K_S,
            ss1.gr,
        ]
        ss2_properties = [
            ss2.gibbs,
            ss2.H,
            ss2.S,
            ss2.V,
            ss2.C_p,
            ss2.C_v,
            ss2.alpha,
            ss2.K_T,
            ss2.K_S,
            ss2.gr,
        ]
        self.assertArraysAlmostEqual(ss1_properties, ss2_properties)

    def test_polynomial_using_subregular_gibbs(self):
        ss1 = two_site_ss_subregular_ternary()
        ss2 = two_site_ss_polynomial_ternary()
        f0 = [0.25, 0.35, 0.4]
        for ss in [ss1, ss2]:
            ss.set_state(1.0e5, 300.0)
            ss.set_composition(f0)

        self.assertArraysAlmostEqual(ss1.excess_partial_gibbs, ss2.excess_partial_gibbs)
        self.assertArraysAlmostEqual(
            ss1.gibbs_hessian.flatten(), ss2.gibbs_hessian.flatten()
        )

    def test_polynomial_using_subregular_entropy(self):
        ss1 = two_site_ss_subregular_ternary()
        ss2 = two_site_ss_polynomial_ternary()
        f0 = [0.25, 0.35, 0.4]
        for ss in [ss1, ss2]:
            ss.set_state(1.0e5, 300.0)
            ss.set_composition(f0)

        self.assertArraysAlmostEqual(
            ss1.excess_partial_entropies, ss2.excess_partial_entropies
        )
        self.assertArraysAlmostEqual(
            ss1.entropy_hessian.flatten(), ss2.entropy_hessian.flatten()
        )

    def test_polynomial_using_subregular_volume(self):
        ss1 = two_site_ss_subregular_ternary()
        ss2 = two_site_ss_polynomial_ternary()
        f0 = [0.25, 0.35, 0.4]
        for ss in [ss1, ss2]:
            ss.set_state(1.0e5, 300.0)
            ss.set_composition(f0)

        self.assertArraysAlmostEqual(
            ss1.excess_partial_volumes, ss2.excess_partial_volumes
        )
        self.assertArraysAlmostEqual(
            ss1.volume_hessian.flatten(), ss2.volume_hessian.flatten()
        )

    def test_polynomial_mbr_mixing(self):
        ss = polynomial_mbr_mixing()
        f0 = [0.25, 0.75]
        ss.set_state(1.0e5, 300.0)
        ss.set_composition(f0)

        self.assertAlmostEqual(
            ss.excess_gibbs, 0.25 * 0.75 * (1200.0 - 300.0 * 5 + 1.0e5 * 1.0e-6)
        )
        self.assertAlmostEqual(ss.excess_entropy, 0.25 * 0.75 * 5)
        self.assertAlmostEqual(ss.excess_volume, 0.25 * 0.75 * 1.0e-6)

    def test_polynomial_ESV_and_mbr_mixing(self):
        ss = polynomial_esv_and_mbr_mixing()
        f0 = [0.25, 0.75]
        ss.set_state(1.0e5, 300.0)
        ss.set_composition(f0)
        self.assertAlmostEqual(
            ss.excess_gibbs, 0.25 * 0.75 * (1200.0 - 300.0 * 5 + 1.0e5 * 1.0e-6)
        )
        self.assertAlmostEqual(ss.excess_entropy, 0.25 * 0.75 * 5)
        self.assertAlmostEqual(ss.excess_volume, 0.25 * 0.75 * 1.0e-6)

    def test_polynomial_ternary_transformed(self):
        ss1 = two_site_ss_polynomial_ternary()
        ss2 = two_site_ss_polynomial_ternary_transformed()

        for ss in [ss1, ss2]:
            ss.set_state(1.0e5, 300.0)
        ss1.set_composition([0.25, 0.35, 0.4])
        ss2.set_composition([0.4, 0.25, 0.35])

        self.assertAlmostEqual(ss1.excess_gibbs, ss2.excess_gibbs)

    def test_subregular(self):
        ss0 = two_site_ss()
        ss1 = two_site_ss_subregular()

        ss0.set_composition([0.3, 0.3, 0.4])
        ss0.set_state(1.0e5, 300.0)

        ss1.set_composition([0.3, 0.3, 0.4])
        ss1.set_state(1.0e5, 300.0)

        self.assertArraysAlmostEqual(ss0.excess_partial_gibbs, ss1.excess_partial_gibbs)

    def test_activities_ideal(self):
        ol = olivine_ideal_ss()
        ol.set_composition(np.array([0.5, 0.5]))
        ol.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(ol.activities, [0.25, 0.25])

    def test_activity_coefficients_ideal(self):
        ol = olivine_ideal_ss()
        ol.set_composition(np.array([0.5, 0.5]))
        ol.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(ol.activity_coefficients, [1.0, 1.0])

    def test_activity_coefficients_non_ideal(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.0e5, 1000.0)
        self.assertArraysAlmostEqual(opx.activity_coefficients, [np.exp(1.0), 1.0])

    def test_formula(self):
        ol = olivine_ideal_ss()
        ol.set_composition([0.5, 0.5])
        self.assertEqual(formula_to_string(ol.formula), "MgFeSiO4")
        self.assertEqual(
            formula_to_string(sum_formulae(ol.endmember_formulae, [0.5, 0.5])),
            "MgFeSiO4",
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

    def test_temkin_activities(self):
        ss = temkin_ss()
        f0 = np.array([0.25, 0.35, 0.3, 0.1])
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)

        S_partial = ss.excess_partial_entropies
        activities = ss.activities

        self.assertArraysAlmostEqual(
            np.exp(-S_partial / burnman.constants.gas_constant), activities
        )

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

    def test_temkin_hessian(self):
        ss = temkin_ss()
        f0 = np.array([0.25, 0.35, 0.3, 0.1])
        ss.set_composition(f0)
        ss.set_state(1.0e5, 300.0)
        H0 = ss.gibbs_hessian

        df = np.array([2.0e-6, -1.5e-6, -0.6e-6, 0.1e-6])
        ss.set_composition(f0 - df / 2.0)
        dGdx1 = ss.partial_gibbs
        ss.set_composition(f0 + df / 2.0)
        dGdx2 = ss.partial_gibbs

        self.assertArraysAlmostEqual(H0.dot(df), dGdx2 - dGdx1)

    def test_function_solution(self):
        ss = [ppv_symmetric(), ppv_function()]
        for s in ss:
            s.set_state(1.0e5, 300.0)
            s.set_composition([0.2, 0.3, 0.5])

        self.assertArraysAlmostEqual(
            ss[0].excess_partial_gibbs, ss[1].excess_partial_gibbs
        )
        self.assertArraysAlmostEqual(ss[0].partial_gibbs, ss[1].partial_gibbs)
        self.assertArraysAlmostEqual(ss[0].partial_entropies, ss[1].partial_entropies)
        self.assertArraysAlmostEqual(ss[0].partial_volumes, ss[1].partial_volumes)
        self.assertArraysAlmostEqual(
            ss[0].gibbs_hessian.flatten(), ss[1].gibbs_hessian.flatten()
        )
        self.assertArraysAlmostEqual(
            ss[0].entropy_hessian.flatten(),
            ss[1].entropy_hessian.flatten(),
            tol_zero=1.0e-12,
        )
        self.assertArraysAlmostEqual(
            ss[0].volume_hessian.flatten(), ss[1].volume_hessian.flatten()
        )
        self.assertArraysAlmostEqual(ss[0].activities, ss[1].activities)
        self.assertArraysAlmostEqual(
            ss[0].activity_coefficients, ss[1].activity_coefficients
        )

    def test_function_solution_low_proportions(self):
        ss = [ppv_symmetric(), ppv_function()]
        for s in ss:
            s.set_state(1.0e5, 300.0)
            s.set_composition([0.0, 0.5, 0.5])

        self.assertArraysAlmostEqual(
            ss[0].excess_partial_gibbs, ss[1].excess_partial_gibbs
        )
        self.assertArraysAlmostEqual(ss[0].partial_gibbs, ss[1].partial_gibbs)
        self.assertArraysAlmostEqual(ss[0].partial_entropies, ss[1].partial_entropies)
        self.assertArraysAlmostEqual(ss[0].partial_volumes, ss[1].partial_volumes)
        self.assertArraysAlmostEqual(
            ss[0].gibbs_hessian.flatten(), ss[1].gibbs_hessian.flatten()
        )
        self.assertArraysAlmostEqual(
            ss[0].entropy_hessian.flatten(),
            ss[1].entropy_hessian.flatten(),
            tol_zero=1.0e-12,
        )
        self.assertArraysAlmostEqual(
            ss[0].volume_hessian.flatten(), ss[1].volume_hessian.flatten()
        )
        self.assertArraysAlmostEqual(ss[0].activities, ss[1].activities)
        self.assertArraysAlmostEqual(
            ss[0].activity_coefficients, ss[1].activity_coefficients
        )


if __name__ == "__main__":
    unittest.main()
