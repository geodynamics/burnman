from __future__ import absolute_import
import unittest
import os
import sys
import warnings
sys.path.insert(1, os.path.abspath('..'))

import numpy as np

import burnman
from burnman.mineral import Mineral
from burnman.combinedmineral import CombinedMineral
from burnman.processchemistry import dictionarize_formula, formula_mass
from util import BurnManTest

class forsterite (Mineral):

    def __init__(self):
        formula = 'Mg2.0Si1.0O4.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'fo',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': -2172590.0,
            'S_0': 95.1,
            'V_0': 4.366e-05,
            'Cp': [233.3, 0.001494, -603800.0, -1869.7],
            'a_0': 2.85e-05,
            'K_0': 1.285e+11,
            'Kprime_0': 3.84,
            'Kdprime_0': -3e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)


class fayalite (Mineral):

    def __init__(self):
        formula = 'Fe2.0Si1.0O4.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'fa',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': -1477720.0,
            'S_0': 151.0,
            'V_0': 4.631e-05,
            'Cp': [201.1, 0.01733, -1960600.0, -900.9],
            'a_0': 2.82e-05,
            'K_0': 1.256e+11,
            'Kprime_0': 4.68,
            'Kdprime_0': -3.7e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)}
        Mineral.__init__(self)


made_forsterite = CombinedMineral([forsterite(), forsterite()], [0.5, 0.5])

# One-mineral solid solution
class forsterite_ss(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'Dummy solid solution'
        self.solution_type = 'symmetric'
        self.endmembers = [[forsterite(), '[Mg]2SiO4']]
        self.energy_interaction = []

        burnman.SolidSolution.__init__(self, molar_fractions)

# Two-mineral solid solution
class forsterite_forsterite_ss(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'Fo-Fo solid solution'
        self.solution_type = 'symmetric'
        self.endmembers = [[forsterite(), '[Mg]2SiO4'], [
                           forsterite(), '[Mg]2SiO4']]
        self.energy_interaction = [[0.]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Ideal solid solution
class olivine_ideal_ss(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'Fo-Fo solid solution'
        self.solution_type = 'ideal'
        self.endmembers = [[
            forsterite(), '[Mg]2SiO4'], [fayalite(), '[Fe]2SiO4']]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Olivine solid solution
class olivine_ss(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'Olivine'
        self.solution_type = 'symmetric'
        self.endmembers = [[
            forsterite(), '[Mg]2SiO4'], [fayalite(), '[Fe]2SiO4']]
        self.energy_interaction = [[8.4e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

        
# Olivine solid solution with combined endmember
class olivine_ss2(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'Olivine'
        self.solution_type = 'symmetric'
        self.endmembers = [[
            made_forsterite, '[Mg]2SiO4'], [fayalite(), '[Fe]2SiO4']]
        self.energy_interaction = [[8.4e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Orthopyroxene solid solution
class orthopyroxene(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        # Name
        self.name = 'orthopyroxene'
        self.solution_type = 'symmetric'
        self.endmembers = [[forsterite(), '[Mg][Mg]Si2O6'], [
                           forsterite(), '[Mg1/2Al1/2][Mg1/2Al1/2]AlSiO6']]
        self.energy_interaction = [[burnman.constants.gas_constant * 1.0e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Three-endmember, two site symmetric solid solution
class two_site_ss(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'two_site_ss'
        self.solution_type = 'symmetric'
        self.endmembers = [[forsterite(), '[Mg]3[Al]2Si3O12'], [
                           forsterite(), '[Fe]3[Al]2Si3O12'], [forsterite(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
        self.energy_interaction = [[10.0e3, 5.0e3], [-10.0e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)
        
# Three-endmember, two site asymmetric solid solution
class two_site_ss_asymmetric(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        self.name = 'two_site_ss (asymmetric)'
        self.solution_type = 'asymmetric'
        self.endmembers = [[forsterite(), '[Mg]3[Al]2Si3O12'], [
                           forsterite(), '[Fe]3[Al]2Si3O12'], [forsterite(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
        self.alphas = [1., 2., 2.]
        self.energy_interaction = [[10.0e3, 5.0e3], [-10.0e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Three-endmember, two site solid solution
class two_site_ss_subregular(burnman.SolidSolution):

    def __init__(self, molar_fractions=None):
        # Name
        self.name = 'two_site_ss (subregular symmetric)'
        self.solution_type = 'subregular'
        self.endmembers = [[forsterite(), '[Mg]3[Al]2Si3O12'], [
                           forsterite(), '[Fe]3[Al]2Si3O12'], [forsterite(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
        # Interaction parameters
        self.energy_interaction = [
            [[10.e3, 10.e3], [5.e3, 5.e3]], [[-10.e3, -10.e3]]]

        burnman.SolidSolution.__init__(self, molar_fractions)


class test_solidsolution(BurnManTest):

    def setup_1min_ss(self):
        P = 1.e5
        T = 1000.
        fo = forsterite()
        fo.set_state(P, T)
        fo_ss = forsterite_ss()
        fo_ss.set_composition([1.0])
        fo_ss.set_state(P, T)
        return fo, fo_ss

    def setup_2min_ss(self):
        P = 1.e5
        T = 1000.
        fo = forsterite()
        fo.set_state(P, T)
        fo_fo_ss = forsterite_forsterite_ss()
        fo_fo_ss.set_composition([0.3, 0.7])
        fo_fo_ss.set_state(P, T)
        return fo, fo_fo_ss

    def setup_ol_ss(self):
        P = 1.e5
        T = 1000.
        fo = forsterite()
        fo.set_state(P, T)

        ol_ss = olivine_ss()
        ol_ss.set_composition([1.0, 0.0])
        ol_ss.set_state(P, T)
        return fo, ol_ss

    def setup_ol_ss2(self):
        P = 1.e5
        T = 1000.
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
            endmember_properties = [fo.gibbs, fo.H, fo.S, fo.V, fo.C_p,
                                    fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr, fo.G]
            ss_properties = [fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p,
                             fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr, fo_ss.G]
            assert len(w) == 3  # we expect to trigger 3 shear modulus warnings
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_2_gibbs(self):
        fo, fo_ss = self.setup_2min_ss()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            endmember_properties = [fo.gibbs, fo.H, fo.S, fo.V, fo.C_p,
                                    fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr, fo.G]
            ss_properties = [fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p,
                             fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr, fo_ss.G]
            assert len(w) == 4  # we expect to trigger 4 shear modulus warnings
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_gibbs(self):
        fo, fo_ss = self.setup_ol_ss()
        endmember_properties = [
            fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr]
        ss_properties = [fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V,
                         fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)
        
    def test_ol_gibbs2(self):
        fo, fo_ss = self.setup_ol_ss2()
        endmember_properties = [
            fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr]
        ss_properties = [fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V,
                         fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_Wh(self):
        ol_ss = olivine_ss()
        H_excess = ol_ss.solution_model.excess_enthalpy(
            1.e5, 1000., [0.5, 0.5]) # Hxs = Exs if Vxs=0
        We = ol_ss.solution_model.We[0][1]
        self.assertArraysAlmostEqual([We / 4.0], [H_excess])

    def test_order_disorder(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.e5, 300.)
        self.assertArraysAlmostEqual([opx.excess_gibbs], [0.])

    def test_site_totals(self):
        ss = two_site_ss()
        ss.set_composition([0.3, 0.3, 0.4])
        ss.set_state(1.e5, 300.)

        site_fractions = np.dot(
            ss.molar_fractions, ss.solution_model.endmember_occupancies)
        i = 0
        site_fill = []
        ones = [1.] * ss.solution_model.n_sites
        for site in ss.solution_model.sites:
            site_fill.append(sum(site_fractions[i:i + len(site)]))
            i += len(site)

        self.assertArraysAlmostEqual(site_fill, ones)

    def test_set_method(self):
        ss = olivine_ss()
        ss.set_method('hp_tmt')

    def test_molar_mass(self):
        ss = olivine_ss()
        ss.set_composition(np.array([0.5, 0.5]))
        self.assertArraysAlmostEqual([ss.molar_mass], [0.5 *
                                     forsterite().params['molar_mass'] + 0.5 * fayalite().params['molar_mass']])

    def test_hessian(self):
        ss = two_site_ss_asymmetric()
        f0 = [0.25, 0.35, 0.4]
        ss.set_composition([0.25, 0.35, 0.4])
        ss.set_state(1.e5, 300.)
        E0 = ss.excess_partial_gibbs
        dp = 1.e-8
        H = np.zeros((3, 3))
        for i in range(3):
            f = np.array(f0)
            f[i] += dp
            ss.set_composition(f)
            H[i,:] = (ss.excess_partial_gibbs - E0)/dp
        
        for i in range(3):
            self.assertArraysAlmostEqual(ss.gibbs_hessian[i], H[i])

    def test_subregular(self):
        ss0 = two_site_ss()
        ss1 = two_site_ss_subregular()

        ss0.set_composition([0.3, 0.3, 0.4])
        ss0.set_state(1.e5, 300.)

        ss1.set_composition([0.3, 0.3, 0.4])
        ss1.set_state(1.e5, 300.)

        self.assertArraysAlmostEqual(
            ss0.excess_partial_gibbs, ss1.excess_partial_gibbs)

    def test_activities_ideal(self):
        ol = olivine_ideal_ss()
        ol.set_composition(np.array([0.5, 0.5]))
        ol.set_state(1.e5, 1000.)
        self.assertArraysAlmostEqual(ol.activities, [0.25, 0.25])

    def test_activity_coefficients_ideal(self):
        ol = olivine_ideal_ss()
        ol.set_composition(np.array([0.5, 0.5]))
        ol.set_state(1.e5, 1000.)
        self.assertArraysAlmostEqual(ol.activity_coefficients, [1., 1.])

    def test_activity_coefficients_non_ideal(self):
        opx = orthopyroxene()
        opx.set_composition(np.array([0.0, 1.0]))
        opx.set_state(1.e5, 1000.)
        self.assertArraysAlmostEqual(
            opx.activity_coefficients, [np.exp(1.), 1.])

if __name__ == '__main__':
    unittest.main()
