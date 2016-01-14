from __future__ import absolute_import
import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.mineral import Mineral
from burnman.processchemistry import *
from util import BurnManTest

atomic_masses=read_masses()

class forsterite (Mineral):
    def __init__(self):
       formula='Mg2.0Si1.0O4.0'
       formula = dictionarize_formula(formula)
       self.params = {
           'name': 'fo',
           'formula': formula,
           'equation_of_state': 'hp_tmt',
           'H_0': -2172590.0 ,
           'S_0': 95.1 ,
           'V_0': 4.366e-05 ,
           'Cp': [233.3, 0.001494, -603800.0, -1869.7] ,
           'a_0': 2.85e-05 ,
           'K_0': 1.285e+11 ,
           'Kprime_0': 3.84 ,
           'Kdprime_0': -3e-11 ,
           'n': sum(formula.values()),
           'molar_mass': formula_mass(formula, atomic_masses)}
       Mineral.__init__(self)

class fayalite (Mineral):
    def __init__(self):
       formula='Fe2.0Si1.0O4.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'fa',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': -1477720.0 ,
            'S_0': 151.0 ,
            'V_0': 4.631e-05 ,
            'Cp': [201.1, 0.01733, -1960600.0, -900.9] ,
            'a_0': 2.82e-05 ,
            'K_0': 1.256e+11 ,
            'Kprime_0': 4.68 ,
            'Kdprime_0': -3.7e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}
       Mineral.__init__(self)

# One-mineral solid solution
class forsterite_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name='Dummy solid solution'
        self.type='symmetric'
        self.endmembers = [[forsterite(), '[Mg]2SiO4']]
        self.enthalpy_interaction=[]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Two-mineral solid solution
class forsterite_forsterite_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name='Fo-Fo solid solution'
        self.type='symmetric'
        self.endmembers = [[forsterite(), '[Mg]2SiO4'], [forsterite(), '[Mg]2SiO4']]
        self.enthalpy_interaction=[[0.]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Olivine solid solution
class olivine_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name='Olivine'
        self.type='symmetric'
        self.endmembers = [[forsterite(), '[Mg]2SiO4'], [fayalite(), '[Fe]2SiO4']]
        self.enthalpy_interaction=[[8.4e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Orthopyroxene solid solution
class orthopyroxene(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        # Name
        self.name='orthopyroxene'
        self.type='symmetric'
        self.endmembers = [[forsterite(), 'Mg[Mg]Si2O6'], [forsterite(), '[Mg1/2Al1/2]2AlSiO6']]
        self.enthalpy_interaction=[[10.0e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Three-endmember, two site solid solution
class two_site_ss(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        self.name='two_site_ss'
        self.type='symmetric'
        self.endmembers = [[forsterite(), '[Mg]3[Al]2Si3O12'],[forsterite(), '[Fe]3[Al]2Si3O12'],[forsterite(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
        self.enthalpy_interaction=[[10.0e3, 5.0e3],[-10.0e3]]

        burnman.SolidSolution.__init__(self, molar_fractions)

# Three-endmember, two site solid solution
class two_site_ss_subregular(burnman.SolidSolution):
    def __init__(self, molar_fractions=None):
        # Name
        self.name='two_site_ss (subregular symmetric)'
        self.type='subregular'
        self.endmembers = [[forsterite(), '[Mg]3[Al]2Si3O12'],[forsterite(), '[Fe]3[Al]2Si3O12'],[forsterite(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
        # Interaction parameters
        self.enthalpy_interaction=[[[10.e3, 10.e3],[5.e3, 5.e3]],[[-10.e3, -10.e3]]]

        burnman.SolidSolution.__init__(self, molar_fractions )


class test_solidsolution(BurnManTest):

    def setup_1min_ss(self):
        P=1.e5
        T=1000.
        fo=forsterite()
        fo.set_state(P,T)
        
        fo_ss=forsterite_ss()
        fo_ss.set_composition([1.0])
        fo_ss.set_state(P,T)
        return fo, fo_ss
        

    def setup_2min_ss(self):
        P=1.e5
        T=1000.
        fo=forsterite()
        fo.set_state(P,T)
        
        fo_fo_ss=forsterite_forsterite_ss()
        fo_fo_ss.set_composition([0.3, 0.7])
        fo_fo_ss.set_state(P,T)
        return fo, fo_fo_ss

    def setup_ol_ss(self):
        P=1.e5
        T=1000.
        fo=forsterite()
        fo.set_state(P,T)
        
        ol_ss=olivine_ss()
        ol_ss.set_composition([1.0, 0.0])
        ol_ss.set_state(P,T)
        return fo, ol_ss
        

    def test_1_gibbs(self):
        fo, fo_ss = self.setup_1min_ss()
        endmember_properties=[fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr, fo.G]
        ss_properties=[fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr, fo_ss.G]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_2_gibbs(self):
        fo, fo_ss = self.setup_2min_ss()
        endmember_properties=[fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr, fo.G]
        ss_properties=[fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr, fo_ss.G]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_gibbs(self):
        fo, fo_ss = self.setup_ol_ss()
        endmember_properties=[fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr]
        ss_properties=[fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

    def test_ol_Wh(self):
        ol_ss=olivine_ss()
        H_excess=ol_ss.solution_model.excess_enthalpy(1.e5, 1000., [0.5,0.5])
        Wh=ol_ss.solution_model.Wh[0][1]
        self.assertArraysAlmostEqual([Wh/4.0], [H_excess])

    def test_order_disorder(self):
        opx = orthopyroxene()
        opx.set_composition( np.array([0.0, 1.0]) )
        opx.set_state(1.e5,300.)
        self.assertArraysAlmostEqual([opx.excess_gibbs], [0.])

    def test_site_totals(self):
        ss=two_site_ss()
        ss.set_composition([0.3,0.3,0.4])
        ss.set_state(1.e5,300.)

        site_fractions=np.dot(ss.molar_fractions, ss.solution_model.endmember_occupancies)
        i=0
        site_fill=[]
        ones=[1.] * ss.solution_model.n_sites
        for site in ss.solution_model.sites:
            site_fill.append(sum(site_fractions[i:i+len(site)]))
            i += len(site)

        self.assertArraysAlmostEqual(site_fill, ones)

    def test_set_method(self):
        ss = olivine_ss()
        ss.set_method('hp_tmt')

    def test_molar_mass(self):
        ss = olivine_ss()
        ss.set_composition( np.array([0.5, 0.5]) )
        self.assertArraysAlmostEqual([ss.molar_mass], [0.5*forsterite().params['molar_mass']+0.5*fayalite().params['molar_mass']])

    def test_subregular(self):
        ss0=two_site_ss()
        ss1=two_site_ss_subregular()

        ss0.set_composition([0.3,0.3,0.4])
        ss0.set_state(1.e5,300.)

        ss1.set_composition([0.3,0.3,0.4])
        ss1.set_state(1.e5,300.)

        self.assertArraysAlmostEqual(ss0.excess_partial_gibbs, ss1.excess_partial_gibbs)


if __name__ == '__main__':
    unittest.main()
