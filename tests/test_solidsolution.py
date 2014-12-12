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
           'equation_of_state': 'mtait',
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

# Configurational entropy
class forsterite_ss(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='Dummy solid solution'

        base_material = [[forsterite(), '[Mg]2SiO4']]

        # Interaction parameters
        enthalpy_interaction=[]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )




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
        
    def test_1_gibbs(self):
        fo, fo_ss = self.setup_1min_ss()
        endmember_properties=[fo.gibbs, fo.H, fo.S, fo.V, fo.C_p, fo.C_v, fo.alpha, fo.K_T, fo.K_S, fo.gr]
        ss_properties=[fo_ss.gibbs, fo_ss.H, fo_ss.S, fo_ss.V, fo_ss.C_p, fo_ss.C_v, fo_ss.alpha, fo_ss.K_T, fo_ss.K_S, fo_ss.gr]
        self.assertArraysAlmostEqual(endmember_properties, ss_properties)

if __name__ == '__main__':
    unittest.main()
