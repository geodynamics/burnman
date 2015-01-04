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


class test_endmembers(BurnManTest):

    def test1_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e5, 298.15)
        volume1=fo.V
        fo.set_method('mt')
        fo.set_state(1.e5, 298.15)
        volume2=fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test2_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e9, 298.15)
        volume1=fo.V
        fo.set_method('mt')
        fo.set_state(1.e9, 298.15)
        volume2=fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test3_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e9, 298.15)
        K1=fo.K_T
        fo.set_method('mt')
        fo.set_state(1.e9, 298.15)
        K2=fo.K_T
        self.assertArraysAlmostEqual([K1], [K2])

if __name__ == '__main__':
    unittest.main()
