from __future__ import absolute_import
import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))
import warnings

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

class bcc_iron (Mineral):
    def __init__(self):
        formula='Fe1.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'BCC iron',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'P_0': 0.,
            'H_0': 9149.0 ,
            'S_0': 36.868 ,
            'V_0': 7.09e-06 ,
            'Cp': [21.09, 0.0101455, -221508., 47.1947] ,
            'a_0': 3.56e-05 ,
            'K_0': 1.64e+11 ,
            'Kprime_0': 5.16 ,
            'Kdprime_0': -3.1e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'curie_temperature': [1043., 0.0] ,
            'magnetic_moment': [2.22, 0.0] ,
            'magnetic_structural_parameter': 0.4 }
        Mineral.__init__(self)

class test_endmembers(BurnManTest):

    def test0_T_ref(self):
        fo=forsterite()
        fo.params['P_0'] = 0.99999999999e5

        fo.set_state(1.e5, 500.)
        data = [fo.H, fo.S]
        T_0 = 2000.
        fo.set_state(1.e5, T_0)
        fo.params['T_0'] = T_0
        fo.params['H_0'] = fo.H
        fo.params['S_0'] = fo.S
        fo.params['V_0'] = fo.V
        fo.params['K_0'] = fo.K_T 
        fo.params['Kdprime_0'] = -fo.params['Kprime_0']/fo.params['K_0']
        fo.params['a_0'] = fo.alpha

        fo.set_state(1.e5, 500.)
        data2 = [fo.H, fo.S]
        self.assertArraysAlmostEqual(data, data2)

    def test1_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e5, 298.15)
        volume1=fo.V
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e5, 298.15)
        volume2=fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test2_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e9, 298.15)
        volume1=fo.V
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e9, 298.15)
        volume2=fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test3_mt_hp_tmt(self):
        fo=forsterite()
        fo.set_state(1.e9, 298.15)
        K1=fo.K_T
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e9, 298.15)
        K2=fo.K_T
        self.assertArraysAlmostEqual([K1], [K2])

    def test4_bcc_iron_hp_tmt(self):
        bcc=bcc_iron()
        bcc.set_state(1.e5, 298.15)
        gibbs_mag=bcc.method._magnetic_gibbs(1.e5, 298.15, bcc.params)
        gibbs1=bcc.gibbs
        del bcc.params['curie_temperature']
        del bcc.params['magnetic_moment']
        del bcc.params['magnetic_structural_parameter']
        bcc.set_state(1.e5, 100.) # reset with different PT
        bcc.set_state(1.e5, 298.15)
        gibbs2=bcc.gibbs
        self.assertArraysAlmostEqual([gibbs1-gibbs_mag], [gibbs2])


if __name__ == '__main__':
    unittest.main()
