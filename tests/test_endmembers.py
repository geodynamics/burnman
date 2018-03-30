from __future__ import absolute_import
import unittest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman.mineral import Mineral
from burnman.processchemistry import dictionarize_formula, formula_mass
from burnman.combinedmineral import CombinedMineral
from util import BurnManTest
import copy



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


class test_endmembers(BurnManTest):

    def test0_T_ref(self):
        fo = forsterite()
        fo.params['P_0'] = 0.99999999999e5
        fo.set_state(1.e5, 500.)
        data = [fo.H, fo.S]
        T_0 = 2000.
        fo.set_state(1.e5, T_0)
        fo2 = copy.deepcopy(fo)
        fo2.params['T_0'] = T_0
        fo2.params['H_0'] = fo.H
        fo2.params['S_0'] = fo.S
        fo2.params['V_0'] = fo.V
        fo2.params['K_0'] = fo.K_T
        fo2.params['Kdprime_0'] = -fo.params['Kprime_0'] / fo2.params['K_0']
        fo2.params['a_0'] = fo.alpha
        fo2.set_state(1.e5, 500.)

        data2 = [fo2.H, fo2.S]

        self.assertArraysAlmostEqual(data, data2)

    def test1_mt_hp_tmt(self):
        fo = forsterite()
        fo.set_state(1.e5, 298.15)
        volume1 = fo.V
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e5, 298.15)
        volume2 = fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test2_mt_hp_tmt(self):
        fo = forsterite()
        fo.set_state(1.e9, 298.15)
        volume1 = fo.V
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e9, 298.15)
        volume2 = fo.V
        self.assertArraysAlmostEqual([volume1], [volume2])

    def test3_mt_hp_tmt(self):
        fo = forsterite()
        fo.set_state(1.e9, 298.15)
        K1 = fo.K_T
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            fo.set_method('mt')
            assert len(w) == 1
        fo.set_state(1.e9, 298.15)
        K2 = fo.K_T
        self.assertArraysAlmostEqual([K1], [K2])

    def test_mgd(self):
        bdg = burnman.minerals.Matas_etal_2007.mg_perovskite()
        bdg.set_state(1.e5, 300.)
        self.assertFloatEqual(
            bdg.grueneisen_parameter, bdg.params['grueneisen_0'])
        self.assertFloatEqual(bdg.isothermal_bulk_modulus, bdg.params['K_0'])
        self.assertFloatEqual(bdg.molar_volume, bdg.params['V_0'])
        self.assertFloatEqual(bdg.shear_modulus, bdg.params['G_0'])
        bdg.set_state(1.e5, 0.)
        self.assertFloatEqual(bdg.molar_heat_capacity_v, 0.)
        self.assertFloatEqual(bdg.molar_heat_capacity_p, 0.)
        self.assertFloatEqual(bdg.thermal_expansivity, 0.)
        self.assertFloatEqual(
            bdg.isothermal_bulk_modulus, bdg.adiabatic_bulk_modulus)

    def test_make_mbr(self):
     bdg = burnman.minerals.SLB_2011.mg_perovskite()
     dS = 1.
     made_bdg = CombinedMineral([bdg, bdg], [2., -1.], [0., dS, 0.])
     bdg.set_state(1.e5, 1000.)
     made_bdg.set_state(1.e5, 1000.)
     self.assertFloatEqual(bdg.S + dS, made_bdg.S)
     
    def test_make_mbr2(self):
     per = burnman.minerals.SLB_2011.periclase()
     stv = burnman.minerals.SLB_2011.stishovite()
     dE = -15000.
     made_bdg = CombinedMineral([burnman.minerals.SLB_2011.periclase(),
                                 burnman.minerals.SLB_2011.stishovite()],
                                [1., 1.], [dE, 0., 0.])
     
     per.set_state(1.e5, 1000.)
     stv.set_state(1.e5, 1000.)
     made_bdg.set_state(1.e5, 1000.)
     self.assertFloatEqual(made_bdg.V, per.V + stv.V)

     
        
if __name__ == '__main__':
    unittest.main()
