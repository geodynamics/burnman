import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.mineral import Mineral
from util import BurnManTest
from burnman.tools import *

class test_tools(BurnManTest):

    def test_eqm_T(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        H_ex = 1.e3
        S_ex = 2.
        fo2.params['H_0'] = fo.params['H_0'] + H_ex
        fo2.params['S_0'] = fo.params['S_0'] + S_ex

        
        T = H_ex/S_ex     
        T_calc = equilibrium_temperature([fo, fo2], [1.0, -1.0], fo.params['P_0'], 1200.)


        self.assertArraysAlmostEqual([T], [T_calc])
        
    def test_eqm_P(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        fo.params['K_0'] = 1.e20
        fo2.params['K_0'] = 1.e20

        H_ex = 1
        V_ex = 1.e-8
        fo2.params['H_0'] = fo.params['H_0'] + H_ex
        fo2.params['V_0'] = fo.params['V_0'] - V_ex

        
        P = fo.params['P_0'] + H_ex/V_ex
        P_calc = equilibrium_pressure([fo, fo2], [1.0, -1.0], fo.params['T_0'])
        self.assertArraysAlmostEqual([P], [P_calc])

        
    def test_fit_PVT_data(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        
        pressures = np.linspace(1.e9, 2.e9, 8)
        temperatures = np.ones_like(pressures)*fo.params['T_0']
        V = np.empty_like(pressures)
        
        PT = [pressures, temperatures]
        for i, P in enumerate(pressures):
            fo.set_state(P, temperatures[i])
            V[i] = fo.V
        
        params = ['V_0', 'K_0', 'Kprime_0']
        popt, pcov = fit_PVT_data(fo, params, PT, V)
        zeros = np.zeros_like(pcov[0])
        
        self.assertArraysAlmostEqual(pcov[0], zeros)


    def test_hugoniot(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        T_ref = 298.15
        P_ref = 1.e5
        pressures = np.array([P_ref])
        temperatures, volumes = hugoniot(fo, P_ref, T_ref, pressures)

        self.assertArraysAlmostEqual(temperatures, np.array([T_ref]))

    def test_fraction_conversion_mass(self):
        per = burnman.minerals.SLB_2011.periclase()
        stv = burnman.minerals.SLB_2011.stishovite()
        fo = burnman.minerals.SLB_2011.forsterite()
        c = burnman.Composite([per, stv, fo], [0.5, 0.25, 0.25])
        mass_fractions = convert_fractions(c, [0.5, 0.25, 0.25], 'molar', 'mass')
        self.assertArraysAlmostEqual([mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]])
        
    def test_fraction_conversion_volume(self):
        per = burnman.minerals.SLB_2011.periclase()
        c = burnman.Composite([per, per, per], [0.25, 0.25, 0.5])
        c.set_state(1.e5, 300.)
        mass_fractions = convert_fractions(c, [0.25, 0.25, 0.5], 'molar', 'volume')
        self.assertArraysAlmostEqual([mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]])
        
if __name__ == '__main__':
    unittest.main()
