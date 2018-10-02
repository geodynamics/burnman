from __future__ import absolute_import
import unittest
from util import BurnManTest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals
from burnman.processchemistry import convert_formula

class test(BurnManTest):

    def test_simple(self):
        bulk_composition_wt = {'Mg': 0.213, 'Fe': 0.0626,
                               'Si': 0.242, 'Ca': 0., 'Al': 0.0}
        bulk_composition_mol = convert_formula(bulk_composition_wt,
                                               to_type='molar')

        bdg = minerals.SLB_2011.mg_fe_bridgmanite()
        per = minerals.SLB_2011.ferropericlase()
        
        assemblage = burnman.Composite([bdg, per])
        assemblage.set_potential_composition_from_bulk(bulk_composition_mol,
                                                       unfitted_elements='O',
                                                       use_solution_guesses=False)

        f_per = assemblage.molar_fractions[1]

        self.assertFloatEqual(f_per, 0.12828483) # value governed only by (Mg+Fe)/(Mg+Fe+Si)

        pressure = 23.83e9 # Pa
        temperature = 2000. # K
        (a, b) = burnman.calculate_nakajima_fp_pv_partition_coefficient(
            pressure, temperature, bulk_composition_mol, 0.5)
        self.assertFloatEqual(a, 0.184533288)
        self.assertFloatEqual(b, 0.102937268)

    def test_guesses(self):
        bulk_composition_wt = {'Mg': 0.213, 'Fe': 0.0626,
                               'Si': 0.242, 'Ca': 0., 'Al': 0.0}
        bulk_composition_mol = convert_formula(bulk_composition_wt,
                                               to_type='molar')

        bdg = minerals.SLB_2011.mg_fe_bridgmanite()
        per = minerals.SLB_2011.ferropericlase()

        bdg.guess = [0.92, 0.08, 0.0]
        
        assemblage = burnman.Composite([bdg, per])
        assemblage.set_potential_composition_from_bulk(bulk_composition_mol,
                                                       unfitted_elements='O',
                                                       use_solution_guesses=True)

        f_per = assemblage.molar_fractions[1]
        self.assertFloatEqual(f_per, 0.12828483) # value governed only by (Mg+Fe)/(Mg+Fe+Si)
        self.assertArraysAlmostEqual(bdg.guess, assemblage.phases[0].molar_fractions)

        pressure = 23.83e9 # Pa
        temperature = 2000. # K
        (a, b) = burnman.calculate_nakajima_fp_pv_partition_coefficient(
            pressure, temperature, bulk_composition_mol, 0.5)
        self.assertFloatEqual(a, 0.184533288)
        self.assertFloatEqual(b, 0.102937268)


    def test_negative_endmembers(self):

        composition = [-0.1, 0.6,  -0.2, 0.21, 0.49]
        g = minerals.JH_2015.garnet()
        g.set_composition(composition)

        g2 = minerals.JH_2015.garnet()
        assemblage = burnman.Composite([g2])
        
        assemblage.set_potential_composition_from_bulk(g.formula,
                                                       unfitted_elements='O',
                                                       use_solution_guesses=False)
        self.assertArraysAlmostEqual(composition, assemblage.phases[0].molar_fractions)
        
if __name__ == '__main__':
    unittest.main()
