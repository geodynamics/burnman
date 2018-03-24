from __future__ import absolute_import
import unittest
from util import BurnManTest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals
from burnman.processchemistry import convert_formula, calculate_potential_phase_amounts

class test(BurnManTest):

    def test_simple(self):
        bulk_composition_wt = {'Mg': 0.213, 'Fe': 0.0626,
                               'Si': 0.242, 'Ca': 0., 'Al': 0.}
        bulk_composition_mol = convert_formula(bulk_composition_wt,
                                               to_type='molar')

        per = minerals.SLB_2011.ferropericlase()
        bdg = minerals.SLB_2011.mg_fe_bridgmanite()
        formulae = per.endmember_formulae
        formulae.extend(bdg.endmember_formulae)

        phase_amounts = calculate_potential_phase_amounts(bulk_composition_mol, formulae)
        f_per = sum(phase_amounts[0:2])/sum(phase_amounts)
        self.assertFloatEqual(f_per, 0.12828483)

        pressure = 23.83e9 # Pa
        temperature = 2000. # K
        (a, b) = burnman.calculate_nakajima_fp_pv_partition_coefficient(
            pressure, temperature, bulk_composition_mol, 0.5)
        self.assertFloatEqual(a, 0.184533288)
        self.assertFloatEqual(b, 0.102937268)

if __name__ == '__main__':
    unittest.main()
