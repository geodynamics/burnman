from __future__ import absolute_import
import unittest
from util import BurnManTest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals


# TODO: test composite that changes number of entries

class test(BurnManTest):

    def test_simple(self):
        inp1 = {'Mg': 0.213, 'Fe': 0.0626,
                'Si': 0.242, 'Ca': 0., 'Al': 0.}  # wt%
        phase_per, rel_mol_per = burnman.calculate_phase_percents(inp1)

        StartP = 23.83  # in GPa
        EndP = 110.0
        deltaP = 1.

        # P,T,a,b,frac_mol_pv,frac_mol_mw    =
        # part_coef_calc(inp2,StartP,EndP,deltaP)

        gt = lambda p: burnman.geotherm.anderson(p)
        pressure = StartP *1.e9
        temperature = 2000
        (a, b) = burnman.calculate_partition_coefficient(
            pressure, temperature, rel_mol_per, 0.5)
        self.assertFloatEqual(a, 0.18453519778)
        self.assertFloatEqual(b, 0.102938439776)


if __name__ == '__main__':
    unittest.main()
