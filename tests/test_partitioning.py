from __future__ import absolute_import
import unittest
from util import BurnManTest

import burnman
from burnman import minerals



class test(BurnManTest):

    def test_simple(self):
        bulk_composition = burnman.Composition({'Mg': 0.213, 'Fe': 0.0626,
                                                'Si': 0.242, 'Ca': 0., 'Al': 0.0},
                                               unit_type='weight')

        bdg = minerals.SLB_2011.mg_fe_bridgmanite()
        per = minerals.SLB_2011.ferropericlase()

        pressure = 23.83e9  # Pa
        temperature = 2000.  # K
        (a, b) = burnman.calculate_nakajima_fp_pv_partition_coefficient(
            pressure, temperature, bulk_composition.molar_composition, 0.5)
        self.assertFloatEqual(a, 0.184533288)
        self.assertFloatEqual(b, 0.102937268)


if __name__ == '__main__':
    unittest.main()
