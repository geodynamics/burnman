from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest

class EosConsistency(BurnManTest):

    def test_HP(self):
        P = 10.e9
        T = 3000.
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            self.assertEqual(burnman.tools.check_eos_consistency(burnman.minerals.HP_2011_ds62.per(), P, T),
                             True)

    def test_SLB(self):
        P = 10.e9
        T = 3000.
        self.assertEqual(burnman.tools.check_eos_consistency(burnman.minerals.SLB_2011.periclase(), P, T),
                         True)

    def test_solution(self):
        P = 10.e9
        T = 3000.
        m = burnman.minerals.SLB_2011.garnet(molar_fractions = [0.2, 0.2, 0.2, 0.2, 0.2])
        self.assertEqual(burnman.tools.check_eos_consistency(m, P, T),
                         True)

if __name__ == '__main__':
    unittest.main()
