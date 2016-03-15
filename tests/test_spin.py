from __future__ import absolute_import
import unittest
from util import BurnManTest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals


class spin_transition(BurnManTest):

    def test_new(self):

        mins = [
            minerals.Murakami_etal_2012.fe_periclase(), minerals.Murakami_etal_2012.fe_periclase_HS(), minerals.Murakami_etal_2012.fe_periclase_LS()]
        for p in mins:
            p.set_method('slb2')

        # print "HS regime: (on/high/low)"
        for p in mins:
            p.set_state(5e9, 300)

        self.assertFloatEqual(mins[0].v_s, mins[1].v_s)

        # print "LS regime: (on/high/low)"
        for p in mins:
            p.set_state(70e9, 300)

        self.assertFloatEqual(mins[0].v_s, mins[2].v_s)

    def test_no_set_state(self):
        m = minerals.Murakami_etal_2012.fe_periclase()
        m.set_state(5e9, 300)
        self.assertArraysAlmostEqual(m.molar_fractions, [0.0, 1.0])
        m.set_state(70e9, 300)
        self.assertArraysAlmostEqual(m.molar_fractions, [1.0, 0.0])

if __name__ == '__main__':
    unittest.main()
