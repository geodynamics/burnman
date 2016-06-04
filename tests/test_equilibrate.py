from __future__ import absolute_import

import os
import sys
sys.path.insert(1, os.path.abspath('..'))
import warnings

import unittest
from util import BurnManTest

from burnman import Composite
from burnman import tools
from burnman import minerals
from burnman.equilibrate import *


fo = minerals.HP_2011_ds62.fo()
fa = minerals.HP_2011_ds62.fa()
sillimanite = minerals.HP_2011_ds62.sill()
andalusite = minerals.HP_2011_ds62.andalusite()
kyanite = minerals.HP_2011_ds62.ky()

assemblage = Composite([sillimanite, andalusite])
assemblage2 = Composite([minerals.SLB_2011.mg_fe_olivine(), minerals.SLB_2011.mg_fe_wadsleyite()])


class Minimization(BurnManTest):    
    def test_equilibrate(self):
        sol = equilibrate(sillimanite.params['formula'], assemblage, [['T', 800.], ['X', [1., 0.], [1., 1.], 0.]])
        self.assertArraysAlmostEqual((sol['c'] - np.array([0., 1.]))/1.e6, [0., 0.])
 
    def test_bulk_equilibrate(self):
        constraints = [['T', 1600.], ['P', 13.e9], ['X', [0., 0., 1., 0.], [1., 0., 1., 0.], 0.]]
        sol = find_equilibrium_composition(fo.params['formula'], fa.params['formula'], 0.71, assemblage2, constraints)
        self.assertFloatEqual(sol['c'][1], sol['X'])
    
    def test_invariant(self):
        PT0 = tools.invariant_point([sillimanite, andalusite], [1., -1.],
                                    [sillimanite, kyanite], [1., -1.],
                                    [0.4e9, 800.])

        PT1 = equilibrate_invariant(sillimanite.params['formula'], [andalusite], [sillimanite, kyanite])[0:2]
        self.assertArraysAlmostEqual(PT0, PT1)
    
    def test_univariant(self):
        P0 = tools.equilibrium_pressure([sillimanite, andalusite], [1., -1.], 800., 1.e9)
        P1 = equilibrate_univariant(sillimanite.params['formula'], [sillimanite], andalusite, 'T', [800.])[0][0]
        self.assertFloatEqual(P0, P1)

if __name__ == '__main__':
    unittest.main()
