from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


m = minerals.HP_2011_ds62.per()

class EosConsistency(BurnManTest):

    def test_HP(self):
        P = 10.e9
        T = 3000.
        dT = 1.
        dP = 100.
        
        m.set_state(P, T)
        G0 = m.gibbs
        S0 = m.S
        V0 = m.V
        
        l0 = [m.gibbs/(m.helmholtz + P*m.V), m.gibbs/(m.H - T*m.S), m.gibbs/(m.internal_energy - T*m.S + P*m.V)]
        l1 = [1., 1., 1.]

        m.set_state(P, T + dT)
        G1 = m.gibbs
        S1 = m.S

        
        m.set_state(P + dP, T)
        G2 = m.gibbs
        V2 = m.V
        
        # T derivatives
        m.set_state(P, T + 0.5*dT)
        l0.extend([(-(G1 - G0)/dT/m.S - 1.)/1000. + 1.,
                   ((T + 0.5*dT)*(S1 - S0)/dT/m.heat_capacity_p - 1.)/1000. + 1.])
        l1.extend([1., 1.])

        # P derivatives
        m.set_state(P + 0.5*dP, T)
        l0.extend([((G2 - G0)/dP/m.V - 1.)/1000. + 1., (-0.5*(V0 + V2)*dP/(V2 - V0)/m.K_T - 1.)/1000. + 1.])
        l1.extend([1., 1.])

        self.assertArraysAlmostEqual(l0, l1)


if __name__ == '__main__':
    unittest.main()
