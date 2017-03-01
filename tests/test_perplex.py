from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest

class PerpleX(BurnManTest):    
    def test_energies(self):
        P = 10.e9
        T = 3000.
        rock = burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')
    
        A, G, V = rock.evaluate(['molar_helmholtz', 'molar_gibbs', 'V'], [P], [T])
        self.assertFloatEqual(A, G - V*P)
        
    def test_grueneisen(self):
        P = 10.e9
        T = 3000.
        rock = burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')

        gr, alpha, V, K_S, cp = rock.evaluate(['grueneisen_parameter', 'alpha', 'V', 'adiabatic_bulk_modulus', 'heat_capacity_p'], [P], [T])
        self.assertFloatEqual(gr, alpha*V*K_S/cp)

        
if __name__ == '__main__':
    unittest.main()
