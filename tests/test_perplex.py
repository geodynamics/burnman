from __future__ import absolute_import
import unittest
import os
import sys
import numpy as np

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest

# Predefine our rock
rock = burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')

class PerpleX(BurnManTest):
        
    def test_energies(self):
        P = [10.e9, 11.e9]
        A, G, V = rock.evaluate(['molar_helmholtz', 'molar_gibbs', 'V'],
                                P, [2000., 2000.])
        self.assertArraysAlmostEqual(A, G - V*P)
        
    def test_grueneisen(self):

        gr, alpha, V, K_S, cp = rock.evaluate(['grueneisen_parameter', 'alpha', 'V', 'adiabatic_bulk_modulus', 'molar_heat_capacity_p'],
                                              [10.e9, 11.e9], [2000., 2000.])
        self.assertArraysAlmostEqual(gr, alpha*V*K_S/cp)

    def test_bounds(self):
        def fnLP():
            return rock.set_state(rock.bounds[0][0] - 1., np.mean(rock.bounds[1]))
        def fnHP():
            return rock.set_state(rock.bounds[0][1] + 1., np.mean(rock.bounds[1]))
        def fnLT():
            return rock.set_state(np.mean(rock.bounds[0]), rock.bounds[1][0] - 1.)
        def fnHT():
            return rock.set_state(np.mean(rock.bounds[0]), rock.bounds[1][1] + 1.)
        
        with np.errstate(all='ignore'):
            self.assertRaises(Exception, fnLP)
            self.assertRaises(Exception, fnHP)
            self.assertRaises(Exception, fnLT)
            self.assertRaises(Exception, fnHT)

        
if __name__ == '__main__':
    unittest.main()
