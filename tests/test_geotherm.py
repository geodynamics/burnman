import unittest
import os, sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


class mypericlase(burnman.Mineral):
    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein 
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8}


class geotherm(BurnManTest):
    def test_adiabat(self):
        rock = mypericlase()
        rock.params['K_0'] = 0.
        pressure = [10e9]
        T0 = 300.
        test_K_adiabat = burnman.geotherm.adiabatic(pressure,T0,rock.params)
        self.assertFloatEqual(test_K_adiabat,T0)




if __name__ == '__main__':
    unittest.main()
