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


class Debye(BurnManTest):
    def test_temperature(self):
        #rock = mypericlase()
        pressure = 0.
        temperature = 300.
        eoses = [burnman.slb.SLB2(), burnman.slb.SLB3(), burnman.birch_murnaghan.BM2(), burnman.birch_murnaghan.BM3()]

        test_debye = burnman.debye.debye_fn(temperature)
        test_debye_cheb = burnman.debye.debye_fn_cheb(temperature)
        self.assertFloatEqual(test_debye, test_debye_cheb)


if __name__ == '__main__':
    unittest.main()
