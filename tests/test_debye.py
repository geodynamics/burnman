from __future__ import absolute_import
import unittest
import os
import sys

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

    def test_same_debye(self):
        x = 300.
        test_debye = burnman.eos.debye.debye_fn(x)
        test_debye_cheb = burnman.eos.debye.debye_fn_cheb(x)
        self.assertFloatEqual(test_debye, test_debye_cheb)

    def test_return_zero(self):
        rock = mypericlase()
        x = 0.
        test_helmholtz = burnman.eos.debye.helmholtz_free_energy(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_helmholtz, 0.)
        test_heat_capacity_v = burnman.eos.debye.molar_heat_capacity_v(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_heat_capacity_v, 0.)
        test_thermal_energy = burnman.eos.debye.thermal_energy(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_thermal_energy, 0.)

    def test_small(self):
        rock = mypericlase()
        x = 1e-16
        test_helmholtz = burnman.eos.debye.helmholtz_free_energy(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_helmholtz, 0.)
        test_heat_capacity_v = burnman.eos.debye.molar_heat_capacity_v(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_heat_capacity_v, 0.)
        test_thermal_energy = burnman.eos.debye.thermal_energy(
            x, rock.params['Debye_0'], rock.params['n'])
        self.assertFloatEqual(test_thermal_energy, 0.)

if __name__ == '__main__':
    unittest.main()
