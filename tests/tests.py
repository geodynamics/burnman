from __future__ import absolute_import
import unittest
from util import BurnManTest
from util import *

from test_anisotropy import *
from test_averaging import *
from test_chemical_potentials import *
from test_composite import *
from test_debye import *
from test_decorators import *
from test_endmembers import *
from test_eos import *
from test_eos_consistency import *
from test_fitting import *
from test_geotherm import *
from test_helper_rock_switcher import *
from test_material import *
from test_minerals import *
from test_model import *
from test_modifiers import *
from test_partitioning import *
from test_perplex import *
from test_planet import *
from test_seismic import *
from test_solidsolution import *
from test_solvers import *
from test_spin import *
from test_tools import *

import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals


class TestRock(BurnManTest):

    def test_rock(self):
        amount_perovskite = 0.3
        rock = burnman.Composite(
            [minerals.SLB_2005.mg_perovskite(), minerals.SLB_2005.periclase()],
            [amount_perovskite, 1.0 - amount_perovskite])
        (phases, fr) = rock.unroll()
        self.assertFloatEqual(fr[0], 0.3)
        self.assertFloatEqual(fr[1], 0.7)


if __name__ == '__main__':
    unittest.main()
