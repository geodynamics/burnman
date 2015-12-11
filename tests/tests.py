from __future__ import absolute_import
import unittest
from util import BurnManTest
from util import *
from test_averaging import *
from test_spin import *
from test_composite import *
from test_model import *
from test_partitioning import *
from test_eos import *
from test_debye import *
from test_geotherm import *
from test_endmembers import *
from test_solidsolution import *
from test_tools import *


import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


class TestRock(BurnManTest):
    def test_rock(self):
        amount_perovskite = 0.3
        rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite], \
            [minerals.SLB_2005.mg_fe_perovskite(0.1), minerals.SLB_2005.ferropericlase(0.2)] )
        (fr,phases)=rock.unroll()
        self.assertFloatEqual(fr[0], 0.3)
        self.assertFloatEqual(fr[1], 0.7)




if __name__ == '__main__':
    unittest.main()
