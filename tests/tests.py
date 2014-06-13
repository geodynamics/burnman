import unittest
from test_vrh import *
from test_spin import *
from test_composite import *

import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


class TestRock(unittest.TestCase):
    def test_rock(self):
        amount_perovskite = 0.3
        rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite], \
            [minerals.SLB_2005.mg_fe_perovskite(0.1), minerals.SLB_2005.ferropericlase(0.2)] )
        rock.set_method('slb2')
        (fr,phases)=rock.unroll()
        self.assertAlmostEqual(fr[0], 0.3, 2)
        self.assertAlmostEqual(fr[1], 0.7, 2)


#class MyTest(unittest.TestCase):
#    def test(self):
#        self.assertEqual(3, 4)



#class CompareL2(unittest.TestCase):    

#class VsVp(unittest.TestCase):


if __name__ == '__main__':
    unittest.main(verbosity=2)
