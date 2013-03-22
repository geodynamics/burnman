import unittest
from test_vrh import *
from test_spin import *

import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


class TestRock(unittest.TestCase):
    def test_rock(self):
        amount_perovskite = 0.3
        rock = burnman.composite( ( ( minerals.mg_fe_perovskite(0.1), amount_perovskite ), 
                                    (minerals.ferropericlase(0.2), 1.0-amount_perovskite) ) )
        rock.set_method('slb')
        self.assertAlmostEqual(rock.phases[0].fraction, 0.3, 2)
        self.assertAlmostEqual(rock.phases[1].fraction, 0.7, 2)


#class MyTest(unittest.TestCase):
#    def test(self):
#        self.assertEqual(3, 4)



#class CompareL2(unittest.TestCase):    

#class VsVp(unittest.TestCase):


if __name__ == '__main__':
    unittest.main()
