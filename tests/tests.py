import unittest
from test_vrh import *
from test_spin import *

import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals




#class MyTest(unittest.TestCase):
#    def test(self):
#        self.assertEqual(3, 4)



#class CompareL2(unittest.TestCase):    

#class VsVp(unittest.TestCase):


if __name__ == '__main__':
    unittest.main()
