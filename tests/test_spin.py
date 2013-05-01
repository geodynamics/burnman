import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals



class spin_transistion(unittest.TestCase):
    def test_new(self):
        
        mins = [minerals.Muretal2012.fe_periclase(), minerals.Muretal2012.fe_periclase_HS(), minerals.Muretal2012.fe_periclase_LS()]
        for p in mins:
            p.set_method('slb2')
        
        #print "HS regime: (on/high/low)"
        for p in mins:
            p.set_state(5e9, 300)
            #print p.v_s()
         
        self.assertAlmostEqual(mins[0].v_s(), mins[1].v_s())
        
        #print "LS regime: (on/high/low)"
        for p in mins:
            p.set_state(70e9, 300)
            #print p.v_s()
        
        self.assertAlmostEqual(mins[0].v_s(), mins[2].v_s())


if __name__ == '__main__':
    unittest.main()
