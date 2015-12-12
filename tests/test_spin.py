from __future__ import absolute_import
import unittest
from util import BurnManTest
import os, sys
sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals



class spin_transition(BurnManTest):
    def test_new(self):
        
        mins = [minerals.Murakami_etal_2012.fe_periclase(), minerals.Murakami_etal_2012.fe_periclase_HS(), minerals.Murakami_etal_2012.fe_periclase_LS()]
        for p in mins:
            p.set_method('slb2')
        
        #print "HS regime: (on/high/low)"
        for p in mins:
            p.set_state(5e9, 300)
            #print p.v_s()
         
        c,f = mins[0].unroll()
        self.assertFloatEqual(c[0].v_s, mins[1].v_s)
        
        #print "LS regime: (on/high/low)"
        for p in mins:
            p.set_state(70e9, 300)
            #print p.v_s()
        
        c,f = mins[0].unroll()
        self.assertFloatEqual(c[0].v_s, mins[2].v_s)

    def test_no_set_state(self):
        m = minerals.Murakami_etal_2012.fe_periclase()
        m.set_state(5e9, 300)
        self.assertIsInstance(m.unroll()[0][0], minerals.Murakami_etal_2012.fe_periclase_HS)
        m.set_state(70e9, 300)
        self.assertIsInstance(m.unroll()[0][0], minerals.Murakami_etal_2012.fe_periclase_LS)


class TestHelperSolidSolution(BurnManTest):
    def test1(self):
        m = minerals.other.ferropericlase(0.1)
        m.set_state(5e9, 300)
        self.assertFloatEqual(m.v_s, 5821.42007777)
        m = minerals.other.ferropericlase(0.7)
        m.set_state(5e9, 300)
        self.assertFloatEqual(m.v_s, 4061.92139873)


class TestHelperFEdep(BurnManTest):
    def test(self):
        weight_percents = {'Mg':0.213, 'Fe': 0.08, 'Si':0.27, 'Ca':0., 'Al':0.}
        Kd_0 = .59 #Fig 5 Nakajima et al 2012

        phase_fractions, relative_molar_percent = burnman. \
            calculate_phase_percents(weight_percents)
        iron_content = lambda p,t: burnman.calculate_partition_coefficient \
                (p,t,relative_molar_percent,Kd_0)

        rock = burnman.Composite([minerals.SLB_2005.mg_fe_perovskite_pt_dependent(iron_content,0),
                                  minerals.SLB_2005.ferropericlase_pt_dependent(iron_content,1)], \
                                 [phase_fractions['pv'], phase_fractions['fp']])
        rock.set_state(5e9, 300)
        mins, fractions = rock.unroll()
        self.assertArraysAlmostEqual(fractions, [0.9428714062806316, 0.057128593719368403])
        self.assertIsInstance(mins[0], minerals.SLB_2005.mg_fe_perovskite)
        self.assertIsInstance(mins[1], minerals.SLB_2005.ferropericlase)
        self.assertFloatEqual(mins[0].molar_mass, 0.101752790682)

        rock.set_state(7e9, 700)
        mins, fractions = rock.unroll()
        self.assertFloatEqual(mins[0].molar_mass, 0.104161162508)



if __name__ == '__main__':
    unittest.main()
