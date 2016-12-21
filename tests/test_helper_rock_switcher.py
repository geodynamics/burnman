from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import numpy as np
import burnman
from burnman import minerals
from burnman import mineral_helpers

from util import BurnManTest




class RockSwitcher(BurnManTest):

    min_lm = burnman.minerals.SLB_2011.mg_bridgmanite()
    min_um = burnman.minerals.SLB_2011.forsterite()

    class min_other(burnman.Mineral):
        def __init__(self):
            self.params = {
                'name': 'something',
                'equation_of_state': 'slb3',
                'T_0': 300.,
                'P_0': 0.,
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
            burnman.Mineral.__init__(self)

    def test_low_high_switches(self):
        rock = mineral_helpers.HelperLowHighPressureRockTransition(25e9,self.min_lm, self.min_um)
        #print rock.to_string()

        T0 = 1500
        rock.set_state(20e9, T0)
        assert(rock.current_rock == self.min_lm)
        rock.set_state(30e9, T0)
        assert(rock.current_rock == self.min_um)

    def test_tripple(self):

        class tripple(mineral_helpers.HelperRockSwitcher):
            def __init__(self, rocks):
                self.rocks = rocks
                mineral_helpers.HelperRockSwitcher.__init__(self)

            def select_rock(self):
                if self._pressure < 20e9:
                    return self.rocks[0]
                elif self._pressure < 30e9:
                    return self.rocks[1]
                else:
                    return self.rocks[2]

            def set_method(self, method):
                for r in self.rocks:
                    r.set_method(method)

        rock = tripple([self.min_other(), self.min_lm, self.min_um])

        pressures = [10e9, 25e9, 35e9]
        molar_volumes = []
        for p in pressures:
            rock.set_state(p, 300)
            molar_volumes.append(rock.molar_volume)

        self.assertArraysAlmostEqual(molar_volumes, [1.062933864e-05, 2.2478565706e-05, 3.62424836698e-05])

    def test_composite(self):
        c1 = burnman.Composite([self.min_other(), self.min_lm], [0.4, 0.6])
        rock = mineral_helpers.HelperLowHighPressureRockTransition(50e9, c1, self.min_um)

        rock.set_state(40e9, 300)
        assert(rock.current_rock == c1)
        rock.set_state(60e9, 300)
        assert(rock.current_rock == self.min_um)

    def test_properties(self):
        rock = mineral_helpers.HelperLowHighPressureRockTransition(25e9,self.min_lm, self.min_um)
        vars = ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G']
        output = rock.evaluate(vars, [20e9, 40e9], [300, 300])

        output1 = self.min_lm.evaluate(vars, [20e9], [300])
        output2 = self.min_um.evaluate(vars, [40e9], [300])
        ref = np.concatenate((output1, output2), axis=1)
        assert(np.array_equal(output, ref))




if __name__ == '__main__':
    unittest.main()
