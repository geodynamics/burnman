from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals
from burnman import seismic
from burnman.planet import Planet
from burnman.layer import Layer
from burnman import mineral_helpers as helpers

from util import BurnManTest
import numpy as np


def make_simple_planet():
    core = Layer("core", np.linspace(0.e3,3480.e3,10))
    core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
    core.set_temperature_mode('user-defined', temperatures= 300.*np.ones_like(core.radii))
    
    mantle = Layer("mantle", np.linspace(3480.e3, 6371.e3, 10))
    mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
    mantle.set_temperature_mode('adiabatic', temperature_top = 1200)
    myplanet = Planet('earth_like',[core, mantle])
    myplanet.make()
    return myplanet, core, mantle


class test_planet(BurnManTest):
    def test_planet_1(self):
        myplanet, core, mantle = make_simple_planet()
        assert(myplanet.get_layer("core") == core)
        assert(myplanet.get_layer_by_radius(2000.e3) == core)
        assert(myplanet.get_layer_by_radius(4000.e3) == mantle)
        assert(myplanet.get_layer_by_radius(5000.e3) == mantle)
        self.assertFloatEqual(myplanet.pressure[0], 445219931774.84, tol=1.e-4)

    def test_sort_layers(self):
        myplanet, core, mantle = make_simple_planet()
        assert(myplanet.layers == [core, mantle])

    def test_temperature(self):
        myplanet, core, mantle = make_simple_planet()
        self.assertArraysAlmostEqual(mantle.S, [mantle.S[-1] for S in mantle.S])

    def test_evaluate(self):
        myplanet, core, mantle = make_simple_planet()
        alpha, rho = myplanet.evaluate(['alpha', 'rho'], myplanet.radii)
        self.assertFloatEqual(alpha[-2], myplanet.alpha[-2])
        self.assertFloatEqual(rho[-2], myplanet.density[-2])

        
if __name__ == '__main__':
    unittest.main()
