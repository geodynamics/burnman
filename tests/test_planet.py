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


class test_planet(BurnManTest):

    def test_planet_1(self):

        core = Layer("core", radius_planet= 6371.e3, min_depth=2890.e3, max_depth =6371.e3, n_slices= 10)
        core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
        core.set_temperature_mode('user_defined', temperatures= 300.*np.ones_like(core.depths))
        
        mantle = Layer("mantle", radius_planet= 6371.e3, min_depth=0e3, max_depth =2890.e3,  n_slices= 10)
        mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        mantle.set_temperature_mode('adiabat', temperature_top = 1200)
        myplanet = Planet('earth_like',[core, mantle])
        myplanet.set_state()
        
        assert(myplanet.get_layer("core") == core)
        assert(myplanet.get_layer_by_radius(2000.e3) == core)
        assert(myplanet.get_layer_by_radius(4000.e3) == mantle)
        assert(myplanet.get_layer_by_radius(5000.e3) == mantle)

        self.assertFloatEqual(myplanet.pressure[-1], 445428905300.0, tol=1.e-4)


    def test_sort_layers(self):
        core = Layer("core", radius_planet= 6371.e3, min_depth=2890.e3, max_depth =6371.e3, n_slices= 10)
        core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
        core.set_temperature_mode('user_defined', temperatures= 300.*np.ones_like(core.depths))
        
        mantle = Layer("mantle", radius_planet= 6371.e3, min_depth=0e3, max_depth =2890.e3,  n_slices= 10)
        mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        mantle.set_temperature_mode('adiabat', temperature_top = 1200)
        myplanet = Planet('earth_like',[core, mantle])
        assert(myplanet.layers == [mantle,core])

    def test_temperature(self):
        core = Layer("core", radius_planet= 6371.e3, min_depth=2890.e3, max_depth =6371.e3, n_slices= 10)
        core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
        core.set_temperature_mode('user-defined', temperatures= burnman.geotherm.brown_shankland(core.depths))
        
        mantle = Layer("mantle", radius_planet= 6371.e3, min_depth=0e3, max_depth =2890.e3,  n_slices= 10)
        mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        mantle.set_temperature_mode('adiabat', temperature_top = 1200)
        myplanet = Planet('earth_like',[core, mantle])
        myplanet.set_state()
        
        ref = [ 1200.,          1304.56273992,  1397.12842763,  1480.88446903,  1557.96838056,
           1629.98958458,  1698.27053255,  1763.96944242,  1828.14624527,  1891.79862543,
           2452.25581395,  2694.29333333,  2902.27777778,  3073.52666667,  3207.49444444,
           3305.10666667,  3369.93333333,  3433.17333333,  3470.79333333,  3484.        ]
        print(myplanet.temperature)
        print(ref)
        self.assertArraysAlmostEqual(myplanet.temperature, ref)


if __name__ == '__main__':
    unittest.main()
