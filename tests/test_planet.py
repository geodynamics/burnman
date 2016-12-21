from __future__ import absolute_import
import unittest
import os
import sys

sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals
from burnman import seismic
from burnman import planet
from burnman import mineral_helpers as helpers

from util import BurnManTest


class test_planet(BurnManTest):

    def test_planet_1(self):

        core = planet.Planet.Layer("core", burnman.minerals.other.Liquid_Fe_Anderson(), 3485e3, 10)
        LM = burnman.minerals.SLB_2011.mg_bridgmanite()
        UM = burnman.minerals.SLB_2011.forsterite()
        mantle_rock = helpers.HelperLowHighPressureRockTransition(25.0e9, UM, LM)
        mantle = planet.Planet.Layer("mantle", mantle_rock, 6371e3, 10)
        myplanet = planet.Planet([core, mantle])

        assert(myplanet.get_layer("core") == core)
        assert(myplanet.get_layer_by_radius(2000e3) == core)
        assert(myplanet.get_layer_by_radius(4000e3) == mantle)
        assert(myplanet.get_layer_by_radius(5000e3) == mantle)

        self.assertFloatEqual(myplanet.pressures[0], 438580390453.7)


    def test_sort_layers(self):
        core = planet.Planet.Layer("core", burnman.minerals.other.Liquid_Fe_Anderson(), 3485e3, 10)
        LM = burnman.minerals.SLB_2011.mg_bridgmanite()
        UM = burnman.minerals.SLB_2011.forsterite()
        mantle_rock = helpers.HelperLowHighPressureRockTransition(25.0e9, UM, LM)
        mantle = planet.Planet.Layer("mantle", mantle_rock, 6371e3, 10)
        myplanet = planet.Planet([mantle, core])
        assert(myplanet.layers == [core, mantle])

    def test_temperature(self):
        core = planet.Planet.Layer("core", burnman.minerals.other.Liquid_Fe_Anderson(), 200.0, 1000.0, 10)
        mantle_lower = planet.Planet.LayerLinearTemperature("lower mantle", burnman.minerals.SLB_2011.mg_bridgmanite(), 800.0, 1000, 300.0, 10)
        mantle_upper = planet.Planet.Layer("upper mantle", burnman.minerals.SLB_2011.forsterite(), 1000.0, 400.0, 10)
        myplanet = planet.Planet([core, mantle_lower, mantle_upper])

        ref = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 998.83333333333337,
               922.22222222222217, 844.44444444444434, 766.66666666666663, 688.88888888888891, 611.11111111111109,
               533.33333333333326, 455.55555555555543, 377.77777777777771, 300.0, 400.0, 400.0, 400.0, 400.0, 400.0,
               400.0, 400.0, 400.0, 400.0, 400.0]

        self.assertArraysAlmostEqual(myplanet.temperatures, ref)


if __name__ == '__main__':
    unittest.main()
