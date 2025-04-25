from __future__ import absolute_import
from util import BurnManTest
import numpy as np
import unittest
import inspect
import burnman
from burnman import calibrants
from burnman.calibrants.tools import pressure_to_pressure


class test_calibrants(BurnManTest):

    def test_calibrant_initialisation(self):
        pressures = []

        libs = dir(calibrants)
        for lib in libs:

            objectgroup = getattr(calibrants, lib)
            if objectgroup.__class__.__name__ == "module":
                for m in dir(objectgroup):
                    obj = getattr(objectgroup, m)
                    if (
                        inspect.isclass(obj)
                        and obj != burnman.Calibrant
                        and issubclass(obj, burnman.Calibrant)
                    ):

                        calibrant = obj()
                        V = calibrant.volume(10.0e9, 1000.0)
                        pressures.append(calibrant.pressure(V, 1000.0))

        pressures = np.array(pressures)
        self.assertArraysAlmostEqual(pressures, 10.0e9 + pressures * 0.0)

    def test_pressure_to_pressure(self):
        gold1 = calibrants.Anderson_1989.Au()
        gold2 = calibrants.Fei_2007.Au()

        temperature = 1000.0
        pressure1 = 20.0e9
        pressure2 = pressure_to_pressure(gold1, gold2, pressure1, temperature)
        pressure1new = pressure_to_pressure(gold2, gold1, pressure2, temperature)

        self.assertAlmostEqual(pressure1 / 1.0e9, pressure1new / 1.0e9, 4)


if __name__ == "__main__":
    unittest.main()
