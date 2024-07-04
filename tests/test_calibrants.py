from __future__ import absolute_import
from util import BurnManTest
import numpy as np
import unittest
import inspect
import burnman
from burnman import calibrants


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


if __name__ == "__main__":
    unittest.main()
