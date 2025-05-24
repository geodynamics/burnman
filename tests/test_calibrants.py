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
        np.testing.assert_allclose(pressures, 10.0e9 + pressures * 0.0, rtol=1e-5)

    def test_pressure_to_pressure(self):
        gold1 = calibrants.Anderson_1989.Au()
        gold2 = calibrants.Fei_2007.Au()

        temperature = 1000.0
        pressure1 = 20.0e9
        pressure2 = pressure_to_pressure(gold1, gold2, pressure1, temperature)
        pressure1new = pressure_to_pressure(gold2, gold1, pressure2, temperature)

        np.testing.assert_allclose(pressure1 / 1.0e9, pressure1new / 1.0e9, rtol=1e-5)

    def test_pressure_to_pressure_cov(self):
        gold1 = calibrants.Anderson_1989.Au()
        gold2 = calibrants.Fei_2007.Au()

        temperature = 1000.0
        pressure1 = 20.0e9
        PTcov1 = np.array([[1.0e18, 1.0e8], [1.0e8, 10.0]])
        pressure2, PTcov2 = pressure_to_pressure(
            gold1, gold2, pressure1, temperature, PTcov1
        )
        pressure1new, PTcov1new = pressure_to_pressure(
            gold2, gold1, pressure2, temperature, PTcov2
        )

        np.testing.assert_allclose(PTcov1, PTcov1new, rtol=1e-5)
        np.testing.assert_allclose(pressure1 / 1.0e9, pressure1new / 1.0e9, rtol=1e-5)


if __name__ == "__main__":
    unittest.main()
