import unittest
import os, sys

sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.mineral import Mineral
from util import BurnManTest
from burnman.tools import *

from burnman.material import material_property



class test_decorators__material_property(BurnManTest):

    class MyCountingMaterial(burnman.Material):

        def __init__(self):
            burnman.Material.__init__(self)
            self.counter = 0

        @material_property
        def some_property(self):
            """my documentation"""
            self.counter += 1
            return 1.0

    def test_caching(self):
        m = self.MyCountingMaterial()
        self.assertEqual(m.counter, 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(m.counter, 1)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(m.counter, 1)

    def test_reset(self):
        m = self.MyCountingMaterial()
        self.assertEqual(len(m._cached), 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(len(m._cached), 1)

        m.reset()

        self.assertEqual(len(m._cached), 0)
        self.assertEqual(m.some_property, 1.0)
        self.assertEqual(len(m._cached), 1)

        self.assertEqual(m.counter, 2)

    def test_doc(self):
        """make sure documentation is passed through with the new decorator"""

        self.assertEqual(self.MyCountingMaterial.some_property.__doc__, "my documentation")

        internal_energy_doc = burnman.Material.internal_energy.__doc__
        self.assertTrue("internal energy of" in internal_energy_doc)

        self.assertEqual(self.MyCountingMaterial.internal_energy.__doc__,
                         burnman.Material.internal_energy.__doc__)

        pressure_doc = burnman.Material.pressure.__doc__
        self.assertTrue("pressure set" in pressure_doc)

        
if __name__ == '__main__':
    unittest.main()
