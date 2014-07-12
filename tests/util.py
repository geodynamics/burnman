import unittest

class BurnManTest(unittest.TestCase):
    def assertArraysAlmostEqual(self, a, b):
        self.assertEqual(len(a), len(b))
        for (i1, i2) in zip(a, b):
            self.assertAlmostEqual(i1, i2)

