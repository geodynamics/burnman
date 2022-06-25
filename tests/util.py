from __future__ import absolute_import
import unittest


class BurnManTest(unittest.TestCase):

    def assertFloatEqual(self, a, b, tol=1e-5, tol_zero=1e-16):
        self.assertAlmostEqual(
            a, b, delta=max(tol_zero, max(abs(a), abs(b)) * tol))

    def assertArraysAlmostEqual(self, a, b, tol=1e-5, tol_zero=1e-16):
        self.assertEqual(len(a), len(b))
        for (i1, i2) in zip(a, b):
            self.assertFloatEqual(i1, i2, tol, tol_zero)


class Huh(BurnManTest):

    def test(self):
        self.assertFloatEqual(5200.01, 5200.015)

    def test_to_zero(self):
        self.assertFloatEqual(0, 1e-20)

if __name__ == '__main__':
    unittest.main()
