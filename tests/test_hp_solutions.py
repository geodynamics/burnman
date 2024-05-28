from __future__ import absolute_import
from __future__ import print_function
import unittest
from util import BurnManTest
from burnman import Mineral, CombinedMineral, Solution
from burnman.minerals import ig50NCKFMASHTOCr, ig50NCKFMASTOCr
import inspect
import numpy as np
from collections import Counter


class hp_solutions(BurnManTest):
    def test_ig50NCKFMASHTOCr(self):
        for m in dir(ig50NCKFMASHTOCr):
            mineral = getattr(ig50NCKFMASHTOCr, m)
            if (
                inspect.isclass(mineral)
                and mineral != Mineral
                and mineral != CombinedMineral
                and issubclass(mineral, Mineral)
            ):
                if issubclass(mineral, Solution) and mineral is not Solution:
                    mbr_names = mineral().endmember_names
                    n_mbrs = len(mbr_names)
                    molar_fractions = np.ones(n_mbrs) / n_mbrs
                    m = mineral()
                    m.set_composition(molar_fractions)
                    self.assertTrue(type(m.formula) is Counter)

    def test_ig50NCKFMASTOCr(self):
        for m in dir(ig50NCKFMASTOCr):
            mineral = getattr(ig50NCKFMASTOCr, m)
            if (
                inspect.isclass(mineral)
                and mineral != Mineral
                and mineral != CombinedMineral
                and issubclass(mineral, Mineral)
            ):
                if issubclass(mineral, Solution) and mineral is not Solution:
                    mbr_names = mineral().endmember_names
                    n_mbrs = len(mbr_names)
                    molar_fractions = np.ones(n_mbrs) / n_mbrs
                    m = mineral()
                    m.set_composition(molar_fractions)
                    self.assertTrue(type(m.formula) is Counter)


if __name__ == "__main__":
    unittest.main()
