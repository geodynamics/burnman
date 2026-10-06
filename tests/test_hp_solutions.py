import unittest
from util import BurnManTest
from burnman import Mineral, CombinedMineral, Solution
from burnman.minerals import (
    ig50NCKFMASHTOCr,
    ig50NCKFMASTOCr,
    mb50NCKFMASHTO,
    mp50KFMASH,
    mp50MnNCKFMASHTO,
    mp50NCKFMASHTO,
)

import inspect
import numpy as np
from collections import Counter


class hp_solutions(BurnManTest):
    def test_mb_acmm_endmembers(self):
        augite = mb50NCKFMASHTO.aug()
        diopside = mb50NCKFMASHTO.dio()
        aug_acmm = augite.endmembers[augite.endmember_names.index("acmm")][0]
        dio_acmm = diopside.endmembers[diopside.endmember_names.index("acmm")][0]
        self.assertIsNot(aug_acmm, dio_acmm)
        self.assertEqual(aug_acmm.property_modifiers[0][1]["delta_E"], -5000.0)
        self.assertEqual(dio_acmm.property_modifiers[0][1]["delta_E"], -7000.0)

    def test_ig_cfm_endmembers(self):
        for dataset in (ig50NCKFMASHTOCr, ig50NCKFMASTOCr):
            with self.subTest(dataset=dataset.__name__):
                olivine = dataset.ol()
                clinopyroxene = dataset.cpx()
                ol_cfm = olivine.endmembers[olivine.endmember_names.index("cfm")][0]
                cpx_cfm = clinopyroxene.endmembers[
                    clinopyroxene.endmember_names.index("cfm")
                ][0]
                self.assertIsNot(ol_cfm, cpx_cfm)
                self.assertEqual(ol_cfm.formula, {"Mg": 1, "Fe": 1, "Si": 1, "O": 4})
                self.assertEqual(cpx_cfm.formula, {"Mg": 1, "Fe": 1, "Si": 2, "O": 6})
                self.assertArraysAlmostEqual(ol_cfm.mixture.molar_fractions, [0.5, 0.5])
                self.assertArraysAlmostEqual(
                    cpx_cfm.mixture.molar_fractions, [0.5, 0.5]
                )
                self.assertEqual(
                    ol_cfm.property_modifiers[0][1],
                    {"delta_E": 0.0, "delta_S": 0.0, "delta_V": 0.0},
                )
                self.assertArraysAlmostEqual(
                    list(cpx_cfm.property_modifiers[0][1].values()),
                    [-1600.0, 2.0, 4.65e-7],
                )

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

    def test_mb50NCKFMASHTO(self):
        for m in dir(mb50NCKFMASHTO):
            mineral = getattr(mb50NCKFMASHTO, m)
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

    def test_mp50KFMASH(self):
        for m in dir(mp50KFMASH):
            mineral = getattr(mp50KFMASH, m)
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

    def test_mp50MnNCKFMASHTO(self):
        for m in dir(mp50MnNCKFMASHTO):
            mineral = getattr(mp50MnNCKFMASHTO, m)
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

    def test_mp50NCKFMASHTO(self):
        for m in dir(mp50NCKFMASHTO):
            mineral = getattr(mp50NCKFMASHTO, m)
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
