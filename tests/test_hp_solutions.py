import unittest
from util import BurnManTest
from burnman import Mineral, CombinedMineral, Solution
from burnman.constants import gas_constant
from burnman.classes.solutionmodel import IdealSolution
from burnman.minerals import (
    HP_2011_fluids,
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
    def test_ig_fluid_activities_and_water_endmember(self):
        fluid = ig50NCKFMASHTOCr.fl()
        water_index = fluid.endmember_names.index("H2O")
        self.assertEqual(fluid.n_endmembers, 11)
        self.assertIsInstance(fluid.endmembers[water_index][0], HP_2011_fluids.H2O)
        proportions = np.full(fluid.n_endmembers, 0.05)
        proportions[water_index] = 0.5
        expected = (1.0 - proportions[water_index]) * proportions
        expected[water_index] = proportions[water_index] ** 2
        self.assertArraysAlmostEqual(
            IdealSolution.activities(fluid.solution_model, 1.0e9, 1000.0, proportions),
            expected,
        )
        for pure in np.eye(fluid.n_endmembers):
            self.assertFloatEqual(
                fluid.solution_model.excess_gibbs_energy(1.0e9, 1000.0, pure),
                0.0,
                tol_zero=1.0e-9,
            )
        fluid.set_composition(proportions)
        fluid.set_state(1.0e9, 1000.0)
        self.assertTrue(np.isfinite(fluid.molar_gibbs))
        self.assertGreater(fluid.molar_volume, 0.0)
        self.assertTrue(np.isfinite(fluid.molar_heat_capacity_p))
        pure_water = np.zeros(fluid.n_endmembers)
        pure_water[water_index] = 1.0
        fluid.set_composition(pure_water)
        water = HP_2011_fluids.H2O()
        water.set_state(1.0e9, 1000.0)
        self.assertFloatEqual(fluid.molar_gibbs, water.molar_gibbs)
        self.assertFloatEqual(fluid.molar_volume, water.molar_volume)

    def test_withdrawn_ig_melts_are_excluded(self):
        for dataset in (ig50NCKFMASHTOCr, ig50NCKFMASTOCr):
            with self.subTest(dataset=dataset.__name__):
                self.assertFalse(hasattr(dataset, "liq"))
                self.assertIn("algebraic error", dataset.__doc__)

    def test_melt_ideal_activities(self):
        for dataset, constructor in (
            (mb50NCKFMASHTO, mb50NCKFMASHTO.L),
            (mp50KFMASH, mp50KFMASH.liq),
            (mp50MnNCKFMASHTO, mp50MnNCKFMASHTO.liq),
            (mp50NCKFMASHTO, mp50NCKFMASHTO.liq),
        ):
            with self.subTest(dataset=dataset.__name__):
                melt = constructor()
                names = melt.endmember_names
                water = names.index("h2oL")
                magnesium = names.index("fo2L")
                iron = names.index("fa2L")
                proportions = np.arange(1.0, melt.n_endmembers + 1.0)
                proportions /= sum(proportions)
                dry_fraction = 1.0 - proportions[water]
                olivine_fraction = proportions[magnesium] + proportions[iron]
                expected = dry_fraction * proportions
                expected[water] = proportions[water] ** 2
                for i in (magnesium, iron):
                    expected[i] = (
                        dry_fraction
                        * olivine_fraction
                        * (proportions[i] / olivine_fraction) ** 5
                    )
                model = melt.solution_model
                self.assertArraysAlmostEqual(
                    IdealSolution.activities(model, 1.0e9, 1500.0, proportions),
                    expected,
                )
                self.assertArraysAlmostEqual(
                    IdealSolution.excess_partial_entropies(
                        model, 1.0e9, 1500.0, proportions
                    ),
                    -gas_constant * np.log(expected),
                )
                # At every pure endmember the mixing Gibbs energy vanishes,
                # including sites whose multiplicity becomes zero.
                for pure in np.eye(melt.n_endmembers):
                    self.assertFloatEqual(
                        model.excess_gibbs_energy(1.0e9, 1500.0, pure),
                        0.0,
                        tol_zero=1.0e-9,
                    )
                melt.set_composition(proportions)
                melt.set_state(1.0e9, 1500.0)
                self.assertTrue(np.isfinite(melt.molar_gibbs))
                self.assertTrue(melt.molar_volume > 0.0)

    def test_metabasite_melt_speciation_correction(self):
        melt = mb50NCKFMASHTO.L()
        anorthite = melt.endmembers[melt.endmember_names.index("anoL")][0]
        self.assertEqual(anorthite.formula, {"Ca": 1, "Al": 2, "Si": 2, "O": 8})
        self.assertArraysAlmostEqual(
            list(anorthite.property_modifiers[0][1].values()),
            [-46500.0, 0.0, -2.5e-6],
        )

    def test_staurolite_dqf_correction(self):
        for dataset in (mp50KFMASH, mp50MnNCKFMASHTO, mp50NCKFMASHTO):
            with self.subTest(dataset=dataset.__name__):
                self.assertEqual(
                    dataset.mstm.property_modifiers[0][1]["delta_E"], -8000.0
                )

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
