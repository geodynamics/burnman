import inspect
import unittest

import numpy as np

from burnman import Solution
from burnman.classes.solutionmodel import IdealSolution
from burnman.constants import gas_constant
from burnman.minerals import HPx_ds636, ig51W24, ig51G25
from burnman.tools.eos import check_eos_consistency


class ig51_solutions(unittest.TestCase):
    def test_complete_solution_sets(self):
        expected = {
            ig51W24: {
                "liq_W24d": 14,
                "ol_H18": 4,
                "g_W24": 6,
                "spl_T21": 8,
                "fsp_H22": 3,
                "cpx_W24": 10,
                "opx_W24": 9,
                "ilm_W24": 5,
                "nph_W24": 6,
                "kals_W24": 2,
                "lct_W24": 2,
                "mel_W24": 5,
            },
            ig51G25: {
                "liq_G25w": 12,
                "fl_G25": 11,
                "ol_H18": 4,
                "g_W24": 6,
                "spl_T21": 8,
                "fsp_H22": 3,
                "cpx_W24": 10,
                "opx_W24": 9,
                "ilm_W24": 5,
                "bi_G25": 6,
                "mu_W14": 6,
                "amp_G16": 11,
                "ep_H11": 3,
                "cd_G25": 3,
            },
        }
        for module, phases in expected.items():
            actual = {
                name: constructor
                for name, constructor in vars(module).items()
                if inspect.isclass(constructor)
                and constructor.__module__ == module.__name__
                and issubclass(constructor, Solution)
            }
            self.assertEqual(set(actual), set(phases))
            for name, constructor in actual.items():
                with self.subTest(dataset=module.__name__, phase=name):
                    solution = constructor()
                    self.assertEqual(solution.n_endmembers, phases[name])
                    solution.set_composition(np.full(phases[name], 1.0 / phases[name]))
                    solution.set_state(1.0e9, 1400.0)
                    self.assertTrue(np.isfinite(solution.molar_gibbs))
                    self.assertGreater(solution.molar_volume, 0.0)
                    self.assertTrue(np.isfinite(solution.molar_heat_capacity_p))

    def test_corrected_melt_activities(self):
        for constructor in (ig51W24.liq_W24d, ig51G25.liq_G25w):
            with self.subTest(melt=constructor.__name__):
                melt = constructor()
                model = melt.solution_model
                # The authors' melt equations mix Mg, Fe, Ca and Al on the
                # microscopic site; olivine has four occupants on this site.
                p = np.arange(1.0, melt.n_endmembers + 1.0)
                p /= sum(p)
                mg, fe, sl, wo = 3, 4, 1, 2
                total = 4.0 * (p[mg] + p[fe]) + p[sl] + p[wo]
                expected = p.copy()
                expected[sl] *= p[sl] / total
                expected[wo] *= p[wo] / total
                for i in (mg, fe):
                    expected[i] = (p[mg] + p[fe]) * (4.0 * p[i] / total) ** 4
                if constructor is ig51G25.liq_G25w:
                    expected *= 1.0 - p[-1]
                    expected[-1] = p[-1] ** 2
                np.testing.assert_allclose(
                    IdealSolution.activities(model, 1.0e9, 1400.0, p),
                    expected,
                    rtol=1.0e-12,
                    atol=1.0e-15,
                )
                np.testing.assert_allclose(
                    IdealSolution.excess_partial_entropies(model, 1.0e9, 1400.0, p),
                    -gas_constant * np.log(expected),
                    rtol=1.0e-12,
                )
                for pure in np.eye(melt.n_endmembers):
                    self.assertAlmostEqual(
                        model.excess_gibbs_energy(1.0e9, 1400.0, pure), 0.0
                    )
                melt.set_composition(p)
                self.assertTrue(
                    check_eos_consistency(melt, P=1.0e9, T=1400.0, tol=1.0e-4)
                )

    def test_g25_fluid_water_and_activities(self):
        fluid = ig51G25.fl_G25()
        water_index = fluid.endmember_names.index("H2O")
        self.assertIsInstance(fluid.endmembers[water_index][0], HPx_ds636.H2O)
        p = np.full(fluid.n_endmembers, 0.05)
        p[water_index] = 0.5
        expected = (1.0 - p[water_index]) * p
        expected[water_index] = p[water_index] ** 2
        np.testing.assert_allclose(
            IdealSolution.activities(fluid.solution_model, 1.0e9, 1200.0, p),
            expected,
        )
        pure = np.zeros(fluid.n_endmembers)
        pure[water_index] = 1.0
        fluid.set_composition(pure)
        fluid.set_state(1.0e9, 1200.0)
        water = HPx_ds636.H2O()
        water.set_state(1.0e9, 1200.0)
        self.assertAlmostEqual(fluid.molar_gibbs, water.molar_gibbs)
        self.assertAlmostEqual(fluid.molar_volume, water.molar_volume)

    def test_temperature_and_pressure_interactions(self):
        # W(py,gr) = 45.4 - 0.010*T + 0.04*P [kJ/mol; P in kbar].
        for module in (ig51W24, ig51G25):
            model = module.g_W24().solution_model
            p = np.array([0.5, 0.0, 0.5, 0.0, 0.0, 0.0])
            # With v(py)=1 and v(gr)=2.5, equal molar proportions have a
            # Van Laar interaction weight of 10/49.
            self.assertEqual(model.alphas[0], 1.0)
            self.assertEqual(model.alphas[2], 2.5)
            weight = 10.0 / 49.0
            for pressure, temperature in ((1.0e5, 1000.0), (1.0e9, 1400.0)):
                excess = model.excess_gibbs_energy(pressure, temperature, p)
                ideal = IdealSolution(model.endmembers).excess_gibbs_energy(
                    pressure, temperature, p
                )
                expected = weight * (45400.0 - 10.0 * temperature + 4.0e-7 * pressure)
                self.assertAlmostEqual(excess - ideal, expected)
                self.assertAlmostEqual(
                    model.excess_entropy(pressure, temperature, p)
                    - IdealSolution(model.endmembers).excess_entropy(
                        pressure, temperature, p
                    ),
                    weight * 10.0,
                )
                self.assertAlmostEqual(
                    model.excess_volume(pressure, temperature, p), weight * 4.0e-7
                )

    def test_buffer_stoichiometry_and_energies(self):
        for module in (ig51W24, ig51G25):
            for label, recipe in (
                ("fmq", (("fa", -3.0), ("q", 3.0), ("mt", 2.0))),
                ("nno", (("NiO", 2.0), ("Ni", -2.0))),
            ):
                with self.subTest(dataset=module.__name__, buffer=label):
                    buffer = getattr(module, label)
                    self.assertEqual(buffer.formula, {"O": 2.0})
                    buffer.set_state(1.0e9, 1200.0)
                    expected = 0.0
                    for endmember, amount in recipe:
                        mineral = getattr(HPx_ds636, endmember)()
                        mineral.set_state(1.0e9, 1200.0)
                        expected += amount * mineral.molar_gibbs
                    self.assertAlmostEqual(buffer.molar_gibbs, expected)

    def test_duplicated_endmember_labels(self):
        for module in (ig51W24, ig51G25):
            olivine = module.ol_H18()
            clinopyroxene = module.cpx_W24()
            ol_cfm = olivine.endmembers[olivine.endmember_names.index("cfm")][0]
            cpx_cfm = clinopyroxene.endmembers[
                clinopyroxene.endmember_names.index("cfm")
            ][0]
            self.assertIsNot(ol_cfm, cpx_cfm)
            self.assertEqual(ol_cfm.formula["Si"], 1.0)
            self.assertEqual(cpx_cfm.formula["Si"], 2.0)
        fluid = ig51G25.fl_G25()
        melt = ig51G25.liq_G25w()
        fluid_fo = fluid.endmembers[fluid.endmember_names.index("fofL")][0]
        melt_fo = melt.endmembers[melt.endmember_names.index("fo2L")][0]
        self.assertIsNot(fluid_fo, melt_fo)
        self.assertEqual(fluid_fo.property_modifiers[0][1]["delta_E"], 7560.0)
        self.assertEqual(melt_fo.property_modifiers[0][1]["delta_E"], 7740.0)


if __name__ == "__main__":
    unittest.main()
