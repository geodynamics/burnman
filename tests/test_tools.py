import unittest
from util import BurnManTest
import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman import Composite
from burnman.tools.chemistry import equilibrium_temperature
from burnman.tools.chemistry import equilibrium_pressure
from burnman.tools.chemistry import hugoniot
from burnman.tools.chemistry import reactions_from_stoichiometric_matrix
from burnman.tools.chemistry import reactions_from_formulae
from burnman.tools.plot import plot_projected_elastic_properties, pretty_plot
from burnman.utils.chemistry import convert_fractions
from burnman.utils.math import bracket
from burnman.utils.math import smooth_array
from burnman.utils.math import interp_smoothed_array_and_derivatives
from burnman.utils.math import _pad_ndarray_inverse_mirror
from burnman.utils.math import is_positive_definite

from burnman.minerals import HP_2011_ds62, mp50MnNCKFMASHTO
from burnman.tools.thermobarometry import estimate_conditions


def setup_assemblage(scale=1.0):
    """
    Set up a mineral assemblage for testing thermobarometry.
    :param scale: Scale factor for compositional noise.
    """
    np.random.seed(42)

    # Define observed mineral assemblage
    mu = mp50MnNCKFMASHTO.mu()
    bi = mp50MnNCKFMASHTO.bi()
    g = mp50MnNCKFMASHTO.g()
    ilmm = mp50MnNCKFMASHTO.ilmm()
    st = mp50MnNCKFMASHTO.st()
    q = HP_2011_ds62.q()

    # These compositions correspond to an equilibrated assemblage at P = 0.4 GPa and T = 873 K
    # In practice, these would be the measured compositions of the minerals
    compositions = [
        np.array(
            [0.55510672, 0.00393874, 0.00365819, 0.43232719, 0.00339884, 0.00157032]
        ),
        np.array(
            [
                0.12024528,
                0.55795845,
                -0.23022757,
                0.2928378,
                0.19457943,
                0.06260454,
                0.00200206,
            ]
        ),
        np.array([0.11996805, 0.75282947, 0.10051702, 0.02013999, 0.00654547]),
        np.array([-0.05461279, 0.84687573, 0.1062428, 0.05404236, 0.04745189]),
        np.array([0.06925327, 0.80706976, 0.01346876, 0.03074615, 0.07946205]),
    ]

    # Assumed compositional uncertainties (1 sigma)
    mu.compositional_covariances = [0.01, 0.001, 0.001, 0.01, 0.001, 0.0005]
    bi.compositional_covariances = [0.01, 0.01, 0.05, 0.01, 0.01, 0.01, 0.0005]
    g.compositional_covariances = [0.01, 0.01, 0.01, 0.01, 0.001]
    ilmm.compositional_covariances = [0.01, 0.01, 0.01, 0.01, 0.01]
    st.compositional_covariances = [0.01, 0.01, 0.01, 0.01, 0.01]

    assemblage = Composite([mu, bi, g, ilmm, st, q], [1.0 / 6.0] * 6)
    for phase, comp in zip(assemblage.phases[:-1], compositions):
        noise = np.random.normal(0, 1.0, len(comp)) * phase.compositional_covariances
        comp_with_noise = comp + noise * scale
        # Normalize to ensure valid molar fractions
        comp_with_noise /= np.sum(comp_with_noise)
        phase.set_composition(comp_with_noise)

    # Convert uncertainties to covariances
    for phase in assemblage.phases[:-1]:
        phase.compositional_covariances = np.diag(
            np.array(phase.compositional_covariances) ** 2
        )

    return assemblage


class test_tools(BurnManTest):
    def test_eqm_T(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        H_ex = 1.0e3
        S_ex = 2.0
        fo2.params["H_0"] = fo.params["H_0"] + H_ex
        fo2.params["S_0"] = fo.params["S_0"] + S_ex

        T = H_ex / S_ex
        T_calc = equilibrium_temperature(
            [fo, fo2], [1.0, -1.0], fo.params["P_0"], 1200.0
        )

        self.assertArraysAlmostEqual([T], [T_calc])

    def test_eqm_P(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        fo2 = burnman.minerals.HP_2011_ds62.fo()

        fo.params["K_0"] = 1.0e20
        fo2.params["K_0"] = 1.0e20

        H_ex = 1
        V_ex = 1.0e-8
        fo2.params["H_0"] = fo.params["H_0"] + H_ex
        fo2.params["V_0"] = fo.params["V_0"] - V_ex

        P = fo.params["P_0"] + H_ex / V_ex
        P_calc = equilibrium_pressure([fo, fo2], [1.0, -1.0], fo.params["T_0"])
        self.assertArraysAlmostEqual([P], [P_calc])

    def test_hugoniot(self):
        fo = burnman.minerals.HP_2011_ds62.fo()
        T_ref = 298.15
        P_ref = 1.0e5
        pressures = np.array([P_ref])
        temperatures, volumes = hugoniot(fo, P_ref, T_ref, pressures)

        self.assertArraysAlmostEqual(temperatures, np.array([T_ref]))

    def test_fraction_conversion_mass(self):
        pv = burnman.minerals.HP_2011_ds62.mpv()
        en = burnman.minerals.HP_2011_ds62.en()
        c = burnman.Composite([pv, en], [0.5, 0.5])
        molar_fractions = convert_fractions(c, [0.5, 0.5], "mass", "molar")
        self.assertArraysAlmostEqual(molar_fractions, [2.0 / 3.0, 1.0 / 3.0])

    def test_fraction_conversion_mass_2(self):
        per = burnman.minerals.SLB_2011.periclase()
        stv = burnman.minerals.SLB_2011.stishovite()
        fo = burnman.minerals.SLB_2011.forsterite()
        c = burnman.Composite([per, stv, fo], [0.5, 0.25, 0.25])
        mass_fractions = convert_fractions(c, [0.5, 0.25, 0.25], "molar", "mass")
        self.assertArraysAlmostEqual(
            [mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]]
        )

    def test_fraction_conversion_volume(self):
        per = burnman.minerals.SLB_2011.periclase()
        c = burnman.Composite([per, per, per], [0.25, 0.25, 0.5])
        c.set_state(1.0e5, 300.0)
        mass_fractions = convert_fractions(c, [0.25, 0.25, 0.5], "molar", "volume")
        self.assertArraysAlmostEqual(
            [mass_fractions[0] + mass_fractions[1]], [mass_fractions[2]]
        )

    def test_bracket(self):
        def fn(x):
            return (x - 1.0) * (x - 2.0)

        # Test that it can find the root at one from the right
        sol = bracket(fn, 1.2, 1.0e-2)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[1] < 1.0)
        self.assertTrue(sol[0] > 1.0)

        # Test that it can find the root at one from the left
        sol = bracket(fn, 0.6, 1.0e-2)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[0] < 1.0)
        self.assertTrue(sol[1] > 1.0)

        # Test that it can find the root at two from the left
        sol = bracket(fn, 1.7, 1.0e-2)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[0] < 2.0)
        self.assertTrue(sol[1] > 2.0)

        # Test that it can find the root at two from the right
        sol = bracket(fn, 2.5, 1.0e-2)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[1] < 2.0)
        self.assertTrue(sol[0] > 2.0)

        # Test that it can find the root at one if we take too large of a dx
        sol = bracket(fn, 1.2, 2.0)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[1] < 1.0)
        self.assertTrue(sol[0] > 1.0)

        # Test that it can find the root at two if we take too large of a dx
        sol = bracket(fn, 1.8, 2.0)
        self.assertTrue(sol[2] * sol[3] < 0.0)
        self.assertTrue(sol[0] < 2.0)
        self.assertTrue(sol[1] > 2.0)

    def test_bracket_failure(self):
        mineral = burnman.minerals.SLB_2011.fayalite()
        # This should be too high pressure for the EoS
        mineral.set_state(300.0e9, 300.0)

        def fn():
            return mineral.molar_volume

        with np.errstate(all="ignore"):
            self.assertRaises(Exception, fn)

    def test_padding_1D(self):
        array = np.array([1.0, 2.0, 3.0, 5.0])
        padding = (1,)
        padded_array = _pad_ndarray_inverse_mirror(array, padding)
        self.assertArraysAlmostEqual(
            padded_array, np.array([0.0, 1.0, 2.0, 3.0, 5.0, 7.0])
        )

    def test_padding_2D(self):
        array = np.array([[1.0, 2.0], [3.0, 5.0]])
        padding = (1, 1)
        padded_array = _pad_ndarray_inverse_mirror(array, padding)
        self.assertArraysAlmostEqual(
            np.ndarray.flatten(padded_array),
            np.ndarray.flatten(
                np.array(
                    [
                        [-3.0, -1.0, -1.0, 1.0],
                        [0.0, 1.0, 2.0, 3.0],
                        [1.0, 3.0, 5.0, 7.0],
                        [4.0, 5.0, 8.0, 9.0],
                    ]
                )
            ),
        )

    def test_smoothing(self):
        array = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        smoothed_array = smooth_array(array, [1.0], [1.0])
        self.assertArraysAlmostEqual(array, smoothed_array)

    def test_interp_smoothing_ij(self):
        array = np.array([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])
        axis_values = np.array([0.0, 1.0, 2.0])
        f, dfdx, dfdy = interp_smoothed_array_and_derivatives(
            array, axis_values, axis_values, 0, 0, indexing="ij"
        )
        self.assertArraysAlmostEqual(
            [f([0.0, 1.0])[0], dfdx([0.0, 1.0])[0], dfdy([0.0, 1.0])[0]],
            [array[0][1], 0.0, 1.0],
        )

    def test_interp_smoothing_xy(self):
        array = np.array([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])
        axis_values = np.array([0.0, 1.0, 2.0])
        f, dfdx, dfdy = interp_smoothed_array_and_derivatives(
            array, axis_values, axis_values, 0, 0, indexing="xy"
        )
        self.assertArraysAlmostEqual(
            [f([0.0, 1.0])[0], dfdx([0.0, 1.0])[0], dfdy([0.0, 1.0])[0]],
            [array[1][0], 1.0, 0.0],
        )

    def test_n_reactions_from_stoichiometry(self):
        M = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 1.0]])
        R = reactions_from_stoichiometric_matrix(M)
        self.assertTrue(len(R) == 6.0)

    def test_reactions_from_formulae(self):
        formulae = ["MgO", "FeO", "Fe2Si2O6", "Fe2SiO4"]
        compound_names = ["per", "fper", "fs", "fa"]
        R = sorted(
            reactions_from_formulae(formulae, compound_names, return_strings=True)
        )
        self.assertTrue(R[0] == "2 fa = 2 fper + 1 fs")

    def test_positive_definite(self):
        arr = np.array([[2.0, 0.0], [0.0, 1.0]])
        self.assertTrue(is_positive_definite(arr))

    def test_plot_elastic(self):

        pretty_plot()

        stishovite_Q1 = burnman.Mineral(
            {
                "name": "Stishovite",
                "formula": {"Si": 1.0, "O": 2.0},
                "equation_of_state": "slb3",
                "F_0": -815061.947383201,
                "V_0": 1.3985103364121215e-05,
                "K_0": 302651021060.0,
                "Kprime_0": 4.231768387000001,
                "Debye_0": 1099.41845635,
                "grueneisen_0": 1.815179405,
                "q_0": 2.2095517465427,
                "G_0": 228000000000.0,
                "Gprime_0": 1.94045,
                "eta_s_0": 4.40394,
                "n": 3.0,
                "Z": 2.0,
                "molar_mass": 0.060085,
                "T_0": 300.0,
                "F_0": 0.0,
                "P_0": 0.0,
            }
        )
        a = 2.754879e-02
        b = 2.839605e-02
        c = 1.3985103364121215e-05 / (a * b)
        cell_parameters = [a, b, c, 90, 90, 90]

        c11 = 7.93215475
        c44 = 3.6433892
        c66 = 4.63773363
        d11 = 0.03176523
        d44 = 0.31665864
        d66 = 0.09046939

        carr = c11 * np.ones((6, 6))
        carr[3, 3] = c44
        carr[4, 4] = c44
        carr[5, 5] = c66
        darr = d11 * np.ones((6, 6))
        darr[3, 3] = d44
        darr[4, 4] = d44
        darr[5, 5] = d66

        anisotropic_parameters = {
            "a": np.array(
                [
                    [0.56655959, -0.08994629, -0.04885663, 0.0, 0.0, 0.0],
                    [-0.08994629, 0.71555477, -0.25117643, 0.0, 0.0, 0.0],
                    [-0.04885663, -0.25117643, 0.49784435, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 1.33005083, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 1.45154562, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, 0.97764314],
                ]
            ),
            "b": np.array(
                [
                    [0.17631275, -0.06562267, -0.04243526, 0.0, 0.0, 0.0],
                    [-0.06562267, 0.17631275, -0.04243526, 0.0, 0.0, 0.0],
                    [-0.04243526, -0.04243526, -0.0516391, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, -0.70121757, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0, -0.70121757, 0.0],
                    [0.0, 0.0, 0.0, 0.0, 0.0, -0.18867047],
                ]
            ),
            "c": carr,
            "d": darr,
        }

        def psi_func(f, Pth, params):
            a = params["a"]
            b = params["b"]
            c = params["c"]
            d = params["d"]

            dPsidf = a + b * np.tanh(c * f + d)
            Psi = a * f + b / c * np.log(np.cosh(c * f + d) / np.cosh(d))
            dPsidPth = np.zeros((6, 6))
            return (Psi, dPsidf, dPsidPth)

        m = burnman.AnisotropicMineral(
            stishovite_Q1,
            cell_parameters,
            anisotropic_parameters,
            psi_function=psi_func,
            orthotropic=True,
        )
        m.set_state(1.0e9, 1000.0)
        plot_types = [
            "vp",
            "vs1",
            "vs2",
            "linear compressibility",
            "minimum poisson ratio",
            "maximum poisson ratio",
            "vp/vs1",
            "vp/vs2",
            "s anisotropy",
            "linear compressibility",
            "youngs modulus",
        ]

        fig = plt.figure(figsize=(12, 7))
        ax = [fig.add_subplot(4, 3, i, projection="polar") for i in range(1, 12)]

        # very low resolution to make sure the function works
        contour_sets, ticks, lines = plot_projected_elastic_properties(
            m, plot_types, ax, 4, 4, 5, 5
        )

        self.assertTrue(len(contour_sets) == 11)
        self.assertTrue(len(ticks) == 11)
        self.assertTrue(len(lines) == 11)

    def test_estimate_conditions(self):
        assemblage = setup_assemblage()
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()
        res = estimate_conditions(
            assemblage, dataset_covariances, guessed_conditions=theta
        )
        self.assertTrue(res.success)
        self.assertEqual(res.n_reactions, 18)

    def test_estimate_conditions_at_equilibrium(self):
        assemblage = setup_assemblage(scale=0.0)
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()
        res = estimate_conditions(
            assemblage, dataset_covariances, guessed_conditions=theta
        )
        self.assertTrue(res.success)
        self.assertArraysAlmostEqual(
            [res.x[0] / 1.0e6, res.x[1]], [400.0, 873.0], tol=0.1
        )

    def test_estimate_conditions_exclude_small_fractions(self):
        assemblage = setup_assemblage()
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()
        res = estimate_conditions(
            assemblage,
            dataset_covariances,
            guessed_conditions=theta,
            small_fraction_tol=0.01,
        )
        self.assertTrue(res.success)
        self.assertEqual(res.n_reactions, 11)

    def test_estimate_conditions_fixed_P(self):
        assemblage = setup_assemblage()
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()
        res = estimate_conditions(
            assemblage,
            dataset_covariances,
            guessed_conditions=theta,
            pressure_bounds=[1.0e9, 1.0e9],
        )
        self.assertTrue(res.success)
        self.assertEqual(res.x[0], 1.0e9)

    def test_estimate_conditions_fixed_T(self):
        assemblage = setup_assemblage()
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()
        res = estimate_conditions(
            assemblage,
            dataset_covariances,
            guessed_conditions=theta,
            temperature_bounds=[500.0, 500.0],
        )
        self.assertTrue(res.success)
        self.assertEqual(res.x[1], 500.0)

    def test_estimate_conditions_fail_fixed_P_and_T(self):
        assemblage = setup_assemblage()
        theta = np.array([0.5e9, 500.0])
        dataset_covariances = HP_2011_ds62.cov()

        # assert that an exception is raised when both P and T are fixed
        with self.assertRaises(Exception):
            _ = estimate_conditions(
                assemblage,
                dataset_covariances,
                guessed_conditions=theta,
                pressure_bounds=[1.0e9, 1.0e9],
                temperature_bounds=[500.0, 500.0],
            )


if __name__ == "__main__":
    unittest.main()
