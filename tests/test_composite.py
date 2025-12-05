from util import BurnManTest
import warnings
import unittest

import burnman
from burnman import minerals, Composite, RelaxedComposite
from burnman.minerals import HP_2011_ds62, mp50MnNCKFMASHTO
from burnman.tools.eos import check_eos_consistency
from burnman.tools.equilibration import equilibrate


def setup_assemblage():
    """
    Set up a mineral assemblage for testing RelaxedComposite.
    """

    # Define observed mineral assemblage
    mu = mp50MnNCKFMASHTO.mu()
    bi = mp50MnNCKFMASHTO.bi()
    g = mp50MnNCKFMASHTO.g()
    ilmm = mp50MnNCKFMASHTO.ilmm()
    st = mp50MnNCKFMASHTO.st()
    q = HP_2011_ds62.q()

    # These compositions correspond to an equilibrated assemblage
    # at P = 0.4 GPa and T = 873 K
    compositions = [
        [0.55510672, 0.00393874, 0.00365819, 0.43232719, 0.00339884, 0.00157032],
        [
            0.12024528,
            0.55795845,
            -0.23022757,
            0.2928378,
            0.19457943,
            0.06260454,
            0.00200206,
        ],
        [0.11996805, 0.75282947, 0.10051702, 0.02013999, 0.00654547],
        [-0.05461279, 0.84687573, 0.1062428, 0.05404236, 0.04745189],
        [0.06925327, 0.80706976, 0.01346876, 0.03074615, 0.07946205],
    ]

    assemblage = Composite([mu, bi, g, ilmm, st, q], [1.0 / 6.0] * 6)
    for phase, composition in zip(assemblage.phases[:-1], compositions):
        phase.set_composition(composition)
    return assemblage


# TODO: test composite that changes number of entries
class composite(BurnManTest):
    def test_unroll(self):
        min1 = minerals.Murakami_etal_2012.fe_periclase()
        min2 = minerals.SLB_2005.periclase()

        c = burnman.Composite([min1], [1.0])
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        self.assertEqual(f, [0.0, 1.0])
        c = burnman.Composite([min1, min2], [0.4, 0.6])
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        self.assertEqual(f, [0.0, 0.4, 0.6])

        c1 = burnman.Composite([min1], [1.0])
        c2 = burnman.Composite([min2], [1.0])
        c = burnman.Composite([min1, c1, c2], [0.1, 0.4, 0.5])
        (m, f) = c.unroll()
        self.assertEqual(f, [0.0, 0.1, 0.0, 0.4, 0.5])

        min1 = burnman.minerals.Murakami_etal_2012.fe_periclase_HS()
        c1 = burnman.Composite([min1, min2], [0.1, 0.9])
        c2 = burnman.Composite([min2], [1.0])
        c = burnman.Composite([min1, c1, c2], [0.3, 0.1, 0.6])
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [0.3, 0.01, 0.09, 0.6])

    def test_changevalues(self):
        class mycomposite(burnman.Material):
            def __init__(self):
                burnman.Material.__init__(self)

            def unroll(self):
                fractions = [0.3, 0.7]
                mins = [
                    minerals.Murakami_etal_2012.fe_periclase(),
                    minerals.SLB_2005.periclase(),
                ]
                if self.temperature > 500:
                    fractions = [0.1, 0.9]
                return (mins, fractions)

        c = mycomposite()
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [0.3, 0.7])
        self.assertEqual(
            mins,
            ",".join(
                [
                    minerals.Murakami_etal_2012.fe_periclase().to_string(),
                    minerals.SLB_2005.periclase().to_string(),
                ]
            ),
        )
        c.set_state(5e9, 3000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [0.1, 0.9])
        self.assertEqual(
            mins,
            ",".join(
                [
                    minerals.Murakami_etal_2012.fe_periclase().to_string(),
                    minerals.SLB_2005.periclase().to_string(),
                ]
            ),
        )

    def test_number(self):
        class mycomposite(burnman.Material):
            def __init__(self):
                burnman.Material.__init__(self)

            def unroll(self):
                if self.temperature > 500:
                    return ([minerals.Murakami_etal_2012.fe_periclase()], [1.0])
                fractions = [0.3, 0.7]
                mins = [
                    minerals.Murakami_etal_2012.fe_periclase(),
                    minerals.SLB_2005.periclase(),
                ]
                return (mins, fractions)

        c = mycomposite()
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [1.0])
        self.assertEqual(mins, minerals.Murakami_etal_2012.fe_periclase().to_string())
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [0.3, 0.7])
        self.assertEqual(
            mins,
            ",".join(
                [
                    minerals.Murakami_etal_2012.fe_periclase().to_string(),
                    minerals.SLB_2005.periclase().to_string(),
                ]
            ),
        )

    def test_nest(self):
        min1 = minerals.Murakami_etal_2012.fe_periclase_LS()
        min2 = minerals.SLB_2005.periclase()
        ca = burnman.Composite([min1], [1.0])
        c = burnman.Composite([ca, min2], [0.4, 0.6])
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [0.4, 0.6])
        self.assertEqual(mins, ",".join([min1.to_string(), min2.to_string()]))

    def test_density_composite(self):
        pyrolite = burnman.Composite(
            [minerals.SLB_2005.mg_perovskite(), minerals.SLB_2005.periclase()],
            [0.95, 0.05],
        )
        pyrolite.set_method("slb3")
        pyrolite.set_state(40.0e9, 2000)

        d1 = pyrolite.phases[0].density
        d2 = pyrolite.phases[1].density
        dmix = pyrolite.density

        self.assertFloatEqual(d1, 4494.483)
        self.assertFloatEqual(d2, 4130.096)
        self.assertFloatEqual(dmix, 4486.263)

    def test_thermodynamic_potentials(self):
        rock1 = burnman.Composite(
            [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.periclase()],
            [1.0, 0.0],
        )
        rock2 = burnman.Composite(
            [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.periclase()],
            [0.0, 1.0],
        )
        rock3 = burnman.Composite(
            [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.periclase()],
            [0.4, 0.6],
        )
        mineral1 = minerals.SLB_2011.mg_perovskite()
        mineral2 = minerals.SLB_2011.periclase()

        rock1.set_state(40.0e9, 2000.0)
        rock2.set_state(40.0e9, 2000.0)
        rock3.set_state(40.0e9, 2000.0)
        mineral1.set_state(40.0e9, 2000.0)
        mineral2.set_state(40.0e9, 2000.0)

        # Gibbs
        gibbs1 = mineral1.molar_gibbs
        gibbs2 = mineral2.molar_gibbs
        G1 = rock1.molar_gibbs
        G2 = rock2.molar_gibbs
        G3 = rock3.molar_gibbs
        self.assertFloatEqual(G1, gibbs1)
        self.assertFloatEqual(G2, gibbs2)
        self.assertFloatEqual(G3, 0.4 * gibbs1 + 0.6 * gibbs2)

        # Helmholtz
        helmholtz1 = mineral1.molar_helmholtz
        helmholtz2 = mineral2.molar_helmholtz
        F1 = rock1.molar_helmholtz
        F2 = rock2.molar_helmholtz
        F3 = rock3.molar_helmholtz
        self.assertFloatEqual(F1, helmholtz1)
        self.assertFloatEqual(F2, helmholtz2)
        self.assertFloatEqual(F3, 0.4 * helmholtz1 + 0.6 * helmholtz2)

        # Enthalpy
        enthalpy1 = mineral1.molar_enthalpy
        enthalpy2 = mineral2.molar_enthalpy
        H1 = rock1.molar_enthalpy
        H2 = rock2.molar_enthalpy
        H3 = rock3.molar_enthalpy
        self.assertFloatEqual(H1, enthalpy1)
        self.assertFloatEqual(H2, enthalpy2)
        self.assertFloatEqual(H3, 0.4 * enthalpy1 + 0.6 * enthalpy2)

        # Energy
        internal_energy1 = mineral1.molar_internal_energy
        internal_energy2 = mineral2.molar_internal_energy
        U1 = rock1.molar_internal_energy
        U2 = rock2.molar_internal_energy
        U3 = rock3.molar_internal_energy
        self.assertFloatEqual(U1, internal_energy1)
        self.assertFloatEqual(U2, internal_energy2)
        self.assertFloatEqual(U3, 0.4 * internal_energy1 + 0.6 * internal_energy2)

        # Entropy
        molar_entropy1 = mineral1.molar_entropy
        molar_entropy2 = mineral2.molar_entropy
        S1 = rock1.molar_entropy
        S2 = rock2.molar_entropy
        S3 = rock3.molar_entropy
        self.assertFloatEqual(S1, molar_entropy1)
        self.assertFloatEqual(S2, molar_entropy2)
        self.assertFloatEqual(S3, 0.4 * molar_entropy1 + 0.6 * molar_entropy2)

    def test_summing_bigger(self):
        min1 = minerals.SLB_2005.periclase()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c = burnman.Composite([min1, min1], [0.8, 0.4])
            assert len(w) == 1
        c.set_method("slb3")
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [2.0 / 3.0, 1.0 / 3.0])

    def test_summing_smaller(self):
        min1 = minerals.SLB_2005.periclase()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c = burnman.Composite([min1, min1], [0.4, 0.2])
            assert len(w) == 1
        c.set_method("slb3")
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [2.0 / 3.0, 1.0 / 3.0])

    def test_summing_slightly_negative(self):
        min1 = minerals.SLB_2005.periclase()
        c = burnman.Composite([min1, min1, min1], [0.8, 0.2, 1.0 - 0.8 - 0.2])
        c.set_method("slb3")
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [0.8, 0.2, 0.0])
        self.assertTrue(f[2] >= 0.0)

    def test_swap_averaging(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.wuestite(), minerals.SLB_2005.mg_perovskite()],
            [0.5, 0.5],
        )

        rock.set_state(5.0e9, 1000.0)
        G1 = rock.effective_shear_modulus
        self.assertFloatEqual(108.686, G1 / 1.0e9)

        rock.set_averaging_scheme("HashinShtrikmanLower")
        rock.set_state(5.0e9, 1000.0)
        G2 = rock.effective_shear_modulus
        self.assertFloatEqual(104.451, G2 / 1.0e9)

        rock.set_averaging_scheme(burnman.averaging_schemes.HashinShtrikmanUpper())
        rock.set_state(5.0e9, 1000.0)
        G3 = rock.effective_shear_modulus
        self.assertFloatEqual(115.155, G3 / 1.0e9)

    def test_mass_to_molar_fractions(self):
        pv = burnman.minerals.HP_2011_ds62.mpv()
        en = burnman.minerals.HP_2011_ds62.en()
        c = burnman.Composite([pv, en], [0.5, 0.5], "mass")
        self.assertArraysAlmostEqual(c.molar_fractions, [2.0 / 3.0, 1.0 / 3.0])

    def test_mass_to_molar_fractions_2(self):
        min1 = minerals.SLB_2005.periclase()
        mass_fractions = [0.8, 0.2]
        c = burnman.Composite([min1, min1], mass_fractions, fraction_type="mass")
        self.assertArraysAlmostEqual(c.molar_fractions, mass_fractions)

    def test_formula(self):
        min1 = minerals.SLB_2005.mg_perovskite()
        min2 = minerals.SLB_2005.periclase()
        molar_fractions = [0.5, 0.5]
        c = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")
        self.assertArraysAlmostEqual(
            [c.formula["Mg"], c.formula["Si"], c.formula["O"]], [1.0, 0.5, 2.0]
        )

    def test_evaluate(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.wuestite(), minerals.SLB_2005.mg_perovskite()],
            [0.5, 0.5],
        )

        K, G = rock.evaluate(["K_eff", "G_eff"], 40.0e9, 2000.0)
        rock.set_state(40.0e9, 2000.0)
        K2 = rock.K_eff
        G2 = rock.G_eff
        self.assertFloatEqual(K, K2)
        self.assertFloatEqual(G, G2)

    def test_query_properties(self):
        min1 = minerals.SLB_2011.periclase()
        rock1 = burnman.Composite([min1, min1], [0.5, 0.5])
        rock1.set_state(min1.params["P_0"], min1.params["T_0"])

        self.assertFloatEqual(rock1.molar_volume, min1.params["V_0"])
        self.assertFloatEqual(rock1.molar_mass, min1.params["molar_mass"])
        self.assertFloatEqual(rock1.isothermal_bulk_modulus_reuss, min1.params["K_0"])
        self.assertFloatEqual(
            rock1.isothermal_compressibility_reuss, 1.0 / min1.params["K_0"]
        )
        self.assertFloatEqual(
            rock1.isentropic_compressibility_reuss,
            1.0 / rock1.isentropic_bulk_modulus_reuss,
        )
        self.assertFloatEqual(rock1.grueneisen_parameter, min1.params["grueneisen_0"])
        self.assertFloatEqual(rock1.thermal_expansivity, min1.alpha)
        self.assertFloatEqual(rock1.molar_heat_capacity_v, min1.C_v)
        self.assertFloatEqual(rock1.molar_heat_capacity_p, min1.C_p)

    def test_set_state_with_volume(self):
        min1 = minerals.SLB_2011.periclase()
        rock1 = burnman.Composite([min1, min1], [0.5, 0.5])
        P = 1.0e9
        P2 = 2.0e9
        T = 1000.0
        rock1.set_state(P, T)
        V = rock1.molar_volume
        rock1.set_state(P2, T)  # reset state
        _ = rock1.molar_volume
        rock1.set_state_with_volume(V, T)  # try to find original state

        self.assertFloatEqual(P, rock1.pressure)

    def test_stoichiometric_matrix_single_mineral(self):
        rock = burnman.Composite([minerals.SLB_2011.quartz()])

        self.assertTrue(len(rock.endmember_names) == 1)
        self.assertTrue(len(rock.elements) == 2)
        self.assertTrue(rock.elements[0] == "Si")
        self.assertTrue(rock.elements[1] == "O")
        self.assertTrue(rock.stoichiometric_matrix.shape == (1, 2))
        self.assertTrue(rock.reaction_basis.shape[0] == 0)
        self.assertArraysAlmostEqual(rock.compositional_null_basis[0], [-2, 1])
        self.assertTrue(rock.independent_element_indices[0] == 0)
        self.assertTrue(rock.dependent_element_indices[0] == 1)
        self.assertTrue(rock.n_reactions == 0)

    def test_stoichiometric_matrix_two_minerals(self):
        rock = burnman.Composite(
            [minerals.HP_2011_ds62.andr(), minerals.HP_2011_ds62.alm()]
        )

        self.assertTrue(len(rock.endmember_names) == 2)
        self.assertTrue(len(rock.elements) == 5)
        self.assertTrue(rock.stoichiometric_matrix.shape == (2, 5))
        self.assertTrue(rock.reaction_basis.shape[0] == 0)
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[0], [4.0 / 9.0, -2.0 / 3.0, 1, 0, 0]
        )
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[1], [-1.0 / 3.0, -1.0, 0, 1, 0]
        )
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[2], [-4.0 / 3.0, -4.0, 0, 0, 1]
        )
        self.assertArraysAlmostEqual(rock.independent_element_indices, [0, 1])
        self.assertArraysAlmostEqual(rock.dependent_element_indices, [2, 3, 4])
        self.assertTrue(rock.n_reactions == 0)

    def test_stoichiometric_matrix_binary_solution(self):
        rock = burnman.Composite([minerals.SLB_2011.mg_fe_olivine()])

        self.assertTrue(len(rock.endmember_names) == 2)
        self.assertTrue(len(rock.elements) == 4)
        self.assertTrue(rock.stoichiometric_matrix.shape == (2, 4))
        self.assertTrue(rock.reaction_basis.shape[0] == 0)
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[0], [-1.0 / 2.0, -1.0 / 2.0, 1, 0]
        )
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[1], [-2.0, -2.0, 0, 1]
        )
        self.assertArraysAlmostEqual(rock.independent_element_indices, [0, 1])
        self.assertArraysAlmostEqual(rock.dependent_element_indices, [2, 3])
        self.assertTrue(rock.n_reactions == 0)

    def test_stoichiometric_matrix_reaction(self):
        rock = burnman.Composite(
            [minerals.SLB_2011.mg_fe_olivine(), minerals.SLB_2011.garnet()]
        )

        self.assertTrue(len(rock.endmember_names) == 7)
        self.assertTrue(len(rock.elements) == 7)
        self.assertTrue(rock.stoichiometric_matrix.shape == (7, 7))
        self.assertArraysAlmostEqual(
            rock.reaction_basis[0], [3.0 / 2.0, -3.0 / 2.0, -1.0, 1.0, 0.0, 0.0, 0]
        )
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[0],
            [-1.0 / 2.0, -1.0, -1.0, -1.0, -3.0 / 2.0, -2.0, 1.0],
        )
        self.assertArraysAlmostEqual(
            rock.independent_element_indices, [0, 1, 2, 3, 4, 5]
        )
        self.assertArraysAlmostEqual(rock.dependent_element_indices, [6])
        self.assertTrue(rock.n_reactions == 1)

    def test_stoichiometric_matrix_reactions(self):
        rock = burnman.Composite(
            [
                minerals.SLB_2011.mg_fe_olivine(),
                minerals.SLB_2011.mg_fe_olivine(),
                minerals.SLB_2011.mg_fe_olivine(),
            ]
        )

        self.assertTrue(len(rock.endmember_names) == 6)
        self.assertTrue(len(rock.elements) == 4)
        self.assertTrue(rock.stoichiometric_matrix.shape == (6, 4))
        self.assertArraysAlmostEqual(rock.reaction_basis[0], [-1, 0, 1, 0, 0, 0])
        self.assertArraysAlmostEqual(rock.reaction_basis[1], [0, -1, 0, 1, 0, 0])
        self.assertArraysAlmostEqual(rock.reaction_basis[2], [-1, 0, 0, 0, 1, 0])
        self.assertArraysAlmostEqual(rock.reaction_basis[3], [0, -1, 0, 0, 0, 1])
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[0], [-1.0 / 2.0, -1.0 / 2.0, 1, 0]
        )
        self.assertArraysAlmostEqual(
            rock.compositional_null_basis[1], [-2.0, -2.0, 0, 1]
        )
        self.assertArraysAlmostEqual(rock.independent_element_indices, [0, 1])
        self.assertArraysAlmostEqual(rock.dependent_element_indices, [2, 3])
        self.assertTrue(rock.n_reactions == 4)

    def test_single_phase_composite(self):
        min1 = minerals.SLB_2011.mg_fe_olivine()
        min2 = minerals.SLB_2011.mg_fe_olivine()
        min1.set_composition([0.8, 0.2])
        min2.set_composition([0.8, 0.2])
        rock1 = burnman.Composite([min1, min2], [0.6, 0.4])

        properties = [
            "gibbs",
            "molar_entropy",
            "molar_volume",
            "molar_heat_capacity_p",
            "molar_heat_capacity_v",
            "thermal_expansivity",
            "density",
            "grueneisen_parameter",
            "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss",
        ]
        rock_prps = rock1.evaluate(properties, 1.0e9, 1400.0)
        min_prps = min1.evaluate(properties, 1.0e9, 1400.0)
        for i in range(len(properties)):
            self.assertFloatEqual(rock_prps[i], min_prps[i])

    def test_single_phase_composite_aliases(self):
        min1 = minerals.SLB_2011.mg_fe_olivine()
        min2 = minerals.SLB_2011.mg_fe_olivine()
        min1.set_composition([0.8, 0.2])
        min2.set_composition([0.8, 0.2])
        rock1 = burnman.Composite([min1, min2], [0.6, 0.4])

        properties = [
            "P",
            "T",
            "internal_energy",
            "helmholtz",
            "gibbs",
            "V",
            "rho",
            "S",
            "H",
            "K_T",
            "K_S",
            "beta_T",
            "beta_S",
            "gr",
            "alpha",
            "C_v",
            "C_p",
        ]
        rock_prps = rock1.evaluate(properties, 1.0e9, 1400.0)
        alias_prps = rock1.evaluate(properties, 1.0e9, 1400.0)
        for i in range(len(properties)):
            self.assertFloatEqual(rock_prps[i], alias_prps[i])

    def test_relaxed_composite(self):
        assemblage = setup_assemblage()
        old_formula = assemblage.formula.copy()
        old_number_of_moles = assemblage.number_of_moles
        relaxed_assemblage = RelaxedComposite(assemblage, assemblage.reaction_basis)
        relaxed_assemblage2 = RelaxedComposite(
            assemblage, assemblage.reaction_basis[:2, :]
        )

        P = 0.4e9
        T = 873.0
        assemblage.set_state(P, T)
        relaxed_assemblage.set_state(P, T)
        relaxed_assemblage2.set_state(P, T)

        # Check that the formula is unchanged by isochemical reaction
        for key in old_formula:
            self.assertAlmostEqual(old_formula[key], relaxed_assemblage.formula[key])
            self.assertAlmostEqual(old_formula[key], relaxed_assemblage2.formula[key])

        # Check that the number of moles has changed by isochemical reaction
        self.assertFalse(old_number_of_moles == relaxed_assemblage.number_of_moles)
        self.assertFalse(old_number_of_moles == relaxed_assemblage2.number_of_moles)

        self.assertTrue(
            check_eos_consistency(
                assemblage, P=P, T=T, including_shear_properties=False
            )
        )

        # The properties calculated for this relaxed assemblage should be EOS consistent
        self.assertTrue(
            check_eos_consistency(
                relaxed_assemblage, P=P, T=T, including_shear_properties=False
            )
        )

        # The properties for this assemblage are calculated by relaxing
        # only a subset of reactions,
        # so its properties should not be EOS consistent
        self.assertFalse(
            check_eos_consistency(
                relaxed_assemblage2, P=P, T=T, including_shear_properties=False
            )
        )

    def test_formula_extensivity_composite(self):
        min1 = minerals.SLB_2011.mg_perovskite()
        min2 = minerals.SLB_2011.periclase()
        molar_fractions = [0.5, 0.5]
        c1 = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")
        c2 = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")
        c2.number_of_moles = 2.0

        self.assertFloatEqual(c1.number_of_moles, c2.number_of_moles / 2.0)
        self.assertFloatEqual(c1.formula["Mg"], c2.formula["Mg"] / 2.0)
        self.assertFloatEqual(c1.formula["Si"], c2.formula["Si"] / 2.0)
        self.assertFloatEqual(c1.formula["O"], c2.formula["O"] / 2.0)
        self.assertFloatEqual(c1.molar_mass, c2.molar_mass)
        self.assertFloatEqual(c1.mass, c2.mass / 2.0)

    def test_property_extensivity_composite(self):
        min1 = minerals.SLB_2011.mg_perovskite()
        min2 = minerals.SLB_2011.periclase()
        molar_fractions = [0.5, 0.5]
        c1 = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")
        c2 = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")
        c2.number_of_moles = 2.0

        c1.set_state(40.0e9, 2000.0)
        c2.set_state(40.0e9, 2000.0)

        properties_extensive = [
            "internal_energy",
            "helmholtz",
            "gibbs",
            "enthalpy",
            "entropy",
            "volume",
            "heat_capacity_p",
            "heat_capacity_v",
            "mass",
            "V",
            "S",
            "H",
            "C_v",
            "C_p",
        ]

        for prop in properties_extensive:
            val1 = getattr(c1, prop)
            val2 = getattr(c2, prop)
            self.assertFloatEqual(val1, val2 / 2.0)

        properties_intensive = [
            "molar_internal_energy",
            "molar_helmholtz",
            "molar_gibbs",
            "molar_enthalpy",
            "molar_entropy",
            "molar_volume",
            "molar_heat_capacity_p",
            "molar_heat_capacity_v",
            "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss",
            "isothermal_compressibility_reuss",
            "isentropic_compressibility_reuss",
            "grueneisen_parameter",
            "thermal_expansivity",
            "molar_mass",
            "density",
        ]

        for prop in properties_intensive:
            val1 = getattr(c1, prop)
            val2 = getattr(c2, prop)
            self.assertFloatEqual(val1, val2)

    def test_composite_state_not_set(self):
        min1 = minerals.SLB_2011.mg_perovskite()
        min2 = minerals.SLB_2011.periclase()
        molar_fractions = [0.5, 0.5]
        c = burnman.Composite([min1, min2], molar_fractions, fraction_type="molar")

        with self.assertRaises(AttributeError):
            _ = c.density

    def test_formula_extensivity_relaxed_composite(self):
        fper = minerals.SLB_2024.ferropericlase()
        bdg = minerals.SLB_2024.bridgmanite()
        c1_unrelaxed = burnman.Composite([fper, bdg], [0.5, 0.5])
        c2_unrelaxed = burnman.Composite([fper, bdg], [0.5, 0.5])
        c2_unrelaxed.number_of_moles = 2.0

        fper.set_composition([0.35, 0.3, 0.2, 0.1, 0.05])
        bdg.set_composition([0.7, 0.1, 0.1, 0.05, 0.03, 0.01, 0.01])

        formula = c1_unrelaxed.formula.copy()
        formula2 = c2_unrelaxed.formula.copy()

        for element in ["Na", "Fe", "Mg", "Al", "Si", "O", "Cr"]:
            self.assertFloatEqual(formula[element], formula2[element] / 2.0)

        c1 = RelaxedComposite(c1_unrelaxed, c1_unrelaxed.reaction_basis)
        c2 = RelaxedComposite(c2_unrelaxed, c2_unrelaxed.reaction_basis)

        for element in ["Na", "Fe", "Mg", "Al", "Si", "O", "Cr"]:
            self.assertFloatEqual(c1.formula[element], c2.formula[element] / 2.0)
            self.assertFloatEqual(c1.unrelaxed.formula[element], c1.formula[element])
            self.assertFloatEqual(c2.unrelaxed.formula[element], c2.formula[element])

        c1.set_state(40.0e9, 2000.0)
        c2.set_state(40.0e9, 2000.0)

        self.assertArraysAlmostEqual(c1.molar_fractions, c2.molar_fractions)

        for element in ["Na", "Fe", "Mg", "Al", "Si", "O", "Cr"]:
            self.assertFloatEqual(c1.formula[element], c2.formula[element] / 2.0)
            self.assertFloatEqual(c1.unrelaxed.formula[element], c1.formula[element])
            self.assertFloatEqual(c2.unrelaxed.formula[element], c2.formula[element])

        self.assertFloatEqual(c1.number_of_moles, c2.number_of_moles / 2.0)
        self.assertFloatEqual(c1.molar_mass, c2.molar_mass)
        self.assertFloatEqual(c1.mass, c2.mass / 2.0)

    def test_property_extensivity_relaxed_composite(self):
        fper = minerals.SLB_2024.ferropericlase()
        bdg = minerals.SLB_2024.bridgmanite()
        c1_unrelaxed = burnman.Composite([fper, bdg], [0.5, 0.5])
        c2_unrelaxed = burnman.Composite([fper, bdg], [0.5, 0.5])
        c2_unrelaxed.number_of_moles = 2.0

        fper.set_composition([0.35, 0.3, 0.2, 0.1, 0.05])
        bdg.set_composition([0.7, 0.1, 0.1, 0.05, 0.03, 0.01, 0.01])

        formula = c1_unrelaxed.formula.copy()
        formula2 = c2_unrelaxed.formula.copy()

        self.assertFloatEqual(formula["Mg"] * 2.0, formula2["Mg"])

        c1 = RelaxedComposite(c1_unrelaxed, c1_unrelaxed.reaction_basis)
        c2 = RelaxedComposite(c2_unrelaxed, c2_unrelaxed.reaction_basis)

        c1.set_state(40.0e9, 2000.0)
        c2.set_state(40.0e9, 2000.0)

        properties_extensive = [
            "internal_energy",
            "helmholtz",
            "gibbs",
            "enthalpy",
            "entropy",
            "volume",
            "heat_capacity_p",
            "heat_capacity_v",
            "mass",
            "H",
            "S",
            "V",
            "C_p",
            "C_v",
        ]

        for prop in properties_extensive:
            val1 = getattr(c1, prop)
            val2 = getattr(c2, prop)
            self.assertFloatEqual(val1, val2 / 2.0)

        properties_intensive = [
            "molar_internal_energy",
            "molar_helmholtz",
            "molar_gibbs",
            "molar_enthalpy",
            "molar_entropy",
            "molar_volume",
            "molar_heat_capacity_p",
            "molar_heat_capacity_v",
            "isothermal_bulk_modulus_reuss",
            "isentropic_bulk_modulus_reuss",
            "isothermal_compressibility_reuss",
            "isentropic_compressibility_reuss",
            "grueneisen_parameter",
            "thermal_expansivity",
            "molar_mass",
            "density",
        ]

        for prop in properties_intensive:
            val1 = getattr(c1, prop)
            val2 = getattr(c2, prop)
            self.assertFloatEqual(val1, val2)

    def test_equilibrate_relaxed_composite(self):

        assemblage = setup_assemblage()
        relaxed_assemblage = RelaxedComposite(assemblage, assemblage.reaction_basis)
        equality_constraints = [["P", 0.5e9], ["T", 900.0]]
        equilibrate(
            relaxed_assemblage.formula, relaxed_assemblage, equality_constraints
        )
        f1 = relaxed_assemblage.molar_fractions.copy()

        assemblage = setup_assemblage()
        relaxed_assemblage = RelaxedComposite(assemblage, assemblage.reaction_basis)
        relaxed_assemblage.set_state(0.5e9, 900.0)
        f2 = relaxed_assemblage.molar_fractions.copy()

        self.assertArraysAlmostEqual(f1, f2)


if __name__ == "__main__":
    unittest.main()
