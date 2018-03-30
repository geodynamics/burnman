from __future__ import absolute_import
import unittest
import os
import sys
sys.path.insert(1, os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


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
                    minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                if (self.temperature > 500):
                    fractions = [0.1, 0.9]
                return (mins, fractions)

        c = mycomposite()
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [0.3, 0.7])
        self.assertEqual(
            mins, ",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(), minerals.SLB_2005.periclase().to_string()]))
        c.set_state(5e9, 3000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [0.1, 0.9])
        self.assertEqual(
            mins, ",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(), minerals.SLB_2005.periclase().to_string()]))

    def test_number(self):
        class mycomposite(burnman.Material):

            def __init__(self):
                burnman.Material.__init__(self)

            def unroll(self):
                if (self.temperature > 500):
                    return ([minerals.Murakami_etal_2012.fe_periclase()], [1.0])
                fractions = [0.3, 0.7]
                mins = [
                    minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                return (mins, fractions)

        c = mycomposite()
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [1.0])
        self.assertEqual(
            mins, minerals.Murakami_etal_2012.fe_periclase().to_string())
        c.set_state(5e9, 300)
        (m, f) = c.unroll()
        mins = ",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f, [0.3, 0.7])
        self.assertEqual(
            mins, ",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(), minerals.SLB_2005.periclase().to_string()]))

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
        pyrolite = burnman.Composite([minerals.SLB_2005.mg_perovskite(),
                                      minerals.SLB_2005.periclase()],
                                     [0.95, 0.05])
        pyrolite.set_method('slb3')
        pyrolite.set_state(40.e9, 2000)

        d1 = pyrolite.phases[0].density
        d2 = pyrolite.phases[1].density
        dmix = pyrolite.density

        self.assertFloatEqual(d1, 4494.483)
        self.assertFloatEqual(d2, 4130.096)
        self.assertFloatEqual(dmix, 4486.263)

    def test_thermodynamic_potentials(self):
        rock1 = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                                   minerals.SLB_2011.periclase()],
                                  [1.0, 0.0])
        rock2 = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                                   minerals.SLB_2011.periclase()],
                                  [0.0, 1.0])
        rock3 = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                                   minerals.SLB_2011.periclase()],
                                  [0.4, 0.6])
        mineral1 = minerals.SLB_2011.mg_perovskite()
        mineral2 = minerals.SLB_2011.periclase()

        rock1.set_state(40.e9, 2000.)
        rock2.set_state(40.e9, 2000.)
        rock3.set_state(40.e9, 2000.)
        mineral1.set_state(40.e9, 2000.)
        mineral2.set_state(40.e9, 2000.)

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
        self.assertFloatEqual(
            U3, 0.4 * internal_energy1 + 0.6 * internal_energy2)

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
        self.assertArraysAlmostEqual(f, [2. / 3., 1. / 3.])

    def test_summing_smaller(self):
        min1 = minerals.SLB_2005.periclase()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c = burnman.Composite([min1, min1], [0.4, 0.2])
            assert len(w) == 1
        c.set_method("slb3")
        c.set_state(5e9, 1000)
        (m, f) = c.unroll()
        self.assertArraysAlmostEqual(f, [2. / 3., 1. / 3.])

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
            [minerals.SLB_2005.wuestite(), minerals.SLB_2005.mg_perovskite()], [0.5, 0.5])

        rock.set_state(5.e9, 1000.)
        G1 = rock.shear_modulus
        self.assertFloatEqual(108.686, G1 / 1.e9)

        rock.set_averaging_scheme('HashinShtrikmanLower')
        rock.set_state(5.e9, 1000.)
        G2 = rock.shear_modulus
        self.assertFloatEqual(104.451, G2 / 1.e9)

        rock.set_averaging_scheme(
            burnman.averaging_schemes.HashinShtrikmanUpper())
        rock.set_state(5.e9, 1000.)
        G3 = rock.shear_modulus
        self.assertFloatEqual(115.155, G3 / 1.e9)

    def test_mass_to_molar_fractions(self):
        pv = burnman.minerals.HP_2011_ds62.mpv()
        en = burnman.minerals.HP_2011_ds62.en()
        c = burnman.Composite([pv, en], [0.5, 0.5], 'mass')
        self.assertArraysAlmostEqual(c.molar_fractions, [2. / 3., 1. / 3.])

    def test_mass_to_molar_fractions_2(self):
        min1 = minerals.SLB_2005.periclase()
        mass_fractions = [0.8, 0.2]
        c = burnman.Composite(
            [min1, min1], mass_fractions, fraction_type='mass')
        self.assertArraysAlmostEqual(c.molar_fractions, mass_fractions)

    def test_velocities_from_rock(self):

        rock = burnman.Composite(
            [minerals.SLB_2005.wuestite(), minerals.SLB_2005.mg_perovskite()], [0.5, 0.5])

        # Use default averaging schemes
        _, _, _, _, K, G = burnman.velocities_from_rock(
            rock, [40.e9, ], [2000., ])
        rock.set_state(40.e9, 2000.)
        K2 = rock.K_S
        G2 = rock.G
        self.assertFloatEqual(K, K2)
        self.assertFloatEqual(G, G2)

        # Now change the averaging schemes
        _, _, _, _, K, G = burnman.velocities_from_rock(
            rock, [40.e9, ], [2000., ], averaging_scheme=burnman.averaging_schemes.HashinShtrikmanAverage())
        rock.set_averaging_scheme('HashinShtrikmanAverage')
        rock.set_state(40.e9, 2000.)
        K2 = rock.K_S
        G2 = rock.G
        self.assertFloatEqual(K, K2)
        self.assertFloatEqual(G, G2)

    def test_query_properties(self):
        min1 = minerals.SLB_2011.periclase()
        rock1 = burnman.Composite([min1, min1], [0.5, 0.5])
        rock1.set_state(min1.params['P_0'], min1.params['T_0'])

        self.assertFloatEqual(rock1.molar_volume, min1.params['V_0'])
        self.assertFloatEqual(rock1.molar_mass, min1.params['molar_mass'])
        self.assertFloatEqual(
            rock1.isothermal_bulk_modulus, min1.params['K_0'])
        self.assertFloatEqual(
            rock1.isothermal_compressibility, 1. / min1.params['K_0'])
        self.assertFloatEqual(
            rock1.adiabatic_compressibility, 1. / rock1.adiabatic_bulk_modulus)
        self.assertFloatEqual(
            rock1.grueneisen_parameter, min1.params['grueneisen_0'])
        self.assertFloatEqual(rock1.thermal_expansivity, min1.alpha)
        self.assertFloatEqual(rock1.molar_heat_capacity_v, min1.C_v)
        self.assertFloatEqual(rock1.molar_heat_capacity_p, min1.C_p)

if __name__ == '__main__':
    unittest.main()
