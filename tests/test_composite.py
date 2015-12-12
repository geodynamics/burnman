from __future__ import absolute_import
import unittest
import os, sys
sys.path.insert(1,os.path.abspath('..'))
import warnings

import burnman
from burnman import minerals

from util import BurnManTest


# TODO: test composite that changes number of entries


class composite(BurnManTest):

    def test_unroll(self):

        min1 = minerals.Murakami_etal_2012.fe_periclase()
        min_hs = burnman.minerals.Murakami_etal_2012.fe_periclase_HS()
        min2 = minerals.SLB_2005.periclase()

        c = burnman.Composite( [min1], [1.0] )
        c.set_state(5e9,300)
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertEqual(f,[1.0])
        self.assertEqual(mins,min_hs.to_string())
        c = burnman.Composite( [min1, min2], [0.4, 0.6] )
        c.set_state(5e9,300)
        (m,f) = c.unroll()
        self.assertEqual(f,[0.4,0.6])

        c1 = burnman.Composite( [min1], [1.0] )
        c2 = burnman.Composite( [min2], [1.0] )
        c = burnman.Composite( [min1, c1, c2], [0.1, 0.4, 0.5] )
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertEqual(f,[0.1,0.4,0.5])
        self.assertEqual(mins,min_hs.to_string()+','+min_hs.to_string()+','+min2.to_string())

        c1 = burnman.Composite( [min1, min2], [0.1, 0.9] )
        c2 = burnman.Composite( [min2], [1.0] )
        c = burnman.Composite( [min1, c1, c2], [0.3, 0.1, 0.6] )
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.01,0.09,0.6])
        self.assertEqual(mins,",".join([min_hs.to_string(),min_hs.to_string(),min2.to_string(),min2.to_string()]))
        
        #c.set_method("slb2")
        #c.set_state(5e9,300)

    def test_changevalues(self):
        class mycomposite(burnman.Material):
            def __init__(self):
                burnman.Material.__init__(self)

            def unroll(self):
                fractions = [0.3, 0.7]
                mins = [minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                if (self.temperature>500):
                    fractions = [0.1, 0.9]
                return (mins, fractions)
        
        c = mycomposite()
        c.set_state(5e9,300)
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.7])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))
        c.set_state(5e9,3000)
        (m,f) = c.unroll()
        self.assertArraysAlmostEqual(f,[0.1,0.9])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))

    def test_number(self):
        class mycomposite(burnman.Material):
            def __init__(self):
                burnman.Material.__init__(self)

            def unroll(self):
                if (self.temperature>500):
                    return ([minerals.Murakami_etal_2012.fe_periclase()],[1.0])
                fractions = [0.3, 0.7]
                mins = [minerals.Murakami_etal_2012.fe_periclase(), minerals.SLB_2005.periclase()]
                return (mins, fractions)
        
        c = mycomposite()
        c.set_state(5e9,1000)
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[1.0])
        self.assertEqual(mins,minerals.Murakami_etal_2012.fe_periclase().to_string())
        c.set_state(5e9,300)
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.3,0.7])
        self.assertEqual(mins,",".join([minerals.Murakami_etal_2012.fe_periclase().to_string(),minerals.SLB_2005.periclase().to_string()]))

    def test_nest(self):
        min1 = minerals.Murakami_etal_2012.fe_periclase_LS()
        min2 = minerals.SLB_2005.periclase()
        ca = burnman.Composite( [min1], [1.0] )
        c = burnman.Composite( [ca, min2], [0.4, 0.6] )
        c.set_state(5e9,1000)
        (m,f) = c.unroll()
        mins=",".join([a.to_string() for a in m])
        self.assertArraysAlmostEqual(f,[0.4, 0.6])
        self.assertEqual(mins,",".join([min1.to_string(),min2.to_string()]))

    def test_density_composite(self):
        pyrolite = burnman.Composite([minerals.SLB_2005.mg_fe_perovskite(0.2), \
                                      minerals.SLB_2005.ferropericlase(0.4)], \
                                     [0.95, 0.05] )
        pyrolite.set_method('slb3')
        pyrolite.set_state(40.e9, 2000)
        
        d1 = int(pyrolite.phases[0].density)
        d2 = int(pyrolite.phases[1].density)
        dmix = int(pyrolite.density)
        assert(d1<dmix)
        assert(dmix<d2)
        assert(d1 == 4732)
        assert(d2 == 5275)
        assert(dmix == 4744)

    def test_thermodynamic_potentials(self):
        rock1 = burnman.Composite([minerals.SLB_2011.mg_perovskite(), \
                                   minerals.SLB_2011.periclase()], \
                                   [1.0, 0.0] )
        rock2 = burnman.Composite([minerals.SLB_2011.mg_perovskite(), \
                                   minerals.SLB_2011.periclase()], \
                                   [0.0, 1.0] )
        rock3 = burnman.Composite([minerals.SLB_2011.mg_perovskite(), \
                                   minerals.SLB_2011.periclase()], \
                                   [0.4, 0.6] )
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
        self.assertFloatEqual(G1,gibbs1)
        self.assertFloatEqual(G2,gibbs2)
        self.assertFloatEqual(G3,0.4*gibbs1+0.6*gibbs2)

        # Helmholtz
        helmholtz1 = mineral1.molar_helmholtz
        helmholtz2 = mineral2.molar_helmholtz
        F1 = rock1.molar_helmholtz
        F2 = rock2.molar_helmholtz
        F3 = rock3.molar_helmholtz
        self.assertFloatEqual(F1,helmholtz1)
        self.assertFloatEqual(F2,helmholtz2)
        self.assertFloatEqual(F3,0.4*helmholtz1+0.6*helmholtz2)

        # Enthalpy
        enthalpy1 = mineral1.molar_enthalpy
        enthalpy2 = mineral2.molar_enthalpy
        H1 = rock1.molar_enthalpy
        H2 = rock2.molar_enthalpy
        H3 = rock3.molar_enthalpy
        self.assertFloatEqual(H1,enthalpy1)
        self.assertFloatEqual(H2,enthalpy2)
        self.assertFloatEqual(H3,0.4*enthalpy1+0.6*enthalpy2)

        # Energy
        internal_energy1 = mineral1.internal_energy
        internal_energy2 = mineral2.internal_energy
        U1 = rock1.internal_energy
        U2 = rock2.internal_energy
        U3 = rock3.internal_energy
        self.assertFloatEqual(U1,internal_energy1)
        self.assertFloatEqual(U2,internal_energy2)
        self.assertFloatEqual(U3,0.4*internal_energy1+0.6*internal_energy2)

        # Entropy
        molar_entropy1 = mineral1.molar_entropy
        molar_entropy2 = mineral2.molar_entropy
        S1 = rock1.molar_entropy
        S2 = rock2.molar_entropy
        S3 = rock3.molar_entropy
        self.assertFloatEqual(S1,molar_entropy1)
        self.assertFloatEqual(S2,molar_entropy2)
        self.assertFloatEqual(S3,0.4*molar_entropy1+0.6*molar_entropy2)

    def test_summing_bigger(self):
        min1 = minerals.SLB_2005.periclase()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c = burnman.Composite( [min1, min1], [0.8, 0.4])
            assert len(w) == 1
        c.set_method("slb3")
        c.set_state(5e9,1000)
        (m,f) = c.unroll()
        self.assertArraysAlmostEqual(f, [2./3., 1./3.])

    def test_summing_smaller(self):
        min1 = minerals.SLB_2005.periclase()
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            c = burnman.Composite( [min1, min1], [0.4, 0.2])
            assert len(w) == 1
        c.set_method("slb3")
        c.set_state(5e9,1000)
        (m,f) = c.unroll()
        self.assertArraysAlmostEqual(f, [2./3., 1./3.])

    def test_summing_slightly_negative(self):
        min1 = minerals.SLB_2005.periclase()
        c = burnman.Composite( [min1, min1, min1], [0.8, 0.2, 1.0-0.8-0.2])
        c.set_method("slb3")
        c.set_state(5e9,1000)
        (m,f) = c.unroll()
        self.assertArraysAlmostEqual(f, [0.8, 0.2, 0.0])
        self.assertTrue(f[2]>=0.0)


if __name__ == '__main__':
    unittest.main()
