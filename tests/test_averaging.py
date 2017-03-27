from __future__ import absolute_import
import unittest
import os
import sys
import warnings
sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals
import burnman.averaging_schemes as avg
from util import BurnManTest


class mypericlase (burnman.Mineral):

    """
    Stixrude & Lithgow-Bertelloni 2005 and references therein
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8}
        burnman.Mineral.__init__(self)


class my_nonrigid_mineral (burnman.Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 0.e9,
            'Gprime_0': 0.0,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 0.0}
        burnman.Mineral.__init__(self)


class VRH_average(BurnManTest):

    def test_one_object(self):
        v = avg.voigt_reuss_hill_function([2.0], [0.123])
        self.assertFloatEqual(0.123, v)

    def test_two_same(self):
        v = avg.voigt_reuss_hill_function([2.0, 2.0], [0.456, 0.456])
        self.assertFloatEqual(0.456, v)

    def test_one_no_volume(self):
        v = avg.voigt_reuss_hill_function([0.0, 2.0], [0.123, 0.456])
        self.assertFloatEqual(0.456, v)

    def test_mix(self):
        v = avg.voigt_reuss_hill_function([1.0, 2.0], [0.1, 0.2])
        self.assertFloatEqual(0.15833333333333, v)


class VRH(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_method('slb3')
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K_vrh[0] / 1.e9)
        self.assertFloatEqual(150.901, G_vrh[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number, [1.0 / number] * number)

        rock.set_method('slb3')
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K_vrh[0] / 1.e9)
        self.assertFloatEqual(150.901, G_vrh[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_method('slb3')
        rock.set_averaging_scheme(avg.VoigtReussHill())
        rho, v_p, v_s, v_phi, K_vrh, G_vrh = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(9957.665, v_p[0])
        self.assertFloatEqual(5606.916, v_s[0])
        self.assertFloatEqual(7565.607, v_phi[0])
        self.assertFloatEqual(279.408, K_vrh[0] / 1.e9)
        self.assertFloatEqual(153.461, G_vrh[0] / 1.e9)


class Reuss(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number,  [1.0 / number] * number)
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.Reuss())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(9887.625, v_p[0])
        self.assertFloatEqual(5606.745, v_s[0])
        self.assertFloatEqual(7473.353, v_phi[0])
        self.assertFloatEqual(272.635, K[0] / 1.e9)
        self.assertFloatEqual(153.452, G[0] / 1.e9)

    def test_non_present_non_rigid_phase(self):
        rock = burnman.Composite(
            [mypericlase(), my_nonrigid_mineral()], [1.0, 0.0])
        rock.set_averaging_scheme(avg.Reuss())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rho, v_p, v_s, v_phi, K, G = \
                                         rock.evaluate(['rho', 'v_p', 'v_s',
                                                        'v_phi', 'K_S', 'G'],
                                                       [10.e9], [300.])
            assert len(w) == 1 # we expect one error to be thrown when p_nonrigid = 0
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_present_non_rigid_phase(self):
        rock = burnman.Composite(
            [mypericlase(), my_nonrigid_mineral()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.Reuss())
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rho, v_p, v_s, v_phi, K, G = \
                                         rock.evaluate(['rho', 'v_p', 'v_s',
                                                        'v_phi', 'K_S', 'G'],
                                                       [10.e9], [300.])
            assert len(w) == 2 # we expect two errors to be thrown when p_nonrigid != 0
        self.assertFloatEqual(0.0, G[0] / 1.e9)


class Voigt(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.Voigt())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(10027.216, v_p[0])
        self.assertFloatEqual(5607.087, v_s[0])
        self.assertFloatEqual(7656.750, v_phi[0])
        self.assertFloatEqual(286.180, K[0] / 1.e9)
        self.assertFloatEqual(153.471, G[0] / 1.e9)


class HSLower(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.HashinShtrikmanLower())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(9951.957, v_p[0])
        self.assertFloatEqual(5606.916, v_s[0])
        self.assertFloatEqual(7558.094, v_phi[0])
        self.assertFloatEqual(278.853, K[0] / 1.e9)
        self.assertFloatEqual(153.461, G[0] / 1.e9)


class HSUpper(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.HashinShtrikmanUpper())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(9952.798, v_p[0])
        self.assertFloatEqual(5606.925, v_s[0])
        self.assertFloatEqual(7559.192, v_phi[0])
        self.assertFloatEqual(278.934, K[0] / 1.e9)
        self.assertFloatEqual(153.462, G[0] / 1.e9)


class HSAverage(BurnManTest):

    def test_1(self):
        rock = burnman.Composite([mypericlase()], [1.0])
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def same(self, number):
        rock = burnman.Composite(
            [mypericlase()] * number, [1.0 / number] * number)
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(3791.392, rho[0])
        self.assertFloatEqual(10285.368, v_p[0])
        self.assertFloatEqual(6308.811, v_s[0])
        self.assertFloatEqual(7260.900, v_phi[0])
        self.assertFloatEqual(199.884, K[0] / 1.e9)
        self.assertFloatEqual(150.901, G[0] / 1.e9)

    def test_same(self):
        self.same(2)
        self.same(3)
        self.same(4)

    def test_two_different(self):
        rock = burnman.Composite(
            [minerals.SLB_2005.periclase(), minerals.SLB_2005.fe_perovskite()], [0.5, 0.5])
        rock.set_averaging_scheme(avg.HashinShtrikmanAverage())
        rho, v_p, v_s, v_phi, K, G = \
            rock.evaluate(
                ['rho', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], [10.e9], [300.])
        self.assertFloatEqual(4881.469, rho[0])
        self.assertFloatEqual(9952.378, v_p[0])
        self.assertFloatEqual(5606.920, v_s[0])
        self.assertFloatEqual(7558.643, v_phi[0])
        self.assertFloatEqual(278.893, K[0] / 1.e9)
        self.assertFloatEqual(153.461, G[0] / 1.e9)

if __name__ == '__main__':
    unittest.main()
