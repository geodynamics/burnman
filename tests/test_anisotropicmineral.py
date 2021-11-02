from __future__ import absolute_import
import unittest
from util import BurnManTest
import numpy as np

from burnman import AnisotropicMineral
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.minerals.SLB_2011 import periclase, forsterite


class test_anisotropic_mineral(BurnManTest):
    def test_isotropic_grueneisen(self):

        per = periclase()
        a = np.cbrt(per.params['V_0'])
        beta_RT = 1./per.params['K_0']
        S44 = 1./per.params['G_0']
        S11 = (beta_RT + 3.*S44)/9.
        S12 = S11 - S44/2.

        cell_parameters = np.array([a, a, a, 90, 90, 90])

        constants = np.zeros((6, 6, 2, 1))
        constants[:3, :3, 1, 0] = S12 / beta_RT

        for i in range(3):
            constants[i, i, 1, 0] = S11 / beta_RT
            constants[i+3, i+3, 1, 0] = S44 / beta_RT

        per2 = AnisotropicMineral(per, cell_parameters, constants)

        P = 1.e9
        T = 1000.
        per.set_state(P, T)
        per2.set_state(P, T)

        gr = per.grueneisen_parameter
        self.assertArraysAlmostEqual(np.diag(per2.grueneisen_tensor),
                                     [gr, gr, gr])

    def test_orthorhombic_consistency(self):
        fo = forsterite()
        cell_lengths = np.array([4.7646, 10.2296, 5.9942])
        cell_lengths *= np.cbrt(fo.params['V_0'] / np.prod(cell_lengths))
        cell_parameters = np.array([cell_lengths[0],
                                    cell_lengths[1],
                                    cell_lengths[2],
                                    90, 90, 90])

        constants = np.zeros((6, 6, 2, 1))
        constants[:, :, 1, 0] = np.array([[0.44,  -0.12, -0.1,  0., 0., 0.],
                                          [-0.12,  0.78, -0.22, 0., 0., 0.],
                                          [-0.1,  -0.22,  0.66, 0., 0., 0.],
                                          [0.,     0.,    0., 1.97, 0., 0.],
                                          [0.,     0.,    0., 0., 1.61, 0.],
                                          [0.,     0.,    0., 0., 0.,   1.55]])

        m = AnisotropicMineral(fo, cell_parameters, constants)
        assert(check_anisotropic_eos_consistency(m))

    def test_triclinic_consistency(self):
        fo = forsterite()
        cell_lengths = np.array([4.7646, 10.2296, 5.9942])
        cell_lengths *= np.cbrt(fo.params['V_0'] / np.prod(cell_lengths))
        cell_lengths *= 1.0710965273664157
        cell_parameters = np.array([cell_lengths[0],
                                    cell_lengths[1],
                                    cell_lengths[2],
                                    60, 70, 80])

        constants = np.zeros((6, 6, 2, 1))
        constants[:, :, 1, 0] = np.array([[0.44,  -0.12, -0.1,  1., 0.5, 0.5],
                                          [-0.12,  0.78, -0.22, 0., 0., 0.],
                                          [-0.1,  -0.22,  0.66, 0., 0., 0.],
                                          [1.,     0.,    0., 1.97, 0., 0.],
                                          [0.5,    0.,    0., 0., 1.61, 0.],
                                          [0.5,    0.,    0., 0., 0.,   1.55]])

        m = AnisotropicMineral(fo, cell_parameters, constants)
        assert(check_anisotropic_eos_consistency(m))


if __name__ == '__main__':
    unittest.main()
