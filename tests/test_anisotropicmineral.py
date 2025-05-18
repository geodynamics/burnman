import unittest
from util import BurnManTest
import numpy as np

from burnman import AnisotropicMineral
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.minerals.SLB_2011 import periclase, forsterite


def make_forsterite(orthotropic=True):
    fo = forsterite()
    cell_lengths = np.array([4.7646, 10.2296, 5.9942])
    cell_lengths *= np.cbrt(fo.params["V_0"] / np.prod(cell_lengths))

    if orthotropic:
        alpha, beta, gamma, a, b, c = [90.0, 90.0, 90.0, 0.0, 0.0, 0.0]
    else:
        alpha, beta, gamma, a, b, c = [85.0, 80.0, 87.0, 0.4, -1.0, -0.6]
        cell_lengths *= 1.006635538793111

    cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], alpha, beta, gamma]
    )

    constants = np.zeros((6, 6, 3, 1))
    constants[:, :, 1, 0] = np.array(
        [
            [0.44, -0.12, -0.1, a, b, c],
            [-0.12, 0.78, -0.22, 0.0, 0.0, 0.0],
            [-0.1, -0.22, 0.66, 0.0, 0.0, 0.0],
            [a, 0.0, 0.0, 1.97, 0.0, 0.0],
            [b, 0.0, 0.0, 0.0, 1.61, 0.0],
            [c, 0.0, 0.0, 0.0, 0.0, 1.55],
        ]
    )
    constants[:, :, 2, 0] = np.array(
        [
            [0.24, -0.12, -0.1, 0.0, 0.0, 0.0],
            [-0.12, 0.38, -0.22, a, a, a],
            [-0.1, -0.22, 0.26, c, b, a],
            [0.0, a, c, 0.0, 0.0, 0.0],
            [0.0, a, b, 0.0, 0.0, 0.0],
            [0.0, a, a, 0.0, 0.0, 0.0],
        ]
    )

    m = AnisotropicMineral(fo, cell_parameters, constants, [0, 1, 2])
    return m


class test_anisotropic_mineral(BurnManTest):
    def test_isotropic_grueneisen(self):
        per = periclase()
        a = np.cbrt(per.params["V_0"])
        beta_RT = 1.0 / per.params["K_0"]
        S44 = 1.0 / per.params["G_0"]
        S11 = (beta_RT + 3.0 * S44) / 9.0
        S12 = S11 - S44 / 2.0

        cell_parameters = np.array([a, a, a, 90, 90, 90])

        constants = np.zeros((6, 6, 2, 1))
        constants[:3, :3, 1, 0] = S12 / beta_RT

        for i in range(3):
            constants[i, i, 1, 0] = S11 / beta_RT
            constants[i + 3, i + 3, 1, 0] = S44 / beta_RT

        per2 = AnisotropicMineral(per, cell_parameters, constants)

        P = 1.0e9
        T = 1000.0
        per.set_state(P, T)
        per2.set_state(P, T)

        gr = per.grueneisen_parameter
        self.assertArraysAlmostEqual(np.diag(per2.grueneisen_tensor), [gr, gr, gr])

    def test_orthotropic_consistency(self):
        m = make_forsterite(orthotropic=True)
        self.assertTrue(check_anisotropic_eos_consistency(m))

    def test_non_orthotropic_consistency(self):
        m = make_forsterite(orthotropic=False)
        self.assertTrue(check_anisotropic_eos_consistency(m, P=2.0e10, tol=1.0e-6))

    def test_stiffness(self):
        for orthotropic in [True, False]:
            m = make_forsterite(orthotropic)
            m.set_state(1.0e9, 300.0)
            Cijkl = m.full_isothermal_stiffness_tensor
            Cij = m.isothermal_stiffness_tensor

            self.assertFloatEqual(Cij[0, 0], Cijkl[0, 0, 0, 0])
            self.assertFloatEqual(Cij[1, 1], Cijkl[1, 1, 1, 1])
            self.assertFloatEqual(Cij[2, 2], Cijkl[2, 2, 2, 2])
            self.assertFloatEqual(Cij[0, 1], Cijkl[0, 0, 1, 1])
            self.assertFloatEqual(Cij[0, 2], Cijkl[0, 0, 2, 2])
            self.assertFloatEqual(Cij[1, 2], Cijkl[1, 1, 2, 2])
            self.assertFloatEqual(Cij[3, 3], Cijkl[1, 2, 1, 2])
            self.assertFloatEqual(Cij[4, 4], Cijkl[0, 2, 0, 2])
            self.assertFloatEqual(Cij[5, 5], Cijkl[0, 1, 0, 1])


if __name__ == "__main__":
    unittest.main()
