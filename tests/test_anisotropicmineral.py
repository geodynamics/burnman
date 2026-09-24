import unittest
from util import BurnManTest
import numpy as np

from burnman import AnisotropicMineral
from burnman.classes.anisotropicmineral import deformation_gradient_alpha_and_compliance
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.utils.anisotropy import (
    voigt_notation_to_stiffness_tensor,
    contract_stiffnesses,
)
from burnman.minerals.SLB_2011 import periclase, forsterite


def make_mineral(crystal_system="orthorhombic"):
    fo = forsterite()
    cell_lengths = np.array([4.7646, 10.2296, 5.9942])
    cell_lengths *= np.cbrt(fo.params["V_0"] / np.prod(cell_lengths))

    if crystal_system == "orthorhombic":
        # for orthorhombic,
        # 14, 15, 16, 24, 25, 26,
        # 34, 35, 36, 45, 46, 56 are all zero,
        # and alpha, beta, gamma are all 90 degrees.
        alpha, beta, gamma = [90.0, 90.0, 90.0]
        a14, a15, a16 = [0.0, 0.0, 0.0]
        a24, a25, a26 = [0.0, 0.0, 0.0]
        a34, a35, a36 = [0.0, 0.0, 0.0]
        a45, a46, a56 = [0.0, 0.0, 0.0]
    elif crystal_system == "monoclinic":
        # for monoclinic, 15, 25, 35, 46 are non-zero,
        # and beta is not 90 degrees.
        alpha, beta, gamma = [90.0, 100.0, 90.0]
        cell_lengths *= 1.0051159860149095
        a14, a15, a16 = [0.0, -1.0, 0.0]
        a24, a25, a26 = [0.0, 0.4, 0.0]
        a34, a35, a36 = [0.0, 0.2, 0.0]
        a45, a46, a56 = [0.0, 0.3, 0.0]
    elif crystal_system == "triclinic":
        # for triclinic, all elements can be non-zero,
        # and alpha, beta, gamma are not 90 degrees.
        alpha, beta, gamma = [85.0, 80.0, 87.0]
        cell_lengths *= 1.006635538793111
        a14, a15, a16 = [0.4, -1.0, -0.6]
        a24, a25, a26 = [0.0, 0.0, 0.0]
        a34, a35, a36 = [0.0, 0.0, 0.0]
        a45, a46, a56 = [1.97, 0.0, 0.0]
    else:
        raise ValueError(
            "crystal_system must be 'orthorhombic', " "'monoclinic' or 'triclinic'"
        )

    cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], alpha, beta, gamma]
    )

    constants = np.zeros((6, 6, 3, 1))
    constants[:, :, 1, 0] = np.array(
        [
            [0.44, -0.12, -0.1, a14, a15, a16],
            [-0.12, 0.78, -0.22, a24, a25, a26],
            [-0.1, -0.22, 0.66, a34, a35, a36],
            [a14, a24, a34, 1.97, a45, a46],
            [a15, a25, a35, a45, 1.61, a56],
            [a16, a26, a36, a46, a56, 1.55],
        ]
    )
    constants[:, :, 2, 0] = np.array(
        [
            [0.24, -0.12, -0.1, a14, a15, a16],
            [-0.12, 0.38, -0.22, a24, a25, a26],
            [-0.1, -0.22, 0.26, a34, a35, a36],
            [a14, a24, a34, 0.0, 0.0, 0.0],
            [a15, a25, a35, 0.0, 0.0, 0.0],
            [a16, a26, a36, 0.0, 0.0, 0.0],
        ]
    )

    m = AnisotropicMineral(fo, cell_parameters, constants, [0, 1, 2])
    return m


class test_anisotropic_mineral(BurnManTest):
    def test_matrix_exponential_derivatives(self):
        """
        Check Frechet derivatives against the closed form for diagonal PsiI.
        """
        eigenvalues = np.array([-0.1, -0.2, 0.15])
        PsiI = np.diag(eigenvalues)
        dPsiIdf = np.array([[0.2, 0.1, -0.04], [0.1, 0.3, 0.02], [-0.04, 0.02, 0.5]])
        dPsiIdT = 1.0e-5 * np.array(
            [[0.4, -0.2, 0.1], [-0.2, 0.2, 0.3], [0.1, 0.3, -0.6]]
        )
        alpha_V = 2.0e-5
        F, dFdf, alpha, _ = deformation_gradient_alpha_and_compliance(
            alpha_V, 1.0e-11, PsiI, np.eye(6), dPsiIdf, dPsiIdT
        )

        # L_exp(diag(lam), E)_ij = E_ij (exp(lam_i)-exp(lam_j))/(lam_i-lam_j).
        divided_difference = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                if i == j:
                    divided_difference[i, j] = np.exp(eigenvalues[i])
                else:
                    delta = eigenvalues[i] - eigenvalues[j]
                    divided_difference[i, j] = (
                        np.exp(eigenvalues[j]) * np.expm1(delta) / delta
                    )
        expected_F = np.diag(np.exp(eigenvalues))
        expected_dFdf = divided_difference * dPsiIdf
        expected_dFdT = divided_difference * (dPsiIdT + alpha_V * dPsiIdf)
        L = expected_dFdT @ np.diag(np.exp(-eigenvalues))
        self.assertArraysAlmostEqual(
            F.flatten(), expected_F.flatten(), tol=1.0e-12, tol_zero=1.0e-14
        )
        self.assertArraysAlmostEqual(
            dFdf.flatten(), expected_dFdf.flatten(), tol=1.0e-12, tol_zero=1.0e-14
        )
        self.assertArraysAlmostEqual(
            alpha.flatten(), (0.5 * (L + L.T)).flatten(), tol=1.0e-12, tol_zero=1.0e-18
        )

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
        self.assertArraysAlmostEqual(
            np.diag(per2.grueneisen_tensor), [gr, gr, gr], tol=1.0e-12
        )

    def test_orthotropic_consistency(self):
        m = make_mineral(crystal_system="orthorhombic")
        self.assertTrue(check_anisotropic_eos_consistency(m, tol=1.0e-6))

    def test_non_orthotropic_consistency(self):
        m = make_mineral(crystal_system="triclinic")
        self.assertTrue(check_anisotropic_eos_consistency(m, P=2.0e10, tol=1.0e-6))

    def test_stiffness(self):
        for crystal_system in ["orthorhombic", "monoclinic", "triclinic"]:
            m = make_mineral(crystal_system=crystal_system)
            m.set_state(1.0e9, 300.0)
            Cijkl = m.full_isothermal_stiffness_tensor
            Cij = m.isothermal_stiffness_tensor

            self.assertEqual(Cij[0, 0], Cijkl[0, 0, 0, 0])
            self.assertEqual(Cij[1, 1], Cijkl[1, 1, 1, 1])
            self.assertEqual(Cij[2, 2], Cijkl[2, 2, 2, 2])
            self.assertEqual(Cij[0, 1], Cijkl[0, 0, 1, 1])
            self.assertEqual(Cij[0, 2], Cijkl[0, 0, 2, 2])
            self.assertEqual(Cij[1, 2], Cijkl[1, 1, 2, 2])
            self.assertEqual(Cij[3, 3], Cijkl[1, 2, 1, 2])
            self.assertEqual(Cij[4, 4], Cijkl[0, 2, 0, 2])
            self.assertEqual(Cij[5, 5], Cijkl[0, 1, 0, 1])

    def test_stiffness_rotation(self):
        """
        Check that cell vectors, elastic stiffnesses and
        thermal expansivity share the crystal frame.
        """
        m = make_mineral(crystal_system="triclinic")
        m.set_state(1.0e9, 1000.0)
        Cij = m.isothermal_stiffness_tensor
        Cij_unrotated = m._unrotated_isothermal_stiffness_tensor

        # Obtain the transformation from the cell vectors, independently of the
        # rotation property. The [0, 1, 2] convention puts a along x and b in xy.
        cell_vectors = m.cell_vectors
        self.assertArraysAlmostEqual(
            cell_vectors[[0, 0, 1], [1, 2, 2]], np.zeros(3), tol_zero=1.0e-14
        )
        self.assertTrue(np.all(np.diag(cell_vectors) > 0.0))
        R = cell_vectors.T @ np.linalg.inv(m.unrotated_cell_vectors.T)
        self.assertArraysAlmostEqual(
            m.rotation_matrix.flatten(), R.flatten(), tol=1.0e-12, tol_zero=1.0e-14
        )

        # make sure R is not the identity matrix
        self.assertFalse(np.allclose(R, np.eye(3)))

        # convert the unrotated stiffness tensor to a full 4th order tensor
        Cijkl_unrotated = voigt_notation_to_stiffness_tensor(Cij_unrotated)

        # rotate the stiffness tensor
        Cijkl_rotated = np.einsum("ip,jq,kr,ls,pqrs->ijkl", R, R, R, R, Cijkl_unrotated)

        # convert back to Voigt notation
        Cij_rotated = contract_stiffnesses(Cijkl_rotated)

        # check that the rotated stiffness tensor is equal to that
        # obtained from the mineral object
        self.assertArraysAlmostEqual(
            Cij.flatten() / np.max(np.abs(Cij)),
            Cij_rotated.flatten() / np.max(np.abs(Cij)),
            tol=1.0e-12,
            tol_zero=1.0e-13,
        )
        self.assertArraysAlmostEqual(
            m.thermal_expansivity_tensor.flatten(),
            (R @ m._unrotated_alpha @ R.T).flatten(),
            tol=1.0e-12,
            tol_zero=1.0e-18,
        )

    def test_elastic_energy_rotation(self):
        """
        Check that the rotation is implemented correctly by ensuring that
        the quadratic elastic energy is unchanged from
        the unrotated to the rotated frame.
        """
        m = make_mineral(crystal_system="triclinic")
        m.set_state(20.0e9, 1000.0)
        strain = 1.0e-3 * np.array(
            [[1.0, 0.2, -0.3], [0.2, -0.4, 0.15], [-0.3, 0.15, -0.6]]
        )
        # Infer the rotation from the cell vectors,
        # not the rotation_matrix property.
        R = m.cell_vectors.T @ np.linalg.inv(m.unrotated_cell_vectors.T)
        strain_rotated = R @ strain @ R.T
        C = voigt_notation_to_stiffness_tensor(m._unrotated_isothermal_stiffness_tensor)
        energy = 0.5 * m.V * np.einsum("ij,ijkl,kl", strain, C, strain)
        energy_rotated = (
            0.5
            * m.V
            * np.einsum(
                "ij,ijkl,kl",
                strain_rotated,
                m.full_isothermal_stiffness_tensor,
                strain_rotated,
            )
        )
        self.assertFloatEqual(energy_rotated, energy, tol=1.0e-12)

    def test_monoclinic_consistency(self):
        """
        Test that the monoclinic mineral has the correct zero elements
        in the stiffness and compliance tensors
        after rotation to the crystallographic frame,
        and that the cell parameters are also consistent.
        """
        m = make_mineral(crystal_system="monoclinic")
        self.assertTrue(check_anisotropic_eos_consistency(m))

        m.set_state(20.0e9, 500.0)

        ST = m.isothermal_compliance_tensor
        SS = m.isentropic_compliance_tensor
        CT = m.isothermal_stiffness_tensor
        CS = m.isentropic_stiffness_tensor
        CT_unrotated = m._unrotated_isothermal_stiffness_tensor

        # 14, 16, 24, 26, 34, 36, 45, 56 should be equal to zero,
        # and alpha and gamma are 90 degrees.
        for i, j in [(0, 3), (0, 5), (1, 3), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]:
            self.assertFloatEqual(ST[i, j], 0.0)
            self.assertFloatEqual(SS[i, j], 0.0)
            self.assertFloatEqual(CT[i, j], 0.0)
            self.assertFloatEqual(CS[i, j], 0.0)
            self.assertFloatEqual(CT_unrotated[i, j], 0.0)

        p = m.cell_parameters
        self.assertFloatEqual(p[3], 90.0)
        self.assertFloatEqual(p[5], 90.0)


if __name__ == "__main__":
    unittest.main()
