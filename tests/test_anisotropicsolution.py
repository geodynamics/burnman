import unittest
from util import BurnManTest
import numpy as np

from burnman.constants import Avogadro
from burnman import AnisotropicMineral, AnisotropicSolution
from burnman import RelaxedAnisotropicSolution
from burnman.classes.solutionmodel import SymmetricRegularSolution
from burnman.tools.eos import check_eos_consistency
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.minerals.SLB_2011 import forsterite
from burnman.utils.unitcell import cell_parameters_to_vectors


def make_nonorthotropic_mineral(a, b, c, alpha, beta, gamma, d, e, f):
    fo = forsterite()
    cell_lengths = np.array([a, b, c])
    cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], alpha, beta, gamma]
    )
    frame_convention = [0, 1, 2]
    fo.params["V_0"] = np.linalg.det(
        cell_parameters_to_vectors(cell_parameters, frame_convention)
    )
    constants = np.zeros((6, 6, 3, 1))
    constants[:, :, 1, 0] = np.array(
        [
            [0.44, -0.12, -0.1, 0.0, d, 0.0],
            [-0.12, 0.78, -0.22, 0.0, e, 0.0],
            [-0.1, -0.22, 0.66, 0.0, f, 0.0],
            [0.0, 0.0, 0.0, 1.97, 0.0, 0.0],
            [d, e, f, 0.0, 1.61, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.55],
        ]
    )
    constants[:, :, 2, 0] = np.array(
        [
            [0.24, -0.12, -0.1, 0.0, d, 0.0],
            [-0.12, 0.38, -0.22, 0.0, d, 0.0],
            [-0.1, -0.22, 0.26, 0.0, e, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, f],
            [d, d, e, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, f, 0.0, 0.0],
        ]
    )

    m = AnisotropicMineral(fo, cell_parameters, constants, frame_convention)
    return m


def make_nonorthotropic_solution(two_fos=False, monoclinic=False):

    cell_lengths_A = np.array([4.7646, 10.2296, 5.9942])
    lth = cell_lengths_A * 1.0e-10 * np.cbrt(Avogadro / 4.0)

    if not monoclinic:
        m1 = make_nonorthotropic_mineral(
            lth[0], lth[1], lth[2], 85.0, 80.0, 87.0, 0.4, -1.0, -0.6
        )
        m2 = make_nonorthotropic_mineral(
            lth[0] * 1.2, lth[1] * 1.4, lth[2], 90.0, 90.0, 90.0, 0.4, -1.0, -0.6
        )
    else:
        m1 = make_nonorthotropic_mineral(
            lth[0], lth[1], lth[2], 90.0, 100.0, 90.0, 0.1, 0.2, 0.4
        )
        m2 = make_nonorthotropic_mineral(
            lth[0] * 1.2, lth[1] * 1.4, lth[2], 90.0, 100.0, 90.0, 0.9, 0.2, 0.1
        )

    if two_fos:
        n_mbrs = 3
        solution_model = SymmetricRegularSolution(
            endmembers=[[m1, "[Mg]2SiO4"], [m1, "[Mg]2SiO4"], [m2, "[Fe]2SiO4"]],
            energy_interaction=[[-1.0e3, 10.0e3], [10.0e3]],
            volume_interaction=[[0.0, 1.0e-6], [1.0e-6]],
        )
    else:
        n_mbrs = 2
        solution_model = SymmetricRegularSolution(
            endmembers=[[m1, "[Mg]2SiO4"], [m2, "[Fe]2SiO4"]],
            energy_interaction=[[10.0e3]],
            volume_interaction=[[1.0e-6]],
        )

    def fn(lnV, Pth, X, params):
        z = np.zeros((6, 6))
        f = np.zeros((3, 3, n_mbrs))
        return (z, z, z, f)

    prm = {}

    ss = AnisotropicSolution(
        name="double forsterite",
        solution_model=solution_model,
        psi_excess_function=fn,
        anisotropic_parameters=prm,
    )

    return ss


class test_two_member_solution(BurnManTest):
    def test_composition_derivatives_at_fixed_volume(self):
        """Compare cell and deformation derivatives with composition changes."""
        ss = make_nonorthotropic_solution()
        fractions = np.array([0.8, 0.2])
        direction = np.array([1.0, -1.0])
        ss.set_composition(fractions)
        ss.set_state(2.0e10, 1000.0)
        volume, temperature = ss.V, ss.T
        dM0 = np.einsum("ijk,k->ij", ss._dM0dp_fixed_VT, direction)
        dF = np.einsum("ijk,k->ij", ss._unrotated_dFdp_fixed_VT, direction)

        # Keep a large enough reference step to avoid cancellation in M0.
        dp = 1.0e-4
        try:
            ss.set_composition(fractions + dp * direction)
            ss.set_state_with_volume(volume, temperature)
            M0_plus = ss.cell_vectors_0.T.copy()
            F_plus = ss._unrotated_F.copy()
            ss.set_composition(fractions - dp * direction)
            ss.set_state_with_volume(volume, temperature)
            M0_minus = ss.cell_vectors_0.T.copy()
            F_minus = ss._unrotated_F.copy()
            self.assertArraysAlmostEqual(
                dM0.flatten(),
                ((M0_plus - M0_minus) / (2.0 * dp)).flatten(),
                tol=1.0e-8,
                tol_zero=1.0e-14,
            )
            self.assertArraysAlmostEqual(
                dF.flatten(),
                ((F_plus - F_minus) / (2.0 * dp)).flatten(),
                tol=1.0e-8,
                tol_zero=1.0e-12,
            )
        finally:
            ss.set_composition(fractions)
            ss.set_state_with_volume(volume, temperature)

    def test_volume(self):
        ss = make_nonorthotropic_solution()
        ps = [0.2, 0.8]
        ss.set_composition(ps)
        Vs = np.array([np.linalg.det(mbr[0].cell_vectors_0) for mbr in ss.endmembers])
        V1 = np.power(Vs[0], ps[0]) * np.power(Vs[1], ps[1])
        V2 = np.linalg.det(ss.cell_vectors_0)
        self.assertFloatEqual(V1, V2, tol=1.0e-12, tol_zero=0.0)

    def test_non_orthotropic_endmember_consistency(self):
        ss = make_nonorthotropic_solution()
        ss.set_composition([1.0, 0.0])
        self.assertTrue(check_anisotropic_eos_consistency(ss, P=2.0e10, tol=1.0e-6))

    def test_non_orthotropic_solution_consistency(self):
        ss = make_nonorthotropic_solution()
        ss.set_composition([0.8, 0.2])
        self.assertTrue(
            check_anisotropic_eos_consistency(ss, P=2.0e10, T=1000.0, tol=2.0e-7)
        )

    def test_relaxed_non_orthotropic_solution_consistency(self):
        ss = make_nonorthotropic_solution(two_fos=True)
        ss = RelaxedAnisotropicSolution(
            ss, [[1.0, -1.0, 0.0]], [[0.5, 0.5, 0.0], [0.0, 0.0, 1.0]]
        )
        ss.set_composition([0.8, 0.2], relaxed=False)
        self.assertTrue(
            check_anisotropic_eos_consistency(ss, P=2.0e10, T=1000.0, tol=2.0e-7)
        )

    def test_non_orthotropic_solution_clone(self):
        ss = make_nonorthotropic_solution()
        ss1 = ss._scalar_solution
        ss1.set_composition([0.8, 0.2])
        self.assertTrue(check_eos_consistency(ss1, P=2.0e10, T=1000.0, tol=1.0e-6))

        ss.set_composition([0.8, 0.2])
        ss.set_state(2.0e10, 1000.0)
        self.assertFloatEqual(
            ss.isothermal_bulk_modulus_reuss,
            ss1.isothermal_bulk_modulus_reuss,
            tol=1.0e-12,
        )

    def test_stiffness(self):
        ss = make_nonorthotropic_solution()
        ss.set_composition([0.8, 0.2])
        ss.set_state(1.0e9, 300.0)
        Cijkl = ss.full_isothermal_stiffness_tensor
        Cij = ss.isothermal_stiffness_tensor

        self.assertEqual(Cij[0, 0], Cijkl[0, 0, 0, 0])
        self.assertEqual(Cij[1, 1], Cijkl[1, 1, 1, 1])
        self.assertEqual(Cij[2, 2], Cijkl[2, 2, 2, 2])
        self.assertEqual(Cij[0, 1], Cijkl[0, 0, 1, 1])
        self.assertEqual(Cij[0, 2], Cijkl[0, 0, 2, 2])
        self.assertEqual(Cij[1, 2], Cijkl[1, 1, 2, 2])
        self.assertEqual(Cij[3, 3], Cijkl[1, 2, 1, 2])
        self.assertEqual(Cij[4, 4], Cijkl[0, 2, 0, 2])
        self.assertEqual(Cij[5, 5], Cijkl[0, 1, 0, 1])

    def test_monoclinic_solution_consistency(self):
        ss = make_nonorthotropic_solution(two_fos=False, monoclinic=True)
        ss.set_composition([0.3, 0.7])
        self.assertTrue(
            check_anisotropic_eos_consistency(ss, P=2.0e10, T=1000.0, tol=1.0e-6)
        )

        ss.set_state(20.0e9, 1000.0)

        ST = ss.isothermal_compliance_tensor
        SS = ss.isentropic_compliance_tensor
        CT = ss.isothermal_stiffness_tensor
        CS = ss.isentropic_stiffness_tensor
        CT_unrotated = ss._unrotated_isothermal_stiffness_tensor

        # 14, 16, 24, 26, 34, 36, 45, 56 should be equal to zero,
        # and alpha and gamma are 90 degrees.
        for i, j in [(0, 3), (0, 5), (1, 3), (1, 5), (2, 3), (2, 5), (3, 4), (4, 5)]:
            self.assertFloatEqual(ST[i, j], 0.0)
            self.assertFloatEqual(SS[i, j], 0.0)
            self.assertFloatEqual(CT[i, j], 0.0)
            self.assertFloatEqual(CS[i, j], 0.0)
            self.assertFloatEqual(CT_unrotated[i, j], 0.0)


if __name__ == "__main__":
    unittest.main()
