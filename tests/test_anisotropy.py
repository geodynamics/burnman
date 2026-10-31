import unittest
from util import BurnManTest
import numpy as np
from scipy.linalg import logm

from burnman.classes import anisotropy
from burnman.utils.unitcell import (
    cell_parameters_to_vectors,
    cell_vectors_to_parameters,
)
from burnman.utils.anisotropy import polar_decomposition


class test_anisotropy(BurnManTest):
    def test_invalid_elastic_inputs(self):
        cijs = anisotropy.IsotropicMaterial(
            3300.0, [40.0e9, 24.0e9]
        ).isentropic_stiffness_tensor
        for density in (0.0, -1.0, np.nan, np.inf):
            with self.subTest(density=density):
                with self.assertRaisesRegex(AssertionError, "density"):
                    anisotropy.AnisotropicMaterial(density, cijs)
        with self.assertRaisesRegex(AssertionError, "6x6"):
            anisotropy.AnisotropicMaterial(3300.0, np.eye(3))
        for value in (np.nan, np.inf):
            invalid = cijs.copy()
            invalid[0, 0] = value
            with self.subTest(stiffness=value):
                with self.assertRaisesRegex(AssertionError, "finite"):
                    anisotropy.AnisotropicMaterial(3300.0, invalid)

    def test_stiffness_requires_exact_input_symmetry(self):
        cijs = anisotropy.IsotropicMaterial(
            3300.0, [40.0e9, 24.0e9]
        ).isentropic_stiffness_tensor.copy()
        cijs[0, 1] = np.nextafter(cijs[0, 1], np.inf)
        self.assertTrue(np.allclose(cijs, cijs.T))
        original = cijs.copy()
        with self.assertRaisesRegex(AssertionError, "symmetric"):
            anisotropy.AnisotropicMaterial(3300.0, cijs)
        self.assertTrue(np.array_equal(cijs, original))

    def test_direction_validation(self):
        m = anisotropy.IsotropicMaterial(3300.0, [40.0e9, 24.0e9])
        invalid = (
            0.0,
            [],
            [0.0, 0.0, 0.0],
            [np.nan, 0.0, 1.0],
            [np.inf, 0.0, 1.0],
            [1.0, 0.0],
            [1.0 + 1.0j, 0.0, 0.0],
            [[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        )
        methods = (
            m.christoffel_tensor,
            m.wave_velocities,
            m.isentropic_linear_compressibility,
            m.isentropic_youngs_modulus,
        )
        for direction in invalid:
            for method in methods:
                with self.subTest(method=method.__name__, direction=direction):
                    with self.assertRaises(AssertionError):
                        method(direction)
            for method in (m.isentropic_shear_modulus, m.isentropic_poissons_ratio):
                for args in (
                    (direction, [0.0, 1.0, 0.0]),
                    ([1.0, 0.0, 0.0], direction),
                ):
                    with self.subTest(method=method.__name__, args=args):
                        with self.assertRaises(AssertionError):
                            method(*args)

    def test_batched_direction_normalization(self):
        lam, shear = 40.0e9, 24.0e9
        m = anisotropy.IsotropicMaterial(3300.0, [lam, shear])
        directions = np.array(
            [[[2.0, 0.0, 0.0], [0.0, 3.0, 0.0]], [[0.0, 0.0, 4.0], [3.0, 4.0, 0.0]]]
        )
        lateral = np.array(
            [[[0.0, 2.0, 0.0], [0.0, 0.0, 3.0]], [[4.0, 0.0, 0.0], [-4.0, 3.0, 0.0]]]
        )
        bulk = lam + 2.0 * shear / 3.0
        young = shear * (3.0 * lam + 2.0 * shear) / (lam + shear)
        poisson = lam / (2.0 * (lam + shear))
        cases = (
            (m.isentropic_linear_compressibility(directions), 1.0 / (3.0 * bulk)),
            (m.isentropic_youngs_modulus(directions), young),
            (m.isentropic_shear_modulus(directions, lateral), shear),
            (m.isentropic_poissons_ratio(directions, lateral), poisson),
        )
        for actual, expected in cases:
            self.assertEqual(actual.shape, (2, 2))
            self.assertArraysAlmostEqual(
                actual.flatten(), np.full(4, expected), tol=1.0e-12
            )

    def test_wave_polarizations_are_eigenvector_columns(self):
        m = anisotropy.OrthorhombicMaterial(
            3300.0,
            np.array([220.0, 80.0, 70.0, 240.0, 65.0, 260.0, 60.0, 68.0, 73.0]) * 1.0e9,
        )
        directions = np.array(
            [[[1.0, 2.0, 3.0], [2.0, -4.0, 1.0]], [[3.0, 0.0, 1.0], [-1.0, 1.0, 5.0]]]
        )
        for direction in (directions[0, 0], directions):
            with self.subTest(shape=direction.shape):
                acoustic = m.christoffel_tensor(direction)
                speeds, vectors = m.wave_velocities(direction)
                self.assertTrue(np.all(np.diff(speeds, axis=-1) <= 0.0))
                scale = np.max(np.abs(acoustic))
                actual = acoustic @ vectors
                expected = vectors * (m.density * speeds**2)[..., None, :]
                self.assertArraysAlmostEqual(
                    (actual / scale).flatten(),
                    (expected / scale).flatten(),
                    tol=1.0e-12,
                    tol_zero=1.0e-12,
                )
                gram = vectors.swapaxes(-1, -2) @ vectors
                self.assertArraysAlmostEqual(
                    gram.flatten(),
                    np.broadcast_to(np.eye(3), gram.shape).flatten(),
                    tol=1.0e-12,
                    tol_zero=1.0e-12,
                )

    def test_stability_assertion_is_opt_in(self):
        stable = anisotropy.IsotropicMaterial(3300.0, [40.0e9, 24.0e9])
        stable.assert_stable()
        unstable = anisotropy.IsotropicMaterial(3300.0, [40.0e9, -1.0e9])
        self.assertLess(unstable.isentropic_stiffness_tensor[3, 3], 0.0)
        with self.assertRaisesRegex(AssertionError, "positive definite"):
            unstable.assert_stable()
        with np.errstate(invalid="raise"):
            with self.assertRaisesRegex(AssertionError, "eigenvalues"):
                unstable.wave_velocities([1.0, 0.0, 0.0])

    def test_rhombohedral_rejects_inconsistent_c66(self):
        for cijs in (
            [220.0, 80.0, 70.0, 8.0, 250.0, 65.0, 75.0],
            [220.0, 80.0, 70.0, 8.0, 6.0, 250.0, 65.0, 75.0],
        ):
            with self.subTest(cijs=cijs):
                with self.assertRaisesRegex(AssertionError, "C66"):
                    anisotropy.RhombohedralMaterial(3300.0, np.array(cijs) * 1.0e9)

    def test_rhombohedral_threefold_symmetry(self):
        angle = 2.0 * np.pi / 3.0
        c, s = np.cos(angle), np.sin(angle)
        rotation = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        direction = np.array([1.0, 2.0, 3.0])
        for cijs in (
            [220.0, 80.0, 70.0, 8.0, 250.0, 65.0, 70.0],
            [220.0, 80.0, 70.0, 8.0, 6.0, 250.0, 65.0, 70.0],
        ):
            with self.subTest(cijs=cijs):
                m = anisotropy.RhombohedralMaterial(3300.0, np.array(cijs) * 1.0e9)
                stiffness = m.full_isentropic_stiffness_tensor
                rotated = np.einsum(
                    "ia,jb,kc,ld,abcd->ijkl",
                    rotation,
                    rotation,
                    rotation,
                    rotation,
                    stiffness,
                )
                scale = np.max(np.abs(stiffness))
                self.assertArraysAlmostEqual(
                    (rotated / scale).flatten(),
                    (stiffness / scale).flatten(),
                    tol=1.0e-12,
                    tol_zero=1.0e-12,
                )
                self.assertArraysAlmostEqual(
                    m.wave_velocities(direction)[0],
                    m.wave_velocities(rotation @ direction)[0],
                    tol=1.0e-12,
                )

    def test_isotropic(self):
        l, G = (0.4e11, 0.24e11)
        elastic_constants = [l, G]
        rho = 2.735e3
        m = anisotropy.IsotropicMaterial(rho, elastic_constants)

        E = G * (3.0 * l + 2.0 * G) / (l + G)
        Vp = np.sqrt(
            (
                m.isentropic_bulk_modulus_reuss
                + 4.0 / 3.0 * m.isentropic_shear_modulus_reuss
            )
            / rho
        )
        Vs = np.sqrt(m.isentropic_shear_modulus_reuss / rho)

        array_1 = [
            m.isentropic_bulk_modulus_reuss,
            m.isentropic_bulk_modulus_vrh,
            m.isentropic_bulk_modulus_voigt,
            m.isentropic_shear_modulus_reuss,
            m.isentropic_shear_modulus_vrh,
            m.isentropic_shear_modulus_voigt,
            m.isentropic_universal_elastic_anisotropy + 6.0,
            m.isentropic_isotropic_poisson_ratio,
            E,
            Vp,
            Vs,
            Vs,
        ]

        d1 = [1.0, 0.0, 0.0]
        d2 = [0.0, 1.0, 0.0]

        ds = np.array([d1, d2])
        beta_100 = m.isentropic_linear_compressibility(direction=ds)[0]
        E_100 = m.isentropic_youngs_modulus(direction=ds)[0]
        G_100_010 = m.isentropic_shear_modulus(plane_normal=d1, shear_direction=d2)
        nu_100_010 = m.isentropic_poissons_ratio(
            axial_direction=d1, lateral_direction=d2
        )
        wave_speeds, _ = m.wave_velocities(propagation_direction=d1)
        Vp, Vs1, Vs2 = wave_speeds

        array_2 = [
            1.0 / 3.0 / beta_100,
            1.0 / 3.0 / beta_100,
            1.0 / 3.0 / beta_100,
            G_100_010,
            G_100_010,
            G,
            6.0,
            nu_100_010,
            E_100,
            Vp,
            Vs1,
            Vs2,
        ]

        self.assertArraysAlmostEqual(array_1, array_2)

    def test_christoffel_under_hydrostatic_stress(self):
        """
        Hydrostatic stress preserves the isotropic acoustic tensor in every direction.
        """
        lam, G = 0.4e11, 0.24e11
        m = anisotropy.IsotropicMaterial(2735.0, [lam, G])
        directions = np.array(
            [[[1.0, 0.0, 0.0], [-1.0, 2.0, 3.0]], [[0.0, 2.0, 0.0], [3.0, 0.0, 4.0]]]
        )
        n = directions / np.linalg.norm(directions, axis=-1, keepdims=True)
        expected = G * np.eye(3) + (lam + G) * np.einsum("...i,...k->...ik", n, n)

        actual = m.christoffel_tensor(directions)
        self.assertEqual(actual.shape, expected.shape)
        self.assertArraysAlmostEqual(
            (actual / G).flatten(),
            (expected / G).flatten(),
            tol=1.0e-12,
            tol_zero=1.0e-12,
        )

    def test_crystal_systems(self):
        m = anisotropy.IsotropicMaterial(3000.0, [1.0, 1.0])
        m = anisotropy.CubicMaterial(3000.0, [1.0, 1.0, 1])
        m = anisotropy.HexagonalMaterial(3000.0, [1.0, 1.0, 1.0, 1.0, 1.0])
        m = anisotropy.TetragonalMaterial(3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        m = anisotropy.TetragonalMaterial(3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        m = anisotropy.RhombohedralMaterial(3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0])
        m = anisotropy.RhombohedralMaterial(
            3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
        )
        m = anisotropy.OrthorhombicMaterial(
            3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        m = anisotropy.MonoclinicMaterial(
            3000.0, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        m = anisotropy.TriclinicMaterial(
            3000.0,
            [
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
                1.0,
            ],
        )

    def test_tetragonal_ii(self):
        m = anisotropy.TetragonalMaterial(3000.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
        self.assertFloatEqual(
            m.isentropic_stiffness_tensor[0][5], -m.isentropic_stiffness_tensor[1][5]
        )

    def test_rhombohedral_i(self):
        m = anisotropy.RhombohedralMaterial(
            3000.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, -0.5]
        )
        self.assertFloatEqual(
            m.isentropic_stiffness_tensor[0][3], -m.isentropic_stiffness_tensor[1][3]
        )

    def test_rhombohedral_ii(self):
        m = anisotropy.RhombohedralMaterial(
            3000.0, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, -0.5]
        )
        self.assertArraysAlmostEqual(
            [m.isentropic_stiffness_tensor[0][3], m.isentropic_stiffness_tensor[0][4]],
            [
                -m.isentropic_stiffness_tensor[1][3],
                -m.isentropic_stiffness_tensor[1][4],
            ],
        )

    def test_cell_parameters_conversions(self):
        cell_parameters = np.array([1.0, 4.2, 2.2, 80.0, 85.0, 88.0])
        lengths_orig = cell_parameters[:3]
        cos_a = np.cos(np.deg2rad(cell_parameters[3:]))
        for convention in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 1, 0],
            [2, 0, 1],
        ]:
            v = cell_parameters_to_vectors(cell_parameters, convention)
            p = cell_vectors_to_parameters(v, convention)
            self.assertArraysAlmostEqual(cell_parameters, p)

            # check angles
            lengths = np.linalg.norm(v, axis=1)
            cos_a2 = np.array(
                [
                    np.dot(v[1], v[2]) / lengths[1] / lengths[2],
                    np.dot(v[0], v[2]) / lengths[0] / lengths[2],
                    np.dot(v[0], v[1]) / lengths[0] / lengths[1],
                ]
            )
            self.assertArraysAlmostEqual(cos_a, cos_a2)
            self.assertArraysAlmostEqual(lengths_orig, lengths)

    def test_cell_parameters_conversions_with_rotation(self):
        cell_parameters = np.array([1.0, 4.2, 2.2, 80.0, 85.0, 88.0])
        lengths_orig = cell_parameters[:3]
        cos_a = np.cos(np.deg2rad(cell_parameters[3:]))
        for convention in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 1, 0],
            [2, 0, 1],
        ]:
            v = cell_parameters_to_vectors(cell_parameters, convention)

            # make a rotation matrix
            theta = 0.34
            c, s = np.cos(theta), np.sin(theta)
            R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
            v = np.dot(v, R)
            p = cell_vectors_to_parameters(v, convention, rotate=True)
            self.assertArraysAlmostEqual(cell_parameters, p)

            # check angles
            lengths = np.linalg.norm(v, axis=1)
            cos_a2 = np.array(
                [
                    np.dot(v[1], v[2]) / lengths[1] / lengths[2],
                    np.dot(v[0], v[2]) / lengths[0] / lengths[2],
                    np.dot(v[0], v[1]) / lengths[0] / lengths[1],
                ]
            )
            self.assertArraysAlmostEqual(cos_a, cos_a2)
            self.assertArraysAlmostEqual(lengths_orig, lengths)

    def test_polar_decomposition(self):
        F = np.array([[1.0, 0.2, 0.3], [0.0, 1.0, 0.4], [0.0, 0.0, 1.0]])
        R1, U = polar_decomposition(F)
        R2, U = polar_decomposition(F, "U")
        R3, V = polar_decomposition(F, "V")
        R4, H = polar_decomposition(F, "H")

        # check that R is a proper rotation matrix
        self.assertArraysAlmostEqual(
            np.dot(R1.T, R1).flatten(), np.eye(3).flatten(), tol_zero=1e-12
        )

        # check that the rotation matrices are the same
        self.assertArraysAlmostEqual(R1.flatten(), R2.flatten(), tol_zero=1e-12)
        self.assertArraysAlmostEqual(R1.flatten(), R3.flatten(), tol_zero=1e-12)
        self.assertArraysAlmostEqual(R1.flatten(), R4.flatten(), tol_zero=1e-12)

        # check that F = R * U = V * R
        self.assertArraysAlmostEqual(
            np.dot(R1, U).flatten(), np.dot(V, R1).flatten(), tol_zero=1e-12
        )

        # check that H = logm(U)
        self.assertArraysAlmostEqual(logm(U).flatten(), H.flatten(), tol_zero=1e-12)


if __name__ == "__main__":
    unittest.main()
