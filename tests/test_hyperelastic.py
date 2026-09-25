import unittest
from util import BurnManTest
import numpy as np
from burnman.minerals import SLB_2011
from burnman import AnisotropicMineral
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.tools.eos import check_hyperelastic_eos_consistency
from burnman import HyperelasticMineral
from scipy.linalg import expm
from scipy.optimize import brentq
from burnman.utils.anisotropy import (
    expand_strains,
    contract_stiffnesses,
    contract_compliances,
    voigt_notation_to_stiffness_tensor,
    voigt_notation_to_compliance_tensor,
)


class test_hyperelastic(BurnManTest):
    def _hydrostatic_property_in_current_frame(self, m, prop):
        """
        Express a property from the hydrostatic AnisotropicMineral
        instance in the current output frame of the HyperelasticMineral
        instance, rotating tensors as necessary.
        """
        value = getattr(m.hydrostatic, prop)

        # If the property is a scalar, just return the value.
        if np.ndim(value) == 0:
            return value

        # Otherwise, rotate the value into the current frame.
        Q = m.rotation_matrix @ m.hydrostatic.rotation_matrix.T

        # 2nd order tensors
        if value.shape == (3, 3):
            return Q @ value @ Q.T

        # Contracted 4th order tensors in Voigt notation
        if value.shape == (6, 6):
            if "compliance" in prop:
                value = voigt_notation_to_compliance_tensor(value)
                contract = contract_compliances
            else:
                value = voigt_notation_to_stiffness_tensor(value)
                contract = contract_stiffnesses
            return contract(np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, value))

        # 4th order tensors
        return np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, value)

    def _get_hyperelastic_mineral(self, orthotropic=False, thermal_shear=False):
        """
        Build a HyperelasticMineral instance with the requested
        symmetry and, optionally, a thermal shear component to the
        anisotropic properties.
        """
        fo = SLB_2011.forsterite()

        cell_lengths = np.array([4.7646, 10.2296, 5.9942])
        cell_lengths *= np.cbrt(fo.params["V_0"] / np.prod(cell_lengths))
        if orthotropic:
            angles = [90.0, 90.0, 90.0]
        else:
            angles = [85.0, 87.0, 92.0]
            cell_lengths *= 1.0019924088272532
        cell_parameters = np.concatenate((cell_lengths, angles))

        constants = np.zeros((6, 6, 2, 2 if thermal_shear else 1))
        constants[:, :, 1, 0] = np.array(
            [
                [0.44, -0.12, -0.10, 0.10, 0.20, 0.30],
                [-0.12, 0.78, -0.22, 0.00, 0.30, 0.00],
                [-0.10, -0.22, 0.66, 0.70, 0.00, 0.40],
                [0.10, 0.00, 0.70, 1.97, 0.10, 0.40],
                [0.20, 0.30, 0.00, 0.10, 1.61, 0.20],
                [0.30, 0.00, 0.40, 0.40, 0.20, 1.55],
            ]
        )

        if orthotropic:
            constants[:3, 3:, 1, 0] = 0.0
            constants[3:, :3, 1, 0] = 0.0
            shear = constants[3:, 3:, 1, 0]
            constants[3:, 3:, 1, 0] = np.diag(np.diag(shear))

        if thermal_shear:
            # A thermal shear changes the polar rotation as T changes at
            # fixed cell vectors.
            constants[0, 3, 0, 1] = constants[3, 0, 0, 1] = 1.0e-12

        m_hydrostatic = AnisotropicMineral(
            fo,
            cell_parameters,
            constants,
            frame_convention=[0, 1, 2],
        )

        self.assertTrue(check_anisotropic_eos_consistency(m_hydrostatic))

        return HyperelasticMineral(m_hydrostatic)

    def _set_finite_strain_state(self, m):
        """
        Apply a fixed Hencky strain to the hydrostatic 20 GPa, 300 K cell
        and return that Hencky strain.
        """
        H = np.array(
            [
                [0.010, 0.004, -0.002],
                [0.004, -0.006, 0.003],
                [-0.002, 0.003, -0.004],
            ]
        )
        U = expm(H)

        P = 20.0e9
        T = 300.0  # K

        m.hydrostatic.set_state(P, T)
        cell_vectors = (U @ m.hydrostatic.unrotated_cell_vectors.T).T
        m.set_state_with_molar_cell_vectors(cell_vectors, T)
        return H

    def test_hydrostatic_geometry(self):
        """
        Check cell geometry and rotation when no extra strain is applied.
        """
        for orthotropic in [True, False]:
            with self.subTest(orthotropic=orthotropic):
                m = self._get_hyperelastic_mineral(orthotropic=orthotropic)
                m.hydrostatic.set_state(20.0e9, 300.0)
                cell_vectors = m.hydrostatic.cell_vectors.copy()
                m.set_state_with_molar_cell_vectors(cell_vectors, 300.0)

                for prop, shape in [
                    ("deformation_gradient_tensor", (3, 3)),
                    ("unrotated_cell_vectors", (3, 3)),
                    ("rotation_matrix", (3, 3)),
                    ("cell_vectors", (3, 3)),
                    ("cell_parameters", (6,)),
                ]:
                    with self.subTest(property=prop):
                        value = getattr(m, prop)
                        self.assertEqual(value.shape, shape)
                        self.assertTrue(np.all(np.isfinite(value)))

                self.assertArraysAlmostEqual(
                    m.hencky_strain.flatten(), np.zeros(9), tol_zero=1.0e-12
                )
                self.assertArraysAlmostEqual(
                    m.rotation_matrix.flatten(),
                    m.hydrostatic.rotation_matrix.flatten(),
                    tol_zero=1.0e-12,
                )
                self.assertArraysAlmostEqual(
                    m.unrotated_cell_vectors.flatten(),
                    m.hydrostatic.unrotated_cell_vectors.flatten(),
                    tol=1.0e-12,
                )
                self.assertArraysAlmostEqual(
                    (m.rotation_matrix @ m.unrotated_cell_vectors.T).flatten(),
                    cell_vectors.T.flatten(),
                    tol=1.0e-12,
                    tol_zero=1.0e-14,
                )

    def test_hydrostatic_limit(self):
        """
        Recover the hydrostatic EoS for each elastic symmetry and cell shape.
        """
        cases = [
            ("orthotropic", {"orthotropic": True}),
            ("nonorthotropic, oblique cell", {}),
        ]
        properties = [
            "helmholtz",
            "internal_energy",
            "entropy",
            "molar_isometric_heat_capacity",
            "molar_isostress_heat_capacity",
            "molar_heat_capacity_v",
            "grueneisen_parameter",
            "isothermal_bulk_modulus_voigt",
            "isothermal_compressibility_reuss",
            "isothermal_compressibility_voigt",
            "isentropic_bulk_modulus_voigt",
            "isentropic_compressibility_reuss",
            "isentropic_compressibility_voigt",
            "isothermal_compliance_tensor",
            "thermal_expansivity_tensor",
            "isothermal_stiffness_tensor",
            "full_isothermal_compliance_tensor",
            "full_isothermal_stiffness_tensor",
            "full_isentropic_compliance_tensor",
            "isentropic_compliance_tensor",
            "isentropic_stiffness_tensor",
            "full_isentropic_stiffness_tensor",
            "grueneisen_tensor",
            "isothermal_compressibility_tensor",
            "isentropic_compressibility_tensor",
            "thermal_stress_tensor",
        ]
        reference_names = {"molar_isostress_heat_capacity": "molar_heat_capacity_p"}
        for name, options in cases:
            m = self._get_hyperelastic_mineral(**options)
            m.hydrostatic.set_state(20.0e9, 300.0)
            m.set_state_with_molar_cell_vectors(
                m.hydrostatic.cell_vectors.copy(), 300.0
            )
            # Snapshot the reference before numerical derivatives perturb it.
            expected_properties = {
                p: self._hydrostatic_property_in_current_frame(
                    m, reference_names.get(p, p)
                )
                for p in properties
            }
            pressure = m.hydrostatic.pressure
            with self.subTest(model=name, property="cauchy_stress"):
                self.assertArraysAlmostEqual(
                    (m.cauchy_stress / 1.0e9).flatten(),
                    (-pressure * np.eye(3) / 1.0e9).flatten(),
                    tol_zero=1.0e-13 if options.get("orthotropic") else 1.0e-4,
                )
            for prop, expected in expected_properties.items():
                with self.subTest(model=name, property=prop):
                    if np.ndim(expected) == 0:
                        self.assertFloatEqual(getattr(m, prop), expected, tol=1.0e-6)
                        continue
                    scale = np.max(np.abs(expected))
                    self.assertGreater(scale, 0.0)
                    self.assertArraysAlmostEqual(
                        np.asarray(getattr(m, prop)).flatten() / scale,
                        expected.flatten() / scale,
                        tol=1.0e-6,
                        # Allow 10 Pa residuals in stiffness components
                        # that vanish by symmetry.
                        tol_zero=10.0 / scale if "stiffness" in prop else 1.0e-7,
                    )

            directions = np.array(
                [
                    [[1.0, 0.0, 0.0], [1.0, 2.0, 3.0]],
                    [[0.0, 2.0, 0.0], [-2.0, 1.0, 4.0]],
                ]
            )
            expected = np.array(
                [m.hydrostatic.christoffel_tensor(n) for n in directions.reshape(-1, 3)]
            ).reshape(2, 2, 3, 3)
            actual = m.christoffel_tensor(directions)
            self.assertEqual(actual.shape, expected.shape)
            # Contractions can nearly cancel, so use the tensor's stiffness scale.
            scale = np.max(np.abs(expected))
            self.assertArraysAlmostEqual(
                (actual / scale).flatten(),
                (expected / scale).flatten(),
                tol=1.0e-6,
                tol_zero=1.0e-7,
            )

    def test_prescribed_isochoric_strain(self):
        """
        Recover the imposed Hencky strain without changing the volume.
        """
        m = self._get_hyperelastic_mineral()
        prescribed_strain = self._set_finite_strain_state(m)
        self.assertArraysAlmostEqual(
            m.hencky_strain.flatten(),
            prescribed_strain.flatten(),
            tol=1.0e-12,
            tol_zero=1.0e-14,
        )
        self.assertFloatEqual(np.trace(m.hencky_strain), 0.0, tol_zero=1.0e-12)
        self.assertFloatEqual(m.V, m.hydrostatic.V, tol=1.0e-12, tol_zero=0.0)

    def test_rigid_rotation(self):
        """
        A rigid rotation preserves energy and rotates stress and stiffness.
        """
        m1 = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m1)
        T = m1.T
        cell_vectors_m1 = m1.cell_vectors.copy()

        # Create another mineral with the same finite-strain state
        # but rotated cell vectors.
        angle = 0.37
        c = np.cos(angle)
        s = np.sin(angle)

        Q = np.array(
            [
                [c, s, 0.0],
                [-s, c, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )

        m2 = self._get_hyperelastic_mineral()
        m2.set_state_with_molar_cell_vectors((Q @ cell_vectors_m1.T).T, T)

        self.assertFloatEqual(m1.helmholtz, m2.helmholtz)

        m1_stress_rotated = Q @ m1.cauchy_stress @ Q.T
        self.assertArraysAlmostEqual(
            m1_stress_rotated.flatten(), m2.cauchy_stress.flatten(), tol_zero=1e-12
        )

        m1_stiffness_rotated = np.einsum(
            "im,jn,ko,lp,mnop->ijkl",
            Q,
            Q,
            Q,
            Q,
            m1.full_isothermal_stiffness_tensor,
        )
        m1_stiffness_rotated_voigt = contract_stiffnesses(m1_stiffness_rotated)

        self.assertArraysAlmostEqual(
            m1_stiffness_rotated_voigt.flatten(),
            m2.isothermal_stiffness_tensor.flatten(),
            tol_zero=1e-12,
            tol=1.0e-4,
        )

    def test_stiffness_from_total_energy(self):
        """
        Check normal and shear stiffness columns using total-energy differences.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)
        T = m.T

        cell_vectors_0 = m.cell_vectors.copy()
        M_0 = m._M_unrotated.copy()

        de_stress = 1.0e-4
        de_tangent = 1.0e-4

        def stress_voigt_from_M(M_base):
            """
            Differentiate total Helmholtz energy in the fixed axes of M_base.
            """
            V_base = np.linalg.det(M_base)
            stress = np.zeros(6)

            for i in range(6):
                deps = np.zeros(6)
                deps[i] = de_stress
                deps_full = expand_strains(deps)

                M_plus = expm(deps_full) @ M_base
                M_minus = expm(-deps_full) @ M_base

                m.set_state_with_molar_cell_vectors(M_plus.T, T)
                helmholtz_plus = m.helmholtz

                m.set_state_with_molar_cell_vectors(M_minus.T, T)
                helmholtz_minus = m.helmholtz

                stress[i] = (helmholtz_plus - helmholtz_minus) / (
                    2.0 * de_stress * V_base
                )

            return stress

        # A normal and a shear column are sufficient to expose most
        # index/order mistakes without making the test unnecessarily slow.
        test_columns = [0, 3]
        stiffness_direct = np.zeros((6, len(test_columns)))

        for k, j in enumerate(test_columns):
            deps = np.zeros(6)
            deps[j] = de_tangent
            deps_full = expand_strains(deps)

            M_plus = expm(deps_full) @ M_0
            M_minus = expm(-deps_full) @ M_0

            stress_plus = stress_voigt_from_M(M_plus)
            stress_minus = stress_voigt_from_M(M_minus)

            stiffness_direct[:, k] = (stress_plus - stress_minus) / (2.0 * de_tangent)

        # Restore the original finite-strain state before obtaining the
        # stiffness from HyperelasticMineral.
        m.set_state_with_molar_cell_vectors(cell_vectors_0, T)

        stiffness = m._isothermal_stiffness_tensor_unrotated

        self.assertArraysAlmostEqual(
            (stiffness[:, test_columns] / 1.0e9).flatten(),
            (stiffness_direct / 1.0e9).flatten(),
            tol_zero=1.0e-3,
        )

    def test_christoffel_from_adiabatic_energy(self):
        """
        Check the seismic tensor against wave energy at constant entropy.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)
        cell_vectors, temperature = m.cell_vectors.copy(), 1000.0
        m.set_state_with_molar_cell_vectors(cell_vectors, temperature)
        entropy, volume = m.entropy, m.V
        direction = np.array([1.0, 2.0, 3.0]) / np.sqrt(14.0)
        christoffel = m.christoffel_tensor(direction)

        def energy(amplitude):
            # A plane wave has displacement gradient amplitude (x) direction.
            deformed = ((np.eye(3) + np.outer(amplitude, direction)) @ cell_vectors.T).T

            def entropy_residual(T):
                m.set_state_with_molar_cell_vectors(deformed, T)
                return m.entropy - entropy

            T = brentq(
                entropy_residual,
                temperature - 10.0,
                temperature + 10.0,
                xtol=1.0e-9,
            )
            m.set_state_with_molar_cell_vectors(deformed, T)
            return m.internal_energy

        step = 2.0e-4
        expected = np.empty((3, 3))
        try:
            for i in range(3):
                for k in range(i + 1):
                    a, b = step * np.eye(3)[i], step * np.eye(3)[k]
                    expected[i, k] = expected[k, i] = (
                        energy(a + b) - energy(a - b) - energy(-a + b) + energy(-a - b)
                    ) / (4.0 * step**2 * volume)
        finally:
            m.set_state_with_molar_cell_vectors(cell_vectors, temperature)

        self.assertArraysAlmostEqual(
            (christoffel / 1.0e9).flatten(),
            (expected / 1.0e9).flatten(),
            tol=1.0e-5,
            tol_zero=1.0e-4,
        )
        self.assertArraysAlmostEqual(
            christoffel.flatten(),
            christoffel.T.flatten(),
            tol=1.0e-12,
        )
        velocities, _ = m.wave_velocities(direction)
        self.assertArraysAlmostEqual(
            velocities,
            np.sqrt(np.linalg.eigvalsh(expected)[::-1] / m.density),
            tol=1.0e-5,
        )

    def test_thermal_stress_at_fixed_cell(self):
        """
        Check thermal stress against the stress derivative at fixed cell vectors.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)

        cell_vectors_0 = m.cell_vectors.copy()
        T_0 = m.T

        dT = 0.1

        m.set_state_with_molar_cell_vectors(
            cell_vectors_0,
            T_0 + dT,
        )
        stress_plus = m.cauchy_stress.copy()

        m.set_state_with_molar_cell_vectors(
            cell_vectors_0,
            T_0 - dT,
        )
        stress_minus = m.cauchy_stress.copy()

        thermal_stress_direct = (stress_plus - stress_minus) / (2.0 * dT)

        m.set_state_with_molar_cell_vectors(
            cell_vectors_0,
            T_0,
        )

        thermal_stress = m.thermal_stress

        # Scale to MPa/K to keep the numerical tolerance meaningful.
        self.assertArraysAlmostEqual(
            (thermal_stress / 1.0e6).flatten(),
            (thermal_stress_direct / 1.0e6).flatten(),
            tol_zero=1.0e-3,
        )

    def test_pressure(self):
        """
        Pressure is finite and equals minus one-third of the stress trace.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)

        self.assertTrue(np.isfinite(m.pressure))

        self.assertAlmostEqual(
            m.pressure / 1.0e9,
            -np.trace(m._cauchy_stress_unrotated) / 3.0e9,
            places=4,
        )

        self.assertAlmostEqual(
            m.pressure / 1.0e9,
            -np.trace(m.cauchy_stress) / 3.0e9,
            places=10,
        )

    def test_finite_strain_heat_capacity_step_convergence(self):
        """
        Check heat capacity against energy curvature at several temperature steps.
        """
        m = self._get_hyperelastic_mineral()
        F = np.array([[0.9, 0.1, 0.05], [0.1, 0.91, 0.04], [0.02, 0.03, 0.89]])
        cell_vectors = (F @ m.hydrostatic.cell_vectors_0.T).T
        T = 300.0
        m.set_state_with_molar_cell_vectors(cell_vectors, T)
        energy = m.helmholtz
        heat_capacity = m.molar_isometric_heat_capacity

        try:
            for dT in [0.5, 0.2, 0.1]:
                with self.subTest(dT=dT):
                    m.set_state_with_molar_cell_vectors(cell_vectors, T + dT)
                    plus = m.helmholtz
                    m.set_state_with_molar_cell_vectors(cell_vectors, T - dT)
                    minus = m.helmholtz
                    curvature = -T * (plus - 2.0 * energy + minus) / dT**2
                    self.assertFloatEqual(curvature, heat_capacity, tol=2.0e-5)
        finally:
            m.set_state_with_molar_cell_vectors(cell_vectors, T)

    def test_derivatives_independent_of_reference_energy(self):
        """
        Adding a constant to the energy must leave all derivatives unchanged.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)
        cells, T = m.cell_vectors.copy(), m.T
        properties = [
            "cauchy_stress",
            "isothermal_stiffness_tensor",
            "entropy",
            "molar_isometric_heat_capacity",
            "thermal_stress",
        ]
        expected = {p: np.array(getattr(m, p), copy=True) for p in properties}
        energy = m.helmholtz
        reference_energy = m.params["F_0"]
        try:
            m.params["F_0"] += 1.0e12
            m.set_state_with_molar_cell_vectors(cells, T)
            self.assertAlmostEqual((m.helmholtz - energy) / 1.0e12, 1.0)
            for prop in properties:
                with self.subTest(property=prop):
                    self.assertArraysAlmostEqual(
                        np.asarray(getattr(m, prop)).ravel(),
                        expected[prop].ravel(),
                        tol=1.0e-12,
                        tol_zero=1.0e-12,
                    )
        finally:
            m.params["F_0"] = reference_energy
            m.set_state_with_molar_cell_vectors(cells, T)

    def test_eos_consistency(self):
        """
        Check thermodynamic consistency across thermal and output-frame cases.
        """
        for thermal_shear, output_frame, temperature in [
            (False, "current", 300.0),
            (True, "current", 1000.0),
            (False, "crystal", 1000.0),
        ]:
            with self.subTest(
                thermal_shear=thermal_shear, frame=output_frame, T=temperature
            ):
                m = self._get_hyperelastic_mineral(thermal_shear=thermal_shear)
                m.set_output_frame(output_frame)
                self.assertTrue(check_hyperelastic_eos_consistency(m, T=temperature))

    def test_output_frame_switch(self):
        """
        Switching output frames rotates cached tensors and preserves the state.
        """
        m = self._get_hyperelastic_mineral()
        self._set_finite_strain_state(m)
        cells = m.cell_vectors.copy()
        Q = m.rotation_to_crystallographic_frame.copy()
        stress = m.cauchy_stress.copy()
        stiffness = m.full_isothermal_stiffness_tensor.copy()
        directions = np.array([[1.0, 2.0, 3.0], [-2.0, 1.0, 4.0]])
        christoffel = m.christoffel_tensor(directions)
        # Populate these caches immediately before changing the frame.
        rotation = m.rotation_matrix.copy()
        m.cell_vectors
        self.assertEqual(m.output_frame, "current")

        m.set_output_frame("crystal")
        self.assertEqual(m.output_frame, "crystal")
        self.assertArraysAlmostEqual(
            m.cell_vectors.flatten(), (Q @ cells.T).T.flatten(), tol=1.0e-12
        )
        self.assertArraysAlmostEqual(
            m.rotation_matrix.flatten(), (Q @ rotation).flatten(), tol=1.0e-12
        )
        self.assertArraysAlmostEqual(
            m.cauchy_stress.flatten(), (Q @ stress @ Q.T).flatten(), tol=1.0e-10
        )
        expected_stiffness = np.einsum("ia,jb,kc,ld,abcd->ijkl", Q, Q, Q, Q, stiffness)
        self.assertArraysAlmostEqual(
            m.full_isothermal_stiffness_tensor.flatten(),
            expected_stiffness.flatten(),
            tol=1.0e-10,
        )
        self.assertArraysAlmostEqual(
            m.christoffel_tensor(directions @ Q.T).flatten(),
            (Q @ christoffel @ Q.T).flatten(),
            tol=1.0e-10,
        )

        m.set_output_frame("current")
        self.assertEqual(m.output_frame, "current")
        self.assertArraysAlmostEqual(
            m.cell_vectors.flatten(), cells.flatten(), tol=0.0, tol_zero=0.0
        )
        self.assertArraysAlmostEqual(
            m.cauchy_stress.flatten(), stress.flatten(), tol=1.0e-10
        )

    def test_invalid_output_frame(self):
        """
        Reject an unknown output frame without changing the selected frame.
        """
        m = self._get_hyperelastic_mineral()
        m.set_output_frame("crystal")
        with self.assertRaises(ValueError):
            m.set_output_frame("invalid")
        self.assertEqual(m.output_frame, "crystal")


if __name__ == "__main__":
    unittest.main()
