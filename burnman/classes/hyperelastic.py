# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.linalg import expm
from burnman.utils.unitcell import cell_vectors_to_parameters
from burnman.utils.unitcell import rotation_to_crystallographic_frame
from burnman.utils.anisotropy import expand_strains, expand_stresses
from burnman.utils.anisotropy import contract_stiffnesses, contract_compliances
from burnman.utils.anisotropy import (
    voigt_notation_to_compliance_tensor,
    voigt_notation_to_stiffness_tensor,
)
from burnman.utils.anisotropy import polar_decomposition
from .mineral import Mineral
from .material import material_property
from .anisotropicmineral import AnisotropicMineral
from enum import Enum


class OutputFrame(Enum):
    crystal = 1
    current = 2


class HyperelasticMineral(AnisotropicMineral):
    """
    This class represents a hyperelastic mineral,
    where the word "hyperelastic" indicates that the mineral
    can undergo large elastic deformations to non-hydrostatic states.

    The hydrostatic state of the mineral is represented by an
    instance of the AnisotropicMineral class. The non-hydrostatic
    elastic energy at fixed volume and temperature
    is then calculated as 1/2 * V * H : C^hyd : H,
    where V is the volume of the mineral,
    H is the right Hencky strain tensor describing the isochoric
    deformation of the mineral away from the hydrostatic state,
    and C^hyd is the *spatial* isothermal stiffness tensor of the mineral
    at the hydrostatic state.
    All other thermodynamic properties are calculated as thermodynamic
    derivatives of the total Helmholtz free energy.

    Despite using the *spatial* stiffness tensor, it can be shown that
    the elastic energy term corresponds to treating the
    isochoric components of the *material* isothermal stiffness tensor
    as constant. The non-isochoric components of the
    material stiffness tensor are not constant,
    and neither is the spatial stiffness tensor.
    Indeed, the spatial stiffness tensor does not even retain the
    major symmetry ijkl = klij under non-hydrostatic conditions.

    Setting the state of a HyperelasticMineral instance requires
    specifying the cell vectors and temperature via
    the set_state_with_molar_cell_vectors method.

    Properties of the mineral can be reported in
    either the current configuration defined by the cell vectors, or in the
    crystallographic frame defined by the cell vectors and frame convention
    of the underlying AnisotropicMineral. By default, properties are reported
    in the current configuration, but the set_output_frame method can
    be used to toggle between "current" and "crystal" output frames.

    The stiffness tensors output by this class are
    derivatives of the Cauchy stress with respect to the strain.
    Under non-hydrostatic conditions, these are not the tensors
    directly relevant for seismic wave propagation.
    The christoffel_tensor method of the AnisotropicMineral class
    applies the required correction, which includes terms
    involving the Cauchy stress
    (Wallace, 1967; doi:10.1103/PhysRev.162.776).
    """

    def __init__(self, hydrostatic_mineral):
        """
        Initialize the hyperelastic mineral with a hydrostatic mineral.
        """
        assert isinstance(
            hydrostatic_mineral, AnisotropicMineral
        ), "The hydrostatic_mineral must be an instance of AnisotropicMineral."
        self.hydrostatic = hydrostatic_mineral

        # minerals in non-hydrostatic states are generally non-orthotropic,
        # even if the hydrostatic mineral is orthotropic.
        self.orthotropic = False
        self.frame_convention = self.hydrostatic.frame_convention
        self._output_frame = OutputFrame.current
        Mineral.__init__(
            self, self.hydrostatic.params, self.hydrostatic.property_modifiers
        )

    def set_output_frame(self, output_frame):
        """
        Set the output frame for the hyperelastic mineral.
        This changes the reported cell vectors and spatial tensor components,
        but does not change the underlying state of the mineral.

        :param output_frame: str
            The output frame to use. Must be one of
            'crystal' or 'current'.
        """
        if output_frame == "crystal":
            self._output_frame = OutputFrame.crystal
        elif output_frame == "current":
            self._output_frame = OutputFrame.current
        else:
            raise ValueError(
                "Invalid output frame. Must be one of 'crystal' or 'current'."
            )
        # Cell vectors, rotations and constitutive tensors are cached material
        # properties and must be recomputed in the newly selected frame.
        self.reset()

    @material_property
    def output_frame(self):
        """
        Get the output frame for the hyperelastic mineral.
        The output frame is used to define the orientation of the
        cell vectors and the other tensor properties.

        :returns: The output frame for the hyperelastic mineral.
        :rtype: str
        """
        return self._output_frame.name

    def set_state_with_molar_cell_vectors(self, cell_vectors, T):
        """
        Set the state of the mineral given cell vectors and temperature.

        :param cell_vectors : array-like (3x3)
            The cell vectors of the mineral in the current spatial frame,
            independent of the selected output_frame.
        :param T : temperature (K)
            The temperature of the mineral.
        """
        self._cell_vectors_current = cell_vectors
        M = cell_vectors.T

        # set the state of the hydrostatic mineral
        self._V = np.linalg.det(M)
        self.hydrostatic.set_state_with_volume(self._V, T)

        # Collect the "unrotated" cell vectors of the hydrostatic mineral
        # Note that these vectors may be rotated relative to the
        # cell vectors at the reference pressure and temperature,
        # with that rotation due to the hydrostatic deformation
        # of the mineral. However, they are neither rotated into the
        # current configuration nor into the crystallographic frame.
        M_hyd = self.hydrostatic.unrotated_cell_vectors.T

        # Compute the deformation gradient F between the
        # hydrostatic reference state and the current state
        F = M @ np.linalg.inv(M_hyd)

        if not np.isclose(np.linalg.det(F), 1.0, atol=1e-6):
            print(
                "Warning: Deformation gradient is not volume-preserving "
                f"(det(F)={np.linalg.det(F):.6f})"
            )

        # Perform polar decomposition of F to obtain the rotation matrix
        # and Hencky strain
        self._isochoric_rotation_matrix, self._isochoric_hencky_strain = (
            polar_decomposition(F, "H")
        )

        # Store the current M with the isochoric rotation matrix removed
        self._M_unrotated = self._isochoric_rotation_matrix.T @ M

        # Check that the Hencky strain is traceless,
        # as it should be for an isochoric deformation
        if not np.isclose(np.trace(self._isochoric_hencky_strain), 0.0, atol=1e-8):
            print(
                "Warning: Hencky strain is not traceless "
                f"(trace(H)={np.trace(self._isochoric_hencky_strain):.6e})"
            )

        # Set the state of the hyperelastic mineral with the temperature.
        # We set the pressure to np.nan; there is a separate method
        # to compute it from the Cauchy stress tensor.
        Mineral.set_state(self, np.nan, T)

    @material_property
    def rotation_to_crystallographic_frame(self):
        """
        :returns: The rotation matrix that rotates the mineral from its
            current configuration into the crystallographic frame.
        :rtype: numpy.array (2D)
        """
        Q = rotation_to_crystallographic_frame(
            self._cell_vectors_current, self.frame_convention
        )
        return Q

    @material_property
    def rotation_matrix(self):
        """
        :returns: The rotation matrix that rotates the mineral from
            its unrotated configuration to the configuration
            indicated by output_frame. If output_frame is 'crystal',
            this is the rotation matrix that rotates the mineral
            into the crystallographic frame. If output_frame is 'current',
            this is the rotation matrix that rotates the mineral
            into the current configuration defined by the cell vectors
            passed as argument to set_state_with_molar_cell_vectors.
        :rtype: numpy.array (2D)
        """
        if self._output_frame == OutputFrame.crystal:
            Q = self.rotation_to_crystallographic_frame
            return Q @ self._isochoric_rotation_matrix
        elif self._output_frame == OutputFrame.current:
            return self._isochoric_rotation_matrix
        else:
            raise ValueError(
                "Invalid output frame. Must be one of 'crystal' or 'current'."
            )

    @material_property
    def unrotated_cell_vectors(self):
        """
        :returns: The cell vectors of the hyperelastic mineral in the
            deformed state, rotated out of the current frame by
            removing the isochoric rotation.

            These are equal to the cell vectors of the hydrostatic mineral
            at the same volume and temperature, then deformed
            by the isochoric Hencky strain without any additional rotation.
        :rtype: numpy.array (2D)
        """
        return self._M_unrotated.T

    @material_property
    def cell_vectors(self):
        """
        :returns: The cell vectors of the deformed mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (3D)
        """
        if self._output_frame == OutputFrame.crystal:
            Q = self.rotation_to_crystallographic_frame
            return (Q @ self._cell_vectors_current.T).T
        elif self._output_frame == OutputFrame.current:
            return self._cell_vectors_current
        else:
            raise ValueError(
                "Invalid output frame. Must be one of 'crystal' or 'current'."
            )

    @material_property
    def deformation_gradient_tensor(self):
        """
        :returns: The deformation gradient tensor describing the deformation
            of the mineral from its undeformed state
            (i.e. the state at the reference pressure and temperature)
            to its current state, in the configuration indicated by
            output_frame (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        M_0 = self.hydrostatic.cell_vectors_0.T
        M = self.cell_vectors.T
        F = M @ np.linalg.inv(M_0)
        return F

    @material_property
    def cell_parameters(self):
        """
        :returns: The molar cell parameters of the mineral, given in standard form:
            [:math:`a`, :math:`b`, :math:`c`,
            :math:`\\alpha`, :math:`\\beta`, :math:`\\gamma`],
            where the first three floats are the lengths of the vectors in [m]
            defining the cell constructed from one mole of formula units.
            The last three floats are angles between vectors
            (given in radians). See the documentation for the function
            :func:`~burnman.utils.unitcell.cell_parameters_to_vectors`
            for the assumed relationships between the cell vectors and
            spatial coordinate axes.
        :rtype: numpy.array (1D)
        """
        Q = self.rotation_to_crystallographic_frame
        vecs_rotated = (Q @ self._cell_vectors_current.T).T
        return cell_vectors_to_parameters(vecs_rotated, self.frame_convention)

    @material_property
    def hencky_strain(self):
        """
        :returns: The right isochoric Hencky strain tensor of the mineral,
            which is a measure of the constant volume, constant temperature
            deformation of the mineral away from a hydrostatic state.
        :rtype: numpy.array (2D)
        """
        return self._isochoric_hencky_strain

    @material_property
    def _helmholtz_elastic(self):
        """
        :returns: The elastic correction to the Helmholtz free energy.
        :rtype: float
        """
        H = self.hencky_strain
        CT = voigt_notation_to_stiffness_tensor(
            self.hydrostatic._unrotated_isothermal_stiffness_tensor
        )

        helmholtz_elastic = (
            0.5 * self.hydrostatic.V * np.einsum("ij, ijkl, kl", H, CT, H)
        )
        return helmholtz_elastic

    @material_property
    def helmholtz(self):
        """
        :returns: The total Helmholtz free energy of the deformed mineral,
        including the nonhydrostatic isochoric elastic component.
        :rtype: float
        """
        return self.hydrostatic.helmholtz + self._helmholtz_elastic

    @material_property
    def _cauchy_stress_unrotated(self):
        """
        :returns: The Cauchy stress tensor in the unrotated configuration;
            i.e., the configuration defined by the unrotated hydrostatic
            configuration modified by the isochoric Hencky strain
            but not rotated into the current configuration or the
            crystallographic frame.
        :rtype: numpy.array (2D)
        """
        # First, store the current cell vectors and temperature
        cell_vectors_current = self._cell_vectors_current.copy()
        M_0 = self._M_unrotated.copy()
        T_0 = self.T
        V = self.hydrostatic.V
        pressure = self.hydrostatic.pressure

        # Apply a small perturbation to the cell vectors
        # to obtain the stress tensor
        de = 1.0e-4  # small strain perturbation
        stress_voigt = np.zeros(6)
        for i in range(6):
            deps = np.zeros(6)
            deps[i] = de  # small strain perturbation
            deps_full = expand_strains(deps)
            M_plus = expm(deps_full) @ M_0
            M_minus = expm(-deps_full) @ M_0

            self.set_state_with_molar_cell_vectors(M_plus.T, T_0)
            helmholtz_plus = self._helmholtz_elastic
            self.set_state_with_molar_cell_vectors(M_minus.T, T_0)
            helmholtz_minus = self._helmholtz_elastic
            stress_component = (helmholtz_plus - helmholtz_minus) / (2 * de * V)
            stress_voigt[i] = stress_component

        # Restore the original state
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0)
        # dF_hydrostatic / dstrain = -P V I is known analytically.
        return expand_stresses(stress_voigt) - pressure * np.eye(3)

    @material_property
    def cauchy_stress(self):
        """
        :returns: The Cauchy stress tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        unrotated_stress = self._cauchy_stress_unrotated
        rotated_stress = (
            self.rotation_matrix @ unrotated_stress @ self.rotation_matrix.T
        )
        return rotated_stress

    @material_property
    def pressure(self):
        """
        :returns: The pressure, computed as the negative one-third
            of the trace of the Cauchy stress.
        :rtype: float
        """
        return -np.trace(self.cauchy_stress) / 3.0

    @material_property
    def _isothermal_stiffness_tensor_unrotated(self):
        """
        :returns: The isothermal stiffness tensor in Voigt notation (6x6),
            computed in the unrotated configuration,
            i.e., the configuration defined by the unrotated hydrostatic
            configuration modified by the isochoric Hencky strain
            but not rotated into the current configuration or the
            crystallographic frame.
        :rtype: numpy.array (2D)
        """

        # First, store the current cell vectors and temperature
        cell_vectors_current = self._cell_vectors_current.copy()
        M_0 = self._M_unrotated.copy()
        T_0 = self.T
        V = self.hydrostatic.V
        pressure = self.hydrostatic.pressure
        bulk_modulus = self.hydrostatic.isothermal_bulk_modulus_reuss

        de = 1.0e-4  # small strain perturbation
        stiffness = np.zeros((6, 6))
        for i in range(6):
            for j in range(6):
                deps_i = np.zeros(6)
                deps_j = np.zeros(6)
                deps_i[i] = de
                deps_j[j] = de

                deps_i_full = expand_strains(deps_i)
                deps_j_full = expand_strains(deps_j)

                # Strains are not commutative
                M_plus_plus = expm(deps_i_full) @ expm(deps_j_full) @ M_0
                M_plus_minus = expm(deps_i_full) @ expm(-deps_j_full) @ M_0
                M_minus_plus = expm(-deps_i_full) @ expm(deps_j_full) @ M_0
                M_minus_minus = expm(-deps_i_full) @ expm(-deps_j_full) @ M_0

                self.set_state_with_molar_cell_vectors(M_plus_plus.T, T_0)
                helmholtz_plus_plus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_plus_minus.T, T_0)
                helmholtz_plus_minus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_minus_plus.T, T_0)
                helmholtz_minus_plus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_minus_minus.T, T_0)
                helmholtz_minus_minus = self._helmholtz_elastic

                stiffness[i, j] = (
                    helmholtz_plus_plus
                    - helmholtz_plus_minus
                    - helmholtz_minus_plus
                    + helmholtz_minus_minus
                ) / (
                    4.0 * de * de * V
                )  # Convert from energy per mole to energy per unit volume

        # Now restore the original state
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0)

        # The hydrostatic tangent is K I (x) I. Only the elastic stress
        # contributes to the remaining derivative of the 1/V prefactor.
        stiffness[:3, :3] += bulk_modulus
        stress = self._cauchy_stress_unrotated + pressure * np.eye(3)
        stiffness2 = -np.einsum("kl, ij -> ijkl", np.eye(3), stress)
        stiffness += contract_stiffnesses(stiffness2)

        return stiffness

    @material_property
    def _full_isothermal_stiffness_tensor_unrotated(self):
        """
        :returns: The full isothermal stiffness tensor of the mineral,
            computed in the unrotated configuration,
            i.e., the configuration defined by the unrotated hydrostatic
            configuration modified by the isochoric Hencky strain
            but not rotated into the current configuration or the
            crystallographic frame.
        :rtype: numpy.array (4D)
        """

        stiffness_voigt = self._isothermal_stiffness_tensor_unrotated
        full_stiffness = voigt_notation_to_stiffness_tensor(stiffness_voigt)
        return full_stiffness

    @material_property
    def full_isothermal_stiffness_tensor(self):
        """
        :returns: The full isothermal stiffness tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (4D)
        """

        unrotated_stiffness = self._full_isothermal_stiffness_tensor_unrotated

        # We must apply the rotation to all four indices of the stiffness tensor, not just the first two.
        # This is done by converting the 3x3x3x3 tensor to a 6x6 Voigt notation and then rotating it.
        rotated_stiffness = np.einsum(
            "im, jn, ko, lp, mnop -> ijkl",
            self.rotation_matrix,
            self.rotation_matrix,
            self.rotation_matrix,
            self.rotation_matrix,
            unrotated_stiffness,
        )
        return rotated_stiffness

    @material_property
    def isothermal_stiffness_tensor(self):
        """
        :returns: The isothermal stiffness tensor of the mineral
            in Voigt notation, reported in the configuration
            indicated by output_frame (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """

        full_stiffness = self.full_isothermal_stiffness_tensor
        stiffness_voigt = contract_stiffnesses(full_stiffness)
        return stiffness_voigt

    @material_property
    def entropy(self):
        """
        :returns: The molar entropy of the mineral.
        :rtype: float
        """

        # First, store the current cell vectors and temperature by explicit copy
        cell_vectors_current = self._cell_vectors_current.copy()
        T_0 = self.T
        hydrostatic_entropy = self.hydrostatic.entropy

        dT = 0.1  # K
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0 + dT)
        helmholtz_plus = self._helmholtz_elastic
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0 - dT)
        helmholtz_minus = self._helmholtz_elastic
        entropy = hydrostatic_entropy - (helmholtz_plus - helmholtz_minus) / (2 * dT)

        # Restore the original state
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0)
        return entropy

    @material_property
    def molar_isometric_heat_capacity(self):
        """
        :returns: The molar isometric heat capacity of the mineral.
        :rtype: float
        """

        # First, store the current cell vectors and temperature by explicit copy
        cell_vectors_current = self._cell_vectors_current.copy()
        T_0 = self.T
        # The scalar contribution is at fixed volume. The elastic correction
        # accounts for constraining cell shape as well as volume.
        hydrostatic_Cv = self.hydrostatic.molar_heat_capacity_v

        dT = 0.1  # K
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0 + dT)
        helmholtz_plus = self._helmholtz_elastic
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0 - dT)
        helmholtz_minus = self._helmholtz_elastic
        self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0)
        helmholtz = self._helmholtz_elastic

        d2F_dT2 = (helmholtz_plus - 2.0 * helmholtz + helmholtz_minus) / (dT * dT)
        return hydrostatic_Cv - T_0 * d2F_dT2

    @material_property
    def _thermal_stress_unrotated(self):
        """
        :return: The thermal stress of the mineral in the
            unrotated configuration,
            i.e., the configuration defined by the unrotated hydrostatic
            configuration modified by the isochoric Hencky strain
            but not rotated into the current configuration or the
            crystallographic frame.

        :rtype: numpy.array (2D)
        """
        cell_vectors_current = self._cell_vectors_current.copy()
        M_0 = self._M_unrotated.copy()
        T_0 = self.T
        V = self.hydrostatic.V
        aK = (
            self.hydrostatic.thermal_expansivity
            * self.hydrostatic.isothermal_bulk_modulus_reuss
        )

        # The mixed Helmholtz derivative is taken at fixed cell vectors in
        # this frame.
        de = 1.0e-4
        dT = 0.1  # K
        thermal_stress_voigt = np.zeros(6)
        try:
            for i in range(6):
                deps = np.zeros(6)
                deps[i] = de
                deps_full = expand_strains(deps)
                M_plus = expm(deps_full) @ M_0
                M_minus = expm(-deps_full) @ M_0

                # Use the same strain directions and cell matrices at both
                # temperatures: pi = (1/V) d^2F/(dstrain dT) = -(1/V) dS/dstrain.
                self.set_state_with_molar_cell_vectors(M_plus.T, T_0 + dT)
                plus_plus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_plus.T, T_0 - dT)
                plus_minus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_minus.T, T_0 + dT)
                minus_plus = self._helmholtz_elastic
                self.set_state_with_molar_cell_vectors(M_minus.T, T_0 - dT)
                minus_minus = self._helmholtz_elastic
                thermal_stress_voigt[i] = (
                    (plus_plus - plus_minus) - (minus_plus - minus_minus)
                ) / (4.0 * de * dT * V)

            # The scalar contribution is -dP/dT at fixed volume.
            return expand_stresses(thermal_stress_voigt) - aK * np.eye(3)
        finally:
            self.set_state_with_molar_cell_vectors(cell_vectors_current, T_0)

    @material_property
    def thermal_stress(self):
        """
        :returns: The thermal stress tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        unrotated_stress = self._thermal_stress_unrotated
        rotated_stress = (
            self.rotation_matrix @ unrotated_stress @ self.rotation_matrix.T
        )
        return rotated_stress

    @material_property
    def isothermal_compliance_tensor(self):
        """
        :returns: The isothermal compliance tensor of the mineral
            in Voigt notation, reported in the configuration indicated
            by output_frame (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        stiffness = self.isothermal_stiffness_tensor
        compliance = np.linalg.inv(stiffness)
        return compliance

    @material_property
    def thermal_expansivity_tensor(self):
        """
        :returns: The thermal expansivity tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        thermal_stress = self.thermal_stress
        full_compliance = self.full_isothermal_compliance_tensor

        thermal_expansivity_tensor = -np.einsum(
            "ijkl, kl -> ij", full_compliance, thermal_stress
        )
        return thermal_expansivity_tensor

    @material_property
    def molar_isostress_heat_capacity(self):
        """
        :returns: The molar isostress heat capacity of the mineral.
        :rtype: float
        """
        alpha = self.thermal_expansivity_tensor
        pi = self.thermal_stress
        Ceps = self.molar_isometric_heat_capacity
        Csigma = Ceps - self.V * self.T * np.einsum("ij, ij", alpha, pi)
        return Csigma

    @material_property
    def full_isothermal_compliance_tensor(self):
        """
        :returns: The full isothermal compliance tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (4D)
        """
        stiffness = self.isothermal_stiffness_tensor
        compliance = np.linalg.inv(stiffness)
        return voigt_notation_to_compliance_tensor(compliance)

    @material_property
    def full_isentropic_compliance_tensor(self):
        """
        :returns: The full isentropic compliance tensor of the mineral,
            reported in the configuration indicated by output_frame
            (either 'crystal' or 'current').
        :rtype: numpy.array (4D)
        """
        full_isothermal_compliance = self.full_isothermal_compliance_tensor
        alpha = self.thermal_expansivity_tensor
        # C_T (and its inverse) need not have major symmetry at finite
        # nonhydrostatic strain. The Sherman-Morrison formula therefore
        # requires distinct left and right contractions for alpha.
        alpha_right = -np.einsum(
            "ij,ijkl->kl", self.thermal_stress, full_isothermal_compliance
        )
        V = self.V
        C_sigma = self.molar_isostress_heat_capacity
        T = self.T

        isentropic_compliance = (
            full_isothermal_compliance
            - np.einsum("ij, kl -> ijkl", alpha, alpha_right) * V * T / C_sigma
        )
        return isentropic_compliance

    @material_property
    def isentropic_compliance_tensor(self):
        """
        :returns: The isentropic compliance tensor of the mineral
            in Voigt notation, reported in the configuration
            indicated by output_frame (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        S_full = self.full_isentropic_compliance_tensor
        return contract_compliances(S_full)

    @material_property
    def isentropic_stiffness_tensor(self):
        """
        :returns: The isentropic stiffness tensor of the mineral
            in Voigt notation, reported in the configuration
            indicated by output_frame (either 'crystal' or 'current').
        :rtype: numpy.array (2D)
        """
        S_S = self.isentropic_compliance_tensor
        stiffness = np.linalg.inv(S_S)
        return stiffness

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        :returns: The Reuss isothermal compressibility of the mineral.
        :rtype: float
        """
        beta = np.sum(
            [
                [self.isothermal_compliance_tensor[i][k] for k in range(3)]
                for i in range(3)
            ]
        )
        return beta

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        :returns: The Reuss isothermal bulk modulus of the mineral.
        :rtype: float
        """
        return 1.0 / self.isothermal_compressibility_reuss

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        :returns: The Reuss isentropic compressibility of the mineral.
        :rtype: float
        """
        beta = np.sum(
            [
                [self.isentropic_compliance_tensor[i][k] for k in range(3)]
                for i in range(3)
            ]
        )
        return beta

    @material_property
    def _molar_volume_unmodified(self):
        """
        :returns: The molar volume of the mineral,
            unmodified by any property modifiers.
        :rtype: float
        """
        return self._V

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        :returns: The Reuss isentropic bulk modulus of the mineral.
        :rtype: float
        """
        return 1.0 / self.isentropic_compressibility_reuss

    @material_property
    def internal_energy(self):
        """
        :returns: The molar internal energy of the mineral.
        :rtype: float
        """
        helmholtz = self.helmholtz
        entropy = self.entropy
        internal_energy = helmholtz + self.T * entropy
        return internal_energy
