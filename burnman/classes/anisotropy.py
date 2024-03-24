# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np

from ..utils.math import unit_normalize
from .material import Material, material_property
from ..utils.anisotropy import (
    voigt_array_from_cijs,
    voigt_notation_to_compliance_tensor,
    voigt_notation_to_stiffness_tensor,
)


class AnisotropicMaterial(Material):
    """
    A base class for anisotropic elastic materials. The base class
    is initialised with a density and a full isentropic stiffness tensor
    in Voigt notation. It can then be interrogated to find the values of
    different properties, such as bounds on seismic velocities.
    There are also several functions which can be called to calculate
    properties along directions oriented with respect to the isentropic
    elastic tensor.

    See :cite:`Mainprice2011`
    and https://docs.materialsproject.org/methodology/elasticity/
    for mathematical descriptions of each function.
    """

    def __init__(self, rho, cijs):
        self._isentropic_stiffness_tensor = cijs
        self._rho = rho

        assert cijs.shape == (6, 6), "cijs must be in Voigt notation (6x6)"
        assert np.allclose(cijs.T, cijs), "stiffness_tensor must be symmetric"

        Material.__init__(self)

    @material_property
    def isentropic_stiffness_tensor(self):
        return self._isentropic_stiffness_tensor

    @material_property
    def full_isentropic_stiffness_tensor(self):
        return voigt_notation_to_stiffness_tensor(self.isentropic_stiffness_tensor)

    @material_property
    def isentropic_compliance_tensor(self):
        return np.linalg.inv(self.isentropic_stiffness_tensor)

    @material_property
    def full_isentropic_compliance_tensor(self):
        return voigt_notation_to_compliance_tensor(self.isentropic_compliance_tensor)

    @material_property
    def density(self):
        return self._rho

    @material_property
    def isentropic_bulk_modulus_voigt(self):
        """
        :returns: The Voigt bound on the isentropic bulk modulus [Pa].
        :rtype: float
        """
        K = (
            np.sum(
                [
                    [self.isentropic_stiffness_tensor[i][k] for k in range(3)]
                    for i in range(3)
                ]
            )
            / 9.0
        )
        return K

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        :returns: The Reuss bound on the isentropic bulk modulus [Pa].
        :rtype: float
        """
        beta = np.sum(
            [
                [self.isentropic_compliance_tensor[i][k] for k in range(3)]
                for i in range(3)
            ]
        )
        return 1.0 / beta

    @material_property
    def isentropic_bulk_modulus_vrh(self):
        """
        :returns: The Voigt-Reuss-Hill average of the isentropic bulk modulus [Pa].
        :rtype: float
        """
        return 0.5 * (
            self.isentropic_bulk_modulus_voigt + self.isentropic_bulk_modulus_reuss
        )

    @material_property
    def isentropic_shear_modulus_voigt(self):
        """
        :returns: The Voigt bound on the isentropic shear modulus [Pa].
        :rtype: float
        """
        G = (
            np.sum([self.isentropic_stiffness_tensor[i][i] for i in [0, 1, 2]])
            + np.sum([self.isentropic_stiffness_tensor[i][i] for i in [3, 4, 5]]) * 3.0
            - (
                self.isentropic_stiffness_tensor[0][1]
                + self.isentropic_stiffness_tensor[1][2]
                + self.isentropic_stiffness_tensor[2][0]
            )
        ) / 15.0
        return G

    @material_property
    def isentropic_shear_modulus_reuss(self):
        """
        :returns: The Reuss bound on the isentropic shear modulus [Pa].
        :rtype: float
        """
        beta = (
            np.sum([self.isentropic_compliance_tensor[i][i] for i in [0, 1, 2]]) * 4.0
            + np.sum([self.isentropic_compliance_tensor[i][i] for i in [3, 4, 5]]) * 3.0
            - (
                self.isentropic_compliance_tensor[0][1]
                + self.isentropic_compliance_tensor[1][2]
                + self.isentropic_compliance_tensor[2][0]
            )
            * 4.0
        ) / 15.0
        return 1.0 / beta

    @material_property
    def isentropic_shear_modulus_vrh(self):
        """
        :returns: The Voigt-Reuss-Hill average of the isentropic shear modulus [Pa].
        :rtype: float
        """
        return 0.5 * (
            self.isentropic_shear_modulus_voigt + self.isentropic_shear_modulus_reuss
        )

    @material_property
    def isentropic_universal_elastic_anisotropy(self):
        """
        :returns: The universal elastic anisotropy [unitless]
        :rtype: float
        """
        return (
            5.0
            * (
                self.isentropic_shear_modulus_voigt
                / self.isentropic_shear_modulus_reuss
            )
            + (self.isentropic_bulk_modulus_voigt / self.isentropic_bulk_modulus_reuss)
            - 6.0
        )

    @material_property
    def isentropic_isotropic_poisson_ratio(self):
        """
        :returns: The isotropic Poisson ratio (mu) [unitless].
            A metric of the lateral response to loading.
        :rtype: float
        """
        return (
            3.0 * self.isentropic_bulk_modulus_vrh
            - 2.0 * self.isentropic_shear_modulus_vrh
        ) / (
            6.0 * self.isentropic_bulk_modulus_vrh
            + 2.0 * self.isentropic_shear_modulus_vrh
        )

    def christoffel_tensor(self, propagation_direction):
        """
        :returns: The Christoffel tensor from an elastic stiffness
            tensor and a propagation direction for a seismic wave
            relative to the stiffness tensor:
            T_ik = C_ijkl n_j n_l.
        :rtype: float
        """
        propagation_direction = unit_normalize(propagation_direction)
        Tik = np.einsum(
            "ijkl, ...j, ...l",
            self.full_isentropic_stiffness_tensor,
            propagation_direction,
            propagation_direction,
        )
        return Tik

    def isentropic_linear_compressibility(self, direction):
        """
        :returns: The linear isentropic compressibility in a given direction
        relative to the stiffness tensor [1/Pa].
        :rtype: float
        """
        direction = unit_normalize(direction)
        beta = np.einsum(
            "ijkk, ...i, ...j",
            self.full_isentropic_compliance_tensor,
            direction,
            direction,
        )
        return beta

    def isentropic_youngs_modulus(self, direction):
        """
        :returns: The isentropic Youngs modulus in a given direction
        relative to the stiffness tensor [Pa].
        :rtype: float
        """
        direction = unit_normalize(direction)
        S = np.einsum(
            "ijkl, ...i, ...j, ...k, ...l",
            self.full_isentropic_compliance_tensor,
            direction,
            direction,
            direction,
            direction,
        )
        return 1.0 / S

    def isentropic_shear_modulus(self, plane_normal, shear_direction):
        """
        :returns: The isentropic shear modulus on a plane in a given
        shear direction relative to the stiffness tensor [Pa].
        :rtype: float
        """
        plane_normal = unit_normalize(plane_normal)
        shear_direction = unit_normalize(shear_direction)

        assert (
            np.abs(plane_normal.dot(shear_direction)) < np.finfo(float).eps
        ), "plane_normal and shear_direction must be orthogonal"

        G = np.einsum(
            "ijkl, ...i, ...j, ...k, ...l",
            self.full_isentropic_compliance_tensor,
            shear_direction,
            plane_normal,
            shear_direction,
            plane_normal,
        )

        return 0.25 / G

    def isentropic_poissons_ratio(self, axial_direction, lateral_direction):
        """
        :returns: The isentropic poisson ratio given loading and response
        directions relative to the stiffness tensor [unitless].
        :rtype: float
        """

        axial_direction = unit_normalize(axial_direction)
        lateral_direction = unit_normalize(lateral_direction)

        if len(axial_direction.shape) == 1:
            assert (
                np.abs(np.einsum("i, i", axial_direction, lateral_direction)) < 1.0e-12
            ), "axial_direction and lateral_direction must be orthogonal"
        else:
            dot = np.einsum("...i, ...i->...", axial_direction, lateral_direction)
            assert np.all(
                np.abs(dot) < 1.0e-12
            ), "axial_direction and lateral_direction must be orthogonal"

        nu = -(
            self.isentropic_youngs_modulus(axial_direction)
            * np.einsum(
                "ijkl, ...i, ...j, ...k, ...l",
                self.full_isentropic_compliance_tensor,
                lateral_direction,
                lateral_direction,
                axial_direction,
                axial_direction,
            )
        )
        return nu

    def wave_velocities(self, propagation_direction):
        """
        :returns: The compressional wave velocity, and two
        shear wave velocities in a given propagation direction [m/s].
        :rtype: list, containing the wave speeds and directions
            of particle motion relative to the stiffness tensor
        """
        Tik = self.christoffel_tensor(propagation_direction)

        eigenvalues, eigenvectors = np.linalg.eig(Tik)

        # sort eigs by decreasing eigenvalue
        idx = np.flip(eigenvalues.argsort(axis=-1), axis=-1)
        vidx = np.expand_dims(idx, -1)
        eigenvalues = np.take_along_axis(eigenvalues, idx, axis=-1)
        eigenvectors = np.take_along_axis(eigenvectors, vidx, axis=-2)
        eigenvalues = np.real(eigenvalues)
        velocities = np.sqrt(eigenvalues / self.density)

        return velocities, eigenvectors


class IsotropicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C12, C44] (i.e. lambda and mu, the Lame parameters)
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 2
        cijs = list(cijs)
        cijs.insert(0, cijs[0] + 2.0 * cijs[1])  # C11 = C12 + 2C44
        index_lists = [
            [(0, 0), (1, 1), (2, 2)],  # C11
            [(0, 1), (0, 2), (1, 2)],  # C12
            [(3, 3), (4, 4), (5, 5)],
        ]  # C44

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class CubicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C44]
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 3
        index_lists = [
            [(0, 0), (1, 1), (2, 2)],  # C11
            [(0, 1), (0, 2), (1, 2)],  # C12
            [(3, 3), (4, 4), (5, 5)],
        ]  # C44

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class HexagonalMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C33, C44]
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 5
        cijs = list(cijs)
        cijs.append((cijs[0] - cijs[1]) / 2.0)  # C66 = (C11-C12)/2.

        index_lists = [
            [(0, 0), (1, 1)],  # C11
            [(0, 1)],  # C12
            [(0, 2), (1, 2)],  # C13
            [(2, 2)],  # C33
            [(3, 3), (4, 4)],  # C44
            [(5, 5)],
        ]  # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class TetragonalMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C33, C44, C66] or
    [C11, C12, C13, C16, C33, C44, C66]
    """

    def __init__(self, rho, cijs):
        if len(cijs) == 6:
            # Tetragonal I / Laue class 4/mmm
            index_lists = [
                [(0, 0), (1, 1)],  # C11
                [(0, 1)],  # C12
                [(0, 2), (1, 2)],  # C13
                [(2, 2)],  # C33
                [(3, 3), (4, 4)],  # C44
                [(5, 5)],
            ]  # C66
        elif len(cijs) == 7:
            # Tetragonal II / Laue class 4/m
            cijs = list(cijs)
            cijs.insert(4, -cijs[3])  # C26 = -C16
            index_lists = [
                [(0, 0), (1, 1)],  # C11
                [(0, 1)],  # C12
                [(0, 2), (1, 2)],  # C13
                [(0, 5)],  # C16
                [(1, 5)],  # C26
                [(2, 2)],  # C33
                [(3, 3), (4, 4)],  # C44
                [(5, 5)],
            ]  # C66
        else:
            raise Exception(
                "Tetragonal materials should have " "either 6 or 7 independent Cijs"
            )

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class RhombohedralMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C14, C33, C44, C66] or
    [C11, C12, C13, C14, C15, C33, C44, C66]
    """

    def __init__(self, rho, cijs):
        cijs = list(cijs)
        if len(cijs) == 7:
            # Rhombohedral I / Laue class \bar{3}m
            cijs.insert(4, -cijs[3])  # C24 = -C14
            index_lists = [
                [(0, 0), (1, 1)],  # C11
                [(0, 1)],  # C12
                [(0, 2), (1, 2)],  # C13
                [(0, 3), (4, 5)],  # C14
                [(1, 3)],  # C24
                [(2, 2)],  # C33
                [(3, 3), (4, 4)],  # C44
                [(5, 5)],
            ]  # C66

        elif len(cijs) == 8:
            # Rhombohedral II / Laue class \bar{3}
            cijs.insert(4, -cijs[3])  # C24 = -C14
            cijs.insert(6, -cijs[5])  # C25 = -C15
            index_lists = [
                [(0, 0), (1, 1)],  # C11
                [(0, 1)],  # C12
                [(0, 2), (1, 2)],  # C13
                [(0, 3), (4, 5)],  # C14
                [(1, 3)],  # C24
                [(0, 4)],  # C15
                [(1, 4), (3, 5)],  # C25
                [(2, 2)],  # C33
                [(3, 3), (4, 4)],  # C44
                [(5, 5)],
            ]  # C66

        else:
            raise Exception(
                "Rhombohedral materials should have " "either 7 or 8 independent Cijs"
            )

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class OrthorhombicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C22, C23, C33, C44, C55, C66]
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 9
        index_lists = [
            [(0, 0)],  # C11
            [(0, 1)],  # C12
            [(0, 2)],  # C13
            [(1, 1)],  # C22
            [(1, 2)],  # C23
            [(2, 2)],  # C33
            [(3, 3)],  # C44
            [(4, 4)],  # C55
            [(5, 5)],
        ]  # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class MonoclinicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66]
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 13
        index_lists = [
            [(0, 0)],  # C11
            [(0, 1)],  # C12
            [(0, 2)],  # C13
            [(0, 4)],  # C15
            [(1, 1)],  # C22
            [(1, 2)],  # C23
            [(1, 4)],  # C25
            [(2, 2)],  # C33
            [(2, 4)],  # C35
            [(3, 3)],  # C44
            [(3, 5)],  # C46
            [(4, 4)],  # C55
            [(5, 5)],
        ]  # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )


class TriclinicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [Cij, where 1<=i<=6 and i<=j<=6]
    """

    def __init__(self, rho, cijs):
        assert len(cijs) == 21
        index_lists = [[(i, j)] for i in range(6) for j in range(i, 6)]

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists)
        )
