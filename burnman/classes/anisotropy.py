# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np

from ..utils.math import unit_normalize
from .material import Material, material_property
from ..utils.anisotropy import (
    voigt_array_from_cijs,
    voigt_notation_to_compliance_tensor,
    voigt_notation_to_stiffness_tensor,
)


def _unit_direction(direction):
    direction = np.asarray(direction)
    assert (
        direction.ndim > 0 and direction.size > 0 and direction.shape[-1] == 3
    ), "Directions must have shape (..., 3) and be nonempty"
    assert np.isrealobj(direction), "Directions must be real"
    return unit_normalize(direction)


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
        cijs = np.asarray(cijs)
        assert (
            np.ndim(rho) == 0 and np.isrealobj(rho) and np.isfinite(rho) and rho > 0.0
        ), "density must be positive and finite"
        assert cijs.shape == (6, 6), "cijs must be in Voigt notation (6x6)"
        assert np.isrealobj(cijs) and np.all(
            np.isfinite(cijs)
        ), "stiffness_tensor must be real and finite"
        assert np.array_equal(cijs.T, cijs), "stiffness_tensor must be symmetric"
        self._isentropic_stiffness_tensor = cijs
        self._rho = rho

        Material.__init__(self)

    def assert_stable(self):
        """Assert positive-definite symmetric isentropic stiffness.

        This is an opt-in elastic stability check for unstressed materials.
        Nonhydrostatic states require their stress-dependent acoustic criteria.
        """
        cijs = self.isentropic_stiffness_tensor
        assert np.isrealobj(cijs) and np.all(
            np.isfinite(cijs)
        ), "stiffness_tensor must be real and finite"
        tolerance = 100.0 * np.finfo(float).eps * np.max(np.abs(cijs))
        assert np.allclose(
            cijs, cijs.T, rtol=0.0, atol=tolerance
        ), "stiffness_tensor must be symmetric to test stability"
        assert np.all(
            np.linalg.eigvalsh(cijs) > 0.0
        ), "stiffness_tensor must be positive definite for elastic stability"

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
        Calculates the adiabatic Christoffel tensor.

        The isentropic stiffness tensors provided in this and other
        anisotropic classes correspond to the isentropic derivative
        of stress with respect to symmetric infinitesimal strain.
        This corresponds to the isentropic elastic tensor B
        of Wallace (1967; doi:10.1103/PhysRev.162.776).

        The Christoffel (acoustic) matrix is the contraction
        T_ik = S_ijkl n_j n_l, where S is Wallace's
        isentropic equation-of-motion tensor
        (not to be confused with the isentropic compliance tensor).
        For a unit direction n, contracting S gives
        T = sym(T_B) + ((n . sigma . n) I - sigma) / 2,
        where T_B_ik = B_ijkl n_j n_l and sigma is the tension-positive
        Cauchy stress. This derivation follows from Wallace (1967),
        using his Equations (2.24), (2.28), and (2.36).

        :param propagation_direction: Direction or array of directions with
            shape (..., 3), expressed in the same frame as the stress and
            stiffness tensors. Each direction is normalized before contraction.
        :returns: Christoffel tensor in Pa, with shape (..., 3, 3).
        :rtype: numpy.ndarray
        """
        n = _unit_direction(propagation_direction)
        Q_B = np.einsum(
            "ijkl,...j,...l->...ik", self.full_isentropic_stiffness_tensor, n, n
        )
        try:
            stress = self.cauchy_stress
            normal_stress = np.einsum("ij,...i,...j->...", stress, n, n)
            correction = 0.5 * (
                normal_stress[..., None, None] * np.eye(3)
                - np.sum(n * n, axis=-1)[..., None, None] * stress
            )
            Tik = 0.5 * (Q_B + Q_B.swapaxes(-1, -2)) + correction
        except AttributeError:
            # If the subclass does not define a cauchy_stress property,
            # assume hydrostatic stress
            # (and therefore a correction equal to zero).
            Tik = Q_B
        return Tik

    def isentropic_linear_compressibility(self, direction):
        """
        :returns: The linear isentropic compressibility in a given direction
            relative to the stiffness tensor [1/Pa].
        :rtype: float
        """
        direction = _unit_direction(direction)
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
        direction = _unit_direction(direction)
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
        plane_normal = _unit_direction(plane_normal)
        shear_direction = _unit_direction(shear_direction)

        assert np.all(
            np.abs(np.einsum("...i,...i->...", plane_normal, shear_direction)) < 1.0e-12
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

        axial_direction = _unit_direction(axial_direction)
        lateral_direction = _unit_direction(lateral_direction)

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

        Polarizations are columns, ordered with the decreasing wave speeds.
        Negative acoustic eigenvalues raise AssertionError.
        """
        Tik = self.christoffel_tensor(propagation_direction)

        assert np.isrealobj(Tik) and np.all(
            np.isfinite(Tik)
        ), "Christoffel tensor must be real and finite"
        tolerance = 100.0 * np.finfo(float).eps * np.max(np.abs(Tik))
        assert np.allclose(
            Tik, Tik.swapaxes(-1, -2), rtol=0.0, atol=tolerance
        ), "Christoffel tensor must be symmetric"
        eigenvalues, eigenvectors = np.linalg.eigh(Tik)
        assert np.all(
            eigenvalues >= 0.0
        ), "Christoffel tensor has negative eigenvalues; real wave speeds are undefined"
        assert (
            np.isfinite(self.density) and self.density > 0.0
        ), "density must be positive and finite"
        velocities = np.sqrt(eigenvalues[..., ::-1] / self.density)

        return velocities, eigenvectors[..., :, ::-1]


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

    Only six or seven coefficients are independent: C66 must equal
    (C11 - C12)/2. Inconsistent supplied values raise AssertionError.
    """

    def __init__(self, rho, cijs):
        cijs = list(cijs)
        assert len(cijs) in (7, 8), "Rhombohedral materials require 7 or 8 Cijs"
        assert np.isclose(
            cijs[-1], (cijs[0] - cijs[1]) / 2.0, rtol=1.0e-12, atol=0.0
        ), "Rhombohedral C66 must equal (C11 - C12)/2"
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
