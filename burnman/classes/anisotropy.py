# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np

from ..utils.math import unit_normalize
from .material import Material, material_property

try:  # numpy.block was new in numpy version 1.13.0.
    block = np.block([[np.ones((3, 3)), 2.*np.ones((3, 3))],
                      [2.*np.ones((3, 3)), 4.*np.ones((3, 3))]])
except:
    block = np.array(np.bmat([[[[1.]*3]*3, [[2.]*3]*3],
                              [[[2.]*3]*3, [[4.]*3]*3]] ))
voigt_compliance_factors = block


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

        assert cijs.shape == (6, 6), 'cijs must be in Voigt notation (6x6)'
        assert np.allclose(cijs.T, cijs), 'stiffness_tensor must be symmetric'

        Material.__init__(self)

    def _voigt_index_to_ij(self, m):
        """
        Returns the ij (or kl) indices of the
        stiffness tensor which correspond to those
        of the Voigt notation m (or n).
        """
        if m == 3:
            return 1, 2
        elif m == 4:
            return 0, 2
        elif m == 5:
            return 0, 1
        else:
            return m, m

    def _voigt_notation_to_stiffness_tensor(self, voigt_notation):
        """
        Converts a stiffness tensor in Voigt notation (6x6 matrix)
        to the full fourth rank tensor (3x3x3x3 matrix).
        """
        stiffness_tensor = np.zeros([3, 3, 3, 3])
        for m in range(6):
            i, j = self._voigt_index_to_ij(m)
            for n in range(6):
                k, l = self._voigt_index_to_ij(n)
                stiffness_tensor[i][j][k][l] = voigt_notation[m][n]
                stiffness_tensor[j][i][k][l] = voigt_notation[m][n]
                stiffness_tensor[i][j][l][k] = voigt_notation[m][n]
                stiffness_tensor[j][i][l][k] = voigt_notation[m][n]
        return stiffness_tensor

    def _voigt_notation_to_compliance_tensor(self, voigt_notation):
        return self._voigt_notation_to_stiffness_tensor(np.divide(voigt_notation,
                                                                  voigt_compliance_factors))

    @material_property
    def isentropic_stiffness_tensor(self):
        return self._isentropic_stiffness_tensor

    @material_property
    def full_isentropic_stiffness_tensor(self):
        return self._voigt_notation_to_stiffness_tensor(self.isentropic_stiffness_tensor)

    @material_property
    def isentropic_compliance_tensor(self):
        return np.linalg.inv(self.isentropic_stiffness_tensor)

    @material_property
    def full_isentropic_compliance_tensor(self):
        return self._voigt_notation_to_compliance_tensor(self.isentropic_compliance_tensor)

    @material_property
    def density(self):
        return self._rho

    @material_property
    def isentropic_bulk_modulus_voigt(self):
        """
        Computes the isentropic bulk modulus (Voigt bound)
        """
        K = np.sum([[self.isentropic_stiffness_tensor[i][k]
                     for k in range(3)]
                    for i in range(3)])/9.
        return K

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Computes the isentropic bulk modulus (Reuss bound)
        """
        beta = np.sum([[self.isentropic_compliance_tensor[i][k] for k in range(3)] for i in range(3)])
        return 1./beta

    @material_property
    def isentropic_bulk_modulus_vrh(self):
        """
        Computes the isentropic bulk modulus (Voigt-Reuss-Hill average)
        """
        return 0.5*(self.isentropic_bulk_modulus_voigt + self.isentropic_bulk_modulus_reuss)

    @material_property
    def isentropic_shear_modulus_voigt(self):
        """
        Computes the isentropic shear modulus (Voigt bound)
        """
        G = ( np.sum([self.isentropic_stiffness_tensor[i][i] for i in [0, 1, 2]]) +
              np.sum([self.isentropic_stiffness_tensor[i][i] for i in [3, 4, 5]])*3. -
              ( self.isentropic_stiffness_tensor[0][1] +
                self.isentropic_stiffness_tensor[1][2] +
                self.isentropic_stiffness_tensor[2][0] )) / 15.
        return G

    @material_property
    def isentropic_shear_modulus_reuss(self):
        """
        Computes the isentropic shear modulus (Reuss bound)
        """
        beta =  ( np.sum([self.isentropic_compliance_tensor[i][i] for i in [0, 1, 2]])*4. +
                  np.sum([self.isentropic_compliance_tensor[i][i] for i in [3, 4, 5]])*3. -
                  ( self.isentropic_compliance_tensor[0][1] +
                    self.isentropic_compliance_tensor[1][2] +
                    self.isentropic_compliance_tensor[2][0])*4. ) / 15.
        return 1./beta

    @material_property
    def isentropic_shear_modulus_vrh(self):
        """
        Computes the shear modulus (Voigt-Reuss-Hill average)
        """
        return 0.5*(self.isentropic_shear_modulus_voigt
                    + self.isentropic_shear_modulus_reuss)

    @material_property
    def isentropic_universal_elastic_anisotropy(self):
        """
        Compute the universal elastic anisotropy
        """
        return ( 5.*(self.isentropic_shear_modulus_voigt/self.isentropic_shear_modulus_reuss) +
                 (self.isentropic_bulk_modulus_voigt/self.isentropic_bulk_modulus_reuss) - 6. )

    @material_property
    def isentropic_isotropic_poisson_ratio(self):
        """
        Compute mu, the isotropic Poisson ratio
        (a description of the laterial response to loading)
        """
        return ((3.*self.isentropic_bulk_modulus_vrh
                 - 2.*self.isentropic_shear_modulus_vrh)
                / (6.*self.isentropic_bulk_modulus_vrh
                   + 2.*self.isentropic_shear_modulus_vrh) )

    def christoffel_tensor(self, propagation_direction):
        """
        Computes the Christoffel tensor from an elastic stiffness
        tensor and a propagation direction for a seismic wave
        relative to the stiffness tensor

        T_ik = C_ijkl n_j n_l
        """
        propagation_direction = unit_normalize(propagation_direction)
        Tik = np.tensordot(np.tensordot(self.full_isentropic_stiffness_tensor,
                                        propagation_direction,
                                        axes=([1],[0])),
                           propagation_direction,
                           axes=([2],[0]))
        return Tik

    def isentropic_linear_compressibility(self, direction):
        """
        Computes the linear isentropic compressibility in a given direction
        relative to the stiffness tensor
        """
        direction = unit_normalize(direction)
        Sijkk = np.einsum('ijkk', self.full_isentropic_compliance_tensor)
        beta = Sijkk.dot(direction).dot(direction)
        return beta

    def isentropic_youngs_modulus(self, direction):
        """
        Computes the isentropic Youngs modulus in a given direction
        relative to the stiffness tensor
        """
        direction = unit_normalize(direction)
        Sijkl = self.full_isentropic_compliance_tensor
        S = Sijkl.dot(direction).dot(direction).dot(direction).dot(direction)
        return 1./S

    def isentropic_shear_modulus(self, plane_normal, shear_direction):
        """
        Computes the isentropic shear modulus on a plane in a given
        shear direction relative to the stiffness tensor
        """
        plane_normal = unit_normalize(plane_normal)
        shear_direction = unit_normalize(shear_direction)

        assert np.abs(plane_normal.dot(shear_direction)) < np.finfo(float).eps, 'plane_normal and shear_direction must be orthogonal'
        Sijkl = self.full_isentropic_compliance_tensor
        G = Sijkl.dot(shear_direction).dot(plane_normal).dot(shear_direction).dot(plane_normal)
        return 0.25/G

    def isentropic_poissons_ratio(self,
                                  axial_direction,
                                  lateral_direction):
        """
        Computes the isentropic poisson ratio given loading and response
        directions relative to the stiffness tensor
        """

        axial_direction = unit_normalize(axial_direction)
        lateral_direction = unit_normalize(lateral_direction)
        assert np.abs(axial_direction.dot(lateral_direction)) < np.finfo(float).eps, 'axial_direction and lateral_direction must be orthogonal'

        Sijkl = self.full_isentropic_compliance_tensor
        x = axial_direction
        y = lateral_direction
        nu = -(Sijkl.dot(y).dot(y).dot(x).dot(x) /
               Sijkl.dot(x).dot(x).dot(x).dot(x) )
        return nu

    def wave_velocities(self, propagation_direction):
        """
        Computes the compressional wave velocity, and two
        shear wave velocities in a given propagation direction

        Returns two lists, containing the wave speeds and
        directions of particle motion relative to the stiffness tensor
        """
        propagation_direction = unit_normalize(propagation_direction)

        Tik = self.christoffel_tensor(propagation_direction)

        eigenvalues, eigenvectors = np.linalg.eig(Tik)

        idx = eigenvalues.argsort()[::-1]
        eigenvalues = np.real(eigenvalues[idx])
        eigenvectors = eigenvectors[:,idx]
        velocities = np.sqrt(eigenvalues/self.density)

        return velocities, eigenvectors

def voigt_array_from_cijs(cijs, index_lists):
    """
    Takes a list of cijs and a list of list of tuples corresponding to
    the positions of each cij in the Voigt form matrix.
    Note that the indices run from 0--5, not 1--6.
    """
    C = np.zeros([6, 6])
    for i, index_list in enumerate(index_lists):
        for indices in index_list:
            C[indices] = cijs[i]
            C[indices[::-1]] = cijs[i]
    return C

class IsotropicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C12, C44] (i.e. lambda and mu, the Lame parameters)
    """
    def __init__(self, rho, cijs):

        assert len(cijs) == 2
        cijs = list(cijs)
        cijs.insert(0, cijs[0] + 2.*cijs[1]) # C11 = C12 + 2C44
        index_lists = [[(0, 0), (1, 1), (2, 2)], # C11
                       [(0, 1), (0, 2), (1, 2)], # C12
                       [(3, 3), (4, 4), (5, 5)]] # C44

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))

class CubicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C44]
    """
    def __init__(self, rho, cijs):

        assert len(cijs) == 3
        index_lists = [[(0, 0), (1, 1), (2, 2)], # C11
                       [(0, 1), (0, 2), (1, 2)], # C12
                       [(3, 3), (4, 4), (5, 5)]] # C44

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))

class HexagonalMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C33, C44]
    """
    def __init__(self, rho, cijs):
        assert len(cijs) == 5
        cijs = list(cijs)
        cijs.append((cijs[0] - cijs[1])/2.) # C66 = (C11-C12)/2.

        index_lists = [[(0, 0), (1, 1)], # C11
                       [(0, 1)], # C12
                       [(0, 2), (1, 2)], # C13
                       [(2, 2)], # C33
                       [(3, 3), (4, 4)], # C44
                       [(5, 5)]] # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))

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
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
        elif len(cijs) == 7:
            # Tetragonal II / Laue class 4/m
            cijs = list(cijs)
            cijs.insert(4, -cijs[3]) # C26 = -C16
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 5)], # C16
                           [(1, 5)], # C26
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66
        else:
            raise Exception('Tetragonal materials should have '
                            'either 6 or 7 independent Cijs')

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))


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
            cijs.insert(4, -cijs[3]) # C24 = -C14
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 3), (4, 5)], # C14
                           [(1, 3)], # C24
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66

        elif len(cijs) == 8:
            # Rhombohedral II / Laue class \bar{3}
            cijs.insert(4, -cijs[3]) # C24 = -C14
            cijs.insert(6, -cijs[5]) # C25 = -C15
            index_lists = [[(0, 0), (1, 1)], # C11
                           [(0, 1)], # C12
                           [(0, 2), (1, 2)], # C13
                           [(0, 3), (4, 5)], # C14
                           [(1, 3)], # C24
                           [(0, 4)], # C15
                           [(1, 4), (3, 5)], # C25
                           [(2, 2)], # C33
                           [(3, 3), (4, 4)], # C44
                           [(5, 5)]] # C66

        else:
            raise Exception('Rhombohedral materials should have '
                            'either 7 or 8 independent Cijs')

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))

class OrthorhombicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C22, C23, C33, C44, C55, C66]
    """
    def __init__(self, rho, cijs):

        assert len(cijs) == 9
        index_lists = [[(0, 0)], # C11
                       [(0, 1)], # C12
                       [(0, 2)], # C13
                       [(1, 1)], # C22
                       [(1, 2)], # C23
                       [(2, 2)], # C33
                       [(3, 3)], # C44
                       [(4, 4)], # C55
                       [(5, 5)]] # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))


class MonoclinicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [C11, C12, C13, C15, C22, C23, C25, C33, C35, C44, C46, C55, C66]
    """
    def __init__(self, rho, cijs):

        assert len(cijs) == 13
        index_lists = [[(0, 0)], # C11
                       [(0, 1)], # C12
                       [(0, 2)], # C13
                       [(0, 4)], # C15
                       [(1, 1)], # C22
                       [(1, 2)], # C23
                       [(1, 4)], # C25
                       [(2, 2)], # C33
                       [(2, 4)], # C35
                       [(3, 3)], # C44
                       [(3, 5)], # C46
                       [(4, 4)], # C55
                       [(5, 5)]] # C66

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))

class TriclinicMaterial(AnisotropicMaterial):
    """
    A class derived from the AnisotropicMaterial base class
    Initialization takes two input parameters; rho and
    [Cij, where 1<=i<=6 and i<=j<=6]
    """
    def __init__(self, rho, cijs):

        assert len(cijs) == 21
        index_lists=[[(i, j)] for i in range(6) for j in range(i, 6)]

        AnisotropicMaterial.__init__(
            self, rho, voigt_array_from_cijs(cijs, index_lists))
