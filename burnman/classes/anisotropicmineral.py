# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.
import numpy as np
from scipy.linalg import expm
from numpy.linalg import cond
from .mineral import Mineral
from .material import Material, material_property
from .anisotropy import AnisotropicMaterial, voigt_compliance_factors
from ..utils.misc import copy_documentation
from ..utils.unitcell import cell_parameters_to_vectors
from ..utils.unitcell import cell_vectors_to_parameters


class AnisotropicMineral(Mineral, AnisotropicMaterial):
    """
    A class implementing the anisotropic mineral equation of state described
    in :cite:`Myhill2022`.
    This class is derived from both Mineral and AnisotropicMaterial,
    and inherits most of the methods from these classes.

    Instantiation of an AnisotropicMineral takes three required arguments;
    a reference Mineral (i.e. a standard isotropic mineral which provides
    volume as a function of pressure and temperature), cell_parameters,
    which give the lengths of the molar cell vectors and the angles between
    them (see :func:`~burnman.utils.unitcell.cell_parameters_to_vectors`),
    and an anisotropic parameters object, which should be either a
    4D array of anisotropic parameters or a dictionary of parameters which
    describe the anisotropic behaviour of the mineral.
    For a description of the physical meaning of the parameters in the
    4D array, please refer to the code
    or to the original paper.

    If the user chooses to define their parameters as a dictionary,
    they must also provide a function to the psi_function argument
    that describes how to compute the tensors Psi, dPsidf and dPsidPth
    (in Voigt form). The function arguments should be f, Pth and params,
    in that order. The output variables Psi, dPsidf and dPsidth
    must be returned in that order in a tuple. The user should
    also explicitly state whether the material is orthotropic or not
    by supplying a boolean to the orthotropic argument.

    States of the mineral can only be queried after setting the
    pressure and temperature using set_state().

    This class is available as ``burnman.AnisotropicMineral``.

    All the material parameters are expected to be in plain SI units.  This
    means that the elastic moduli should be in Pascals and NOT Gigapascals.
    Additionally, the cell parameters should be in m/(mol formula unit)
    and not in unit cell lengths. To convert unit cell lengths given in
    Angstrom to molar cell parameters you should multiply by 10^(-10) *
    (N_a / Z)^1/3, where N_a is Avogadro's number
    and Z is the number of formula units per unit cell.
    You can look up Z in many places, including www.mindat.org.

    Finally, it is assumed that the unit cell of the anisotropic material
    is aligned in a particular way relative to the coordinate axes
    (the anisotropic_parameters are defined relative to the coordinate axes).
    The crystallographic a-axis is assumed to be parallel to the first
    spatial coordinate axis, and the crystallographic b-axis is assumed to
    be perpendicular to the third spatial coordinate axis.
    """

    def __init__(self, isotropic_mineral, cell_parameters,
                 anisotropic_parameters,
                 psi_function=None, orthotropic=None):

        if psi_function is None:
            self.check_standard_parameters(anisotropic_parameters)
            self.anisotropic_params = {'c': anisotropic_parameters}
            self.psi_function = self.standard_psi_function
        else:
            if not isinstance(orthotropic, bool):
                raise Exception('If the Psi function is provided, '
                                'you must specify whether your material is '
                                'orthotropic as a boolean.')
            self.orthotropic = orthotropic
            self.anisotropic_params = anisotropic_parameters
            self.psi_function = psi_function

        self.cell_vectors_0 = cell_parameters_to_vectors(cell_parameters)

        if (np.abs(np.linalg.det(self.cell_vectors_0)
                   - isotropic_mineral.params['V_0']) > np.finfo(float).eps):
            factor = np.cbrt(isotropic_mineral.params['V_0']
                             / np.linalg.det(self.cell_vectors_0))
            raise Exception('The standard state unit vectors are inconsistent '
                            'with the volume. Suggest multiplying each '
                            f'by {factor}.')

        if 'name' in isotropic_mineral.params:
            self.name = isotropic_mineral.params['name']

        Mineral.__init__(self, isotropic_mineral.params,
                         isotropic_mineral.property_modifiers)

    def standard_psi_function(self, f, Pth, params):
        # Compute Psi, dPsidPth, dPsidf, needed by most anisotropic properties
        c = params['c']
        ns = np.arange(c.shape[-1])
        x = c[:, :, 0, :] + c[:, :, 1, :]*f
        dPsidf = c[:, :, 1, :]

        for i in list(range(2, c.shape[2])):
            # non-intuitively, the += operator doesn't simply add in-place,
            # so here we overwrite the arrays with new ones
            x = x + c[:, :, i, :]*np.power(f, float(i))/float(i)
            dPsidf = dPsidf + c[:, :, i, :] * np.power(f, float(i)-1.)

        Psi = np.einsum('ikn, n->ik', x, np.power(Pth, ns))

        dPsidPth = np.einsum('ikn, n->ik',
                             x[:, :, 1:],
                             ns[1:]*np.power(Pth, ns[1:]-1))

        dPsidf = np.einsum('ikn, n->ik', dPsidf, np.power(Pth, ns))
        return (Psi, dPsidf, dPsidPth)

    @copy_documentation(Material.set_state)
    def set_state(self, pressure, temperature):
        # 1) Compute dPthdf|T
        # relatively large dP needed for accurate estimate of dPthdf
        Mineral.set_state(self, pressure, temperature)
        dP = self.isothermal_bulk_modulus_reuss*1.e-5

        Mineral.set_state(self, pressure-dP/2., temperature)
        V1 = self.V
        Pth1 = pressure-dP/2. - self.method.pressure(self.params['T_0'],
                                                     V1, self.params)

        Mineral.set_state(self, pressure+dP/2., temperature)
        V2 = self.V
        Pth2 = pressure+dP/2. - self.method.pressure(self.params['T_0'],
                                                     V2, self.params)

        self.dPthdf = (Pth2 - Pth1) / np.log(V2/V1)
        Mineral.set_state(self, pressure, temperature)

        # 2) Compute other properties needed for anisotropic equation of state
        V = self.V
        V_0 = self.params['V_0']
        Vrel = V/V_0
        f = np.log(Vrel)
        self._Vrel = Vrel
        self._f = f

        Pth = pressure - self.method.pressure(self.params['T_0'], self.V,
                                              self.params)

        out = self.psi_function(f, Pth, self.anisotropic_params)

        self.Psi_Voigt = out[0]
        self.dPsidf_Voigt = out[1]
        self.dPsidPth_Voigt = out[2]

    def _contract_compliances(self, compliances):
        """
        Takes a compliance tensor in standard (3x3x3x3) form
        and returns the Voigt form (6x6). Note the compliance factors
        which are required to maintain the inverse relationship with the
        corresponding stiffness tensor.
        """
        voigt_notation = np.zeros((6, 6))
        for p in range(6):
            i, j = self._voigt_index_to_ij(p)
            for q in range(6):
                m, n = self._voigt_index_to_ij(q)
                voigt_notation[p, q] = compliances[i, j, m, n]
        return np.multiply(voigt_notation, voigt_compliance_factors)

    @material_property
    def deformation_gradient_tensor(self):
        """
        Returns
        -------
        deformation_gradient_tensor : 2D numpy array
            The deformation gradient tensor describing the deformation of the
            mineral from its undeformed state
            (i.e. the state at the reference pressure and temperature).
        """
        Psi_full = self._voigt_notation_to_compliance_tensor(self.Psi_Voigt)
        F = expm(np.einsum('ijkl, kl', Psi_full, np.eye(3)))
        return F

    @material_property
    def unrotated_cell_vectors(self):
        """
        Returns
        -------
        unrotated_cell_vectors : 2D numpy array
            The vectors of the cell constructed from one mole of formula units
            after deformation of the mineral from its undeformed state
            (i.e. the state at the reference pressure and temperature).
            Each vector is given in [m]. See the documentation for the function
            :func:`~burnman.utils.unitcell.cell_parameters_to_vectors`
            for the assumed relationships between the cell vectors and
            spatial coordinate axes.
        """
        return self.deformation_gradient_tensor.dot(self.cell_vectors_0)

    @material_property
    def deformed_coordinate_frame(self):
        """
        Returns
        -------
        deformed_coordinate_frame : 2D numpy array
            The orientations of the three spatial coordinate axes
            after deformation of the mineral. For orthotropic minerals,
            this is equal to the identity matrix, as hydrostatic stresses only
            induce rotations in monoclinic and triclinic crystals.
        """
        if self.orthotropic:
            return np.eye(3)
        else:
            M = self.unrotated_cell_vectors
            Q = np.empty((3, 3))
            Q[0] = M[0] / np.linalg.norm(M[0])
            Q[2] = np.cross(M[0], M[1])/np.linalg.norm(np.cross(M[0], M[1]))
            Q[1] = np.cross(Q[2], Q[0])
            return Q

    @material_property
    def rotation_matrix(self):
        """
        Returns
        -------
        rotation_matrix : 2D numpy array
            The matrix required to rotate the properties of the deformed
            mineral into the deformed coordinate frame. For orthotropic
            minerals, this is equal to the identity matrix.
        """
        return self.deformed_coordinate_frame.T

    @material_property
    def cell_vectors(self):
        """
        Returns
        -------
        cell_vectors : 2D numpy array
            The vectors of the cell constructed from one mole of formula units.
            Each vector is given in [m]. See the documentation for the function
            :func:`~burnman.utils.unitcell.cell_parameters_to_vectors`
            for the assumed relationships between the cell vectors and
            spatial coordinate axes.
        """
        if self.orthotropic:
            return self.unrotated_cell_vectors
        else:
            return np.einsum('ij, jk->ik',
                             self.unrotated_cell_vectors,
                             self.rotation_matrix)

    @material_property
    def cell_parameters(self):
        """
        Returns
        -------
        cell_parameters : 1D numpy array
            The molar cell parameters of the mineral, given in standard form:
            [:math:`a`, :math:`b`, :math:`c`,
            :math:`\\alpha`, :math:`\\beta`, :math:`\\gamma`],
            where the first three floats are the lengths of the vectors in [m]
            defining the cell constructed from one mole of formula units.
            The last three floats are angles between vectors
            (given in radians). See the documentation for the function
            :func:`~burnman.utils.unitcell.cell_parameters_to_vectors`
            for the assumed relationships between the cell vectors and
            spatial coordinate axes.
        """
        return cell_vectors_to_parameters(self.cell_vectors)

    @material_property
    def shear_modulus(self):
        """
        Anisotropic minerals do not (in general) have a single shear modulus.
        This function returns a NotImplementedError. Users should instead
        consider directly querying the elements in the
        isothermal_stiffness_tensor or isentropic_stiffness_tensor.
        """
        raise NotImplementedError("Anisotropic minerals do not have a shear "
                                  "modulus property. Query "
                                  "the isentropic or isothermal stiffness "
                                  "tensors directory, or use"
                                  "isentropic_shear_modulus_reuss or "
                                  "isentropic_shear_modulus_voigt.")

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Anisotropic minerals do not have a single isothermal bulk modulus.
        This function returns a NotImplementedError. Users should instead
        consider either using isothermal_bulk_modulus_reuss,
        isothermal_bulk_modulus_voigt,
        or directly querying the elements in the isothermal_stiffness_tensor.
        """
        raise NotImplementedError("isothermal_bulk_modulus is not "
                                  "sufficiently explicit for an "
                                  "anisotropic mineral. Did you mean "
                                  "isothermal_bulk_modulus_reuss?")

    @material_property
    def isentropic_bulk_modulus(self):
        """
        Anisotropic minerals do not have a single isentropic bulk modulus.
        This function returns a NotImplementedError. Users should instead
        consider either using isentropic_bulk_modulus_reuss,
        isentropic_bulk_modulus_voigt (both derived from Anisotropicmineral),
        or directly querying the elements in the isentropic_stiffness_tensor.
        """
        raise NotImplementedError("isentropic_bulk_modulus is not "
                                  "sufficiently explicit for an "
                                  "anisotropic mineral. Did you mean "
                                  "isentropic_bulk_modulus_reuss?")

    isothermal_bulk_modulus_reuss = Mineral.isothermal_bulk_modulus

    @material_property
    def isothermal_bulk_modulus_voigt(self):
        """
        Returns
        -------
        isothermal_bulk_modulus_voigt : float
            The Voigt bound on the isothermal bulk modulus in [Pa].
        """
        K = np.sum([[self.isothermal_stiffness_tensor[i][k]
                     for k in range(3)]
                    for i in range(3)])/9.
        return K

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns
        -------
        isothermal_compressibility_reuss : float
            The Reuss bound on the isothermal compressibility in [1/Pa].
        """
        return 1./self.isothermal_bulk_modulus_reuss

    @material_property
    def isothermal_compliance_tensor(self):
        """
        Returns
        -------
        isothermal_compliance_tensor : 2D numpy array
            The isothermal compliance tensor [1/Pa]
            in Voigt form (:math:`\\mathbb{S}_{\\text{T} pq}`).
        """
        S_T = (self.isothermal_compressibility_reuss
               * (self.dPsidf_Voigt + self.dPsidPth_Voigt * self.dPthdf))
        if self.orthotropic:
            return S_T
        else:
            R = self.rotation_matrix
            S = self._voigt_notation_to_compliance_tensor(S_T)
            S_rotated = np.einsum('mi, nj, ok, pl, ijkl->mnop', R, R, R, R, S)
            return self._contract_compliances(S_rotated)

    @material_property
    def thermal_expansivity_tensor(self):
        """
        Returns
        -------
        thermal_expansivity_tensor : 2D numpy array
            The tensor of thermal expansivities [1/K].
        """
        a = self.alpha * (self.dPsidf_Voigt
                          + self.dPsidPth_Voigt
                          * (self.dPthdf
                             + 1./self.isothermal_compressibility_reuss))
        alpha = np.einsum('ijkl, kl',
                          self._voigt_notation_to_compliance_tensor(a),
                          np.eye(3))

        if self.orthotropic:
            return alpha
        else:
            R = self.rotation_matrix
            return np.einsum('mi, nj, ij->mn', R, R, alpha)

    # Derived properties start here
    @material_property
    def isothermal_stiffness_tensor(self):
        """
        Returns
        -------
        isothermal_stiffness_tensor : 2D numpy array
            The isothermal stiffness tensor [Pa]
            in Voigt form (:math:`\\mathbb{C}_{\\text{T} pq}`).
        """
        return np.linalg.inv(self.isothermal_compliance_tensor)

    @material_property
    def full_isothermal_compliance_tensor(self):
        """
        Returns
        -------
        full_isothermal_stiffness_tensor : 4D numpy array
            The isothermal compliance tensor [1/Pa]
            in standard form (:math:`\\mathbb{S}_{\\text{T} ijkl}`).
        """
        S_Voigt = self.isothermal_compliance_tensor
        return self._voigt_notation_to_compliance_tensor(S_Voigt)

    @material_property
    def full_isothermal_stiffness_tensor(self):
        """
        Returns
        -------
        full_isothermal_stiffness_tensor : 4D numpy array
            The isothermal stiffness tensor [Pa]
            in standard form (:math:`\\mathbb{C}_{\\text{T} ijkl}`).
        """
        CT = self.isothermal_stiffness_tensor
        return self._voigt_notation_to_stiffness_tensor(CT)

    @material_property
    def full_isentropic_compliance_tensor(self):
        """
        Returns
        -------
        full_isentropic_stiffness_tensor : 4D numpy array
            The isentropic compliance tensor [1/Pa]
            in standard form (:math:`\\mathbb{S}_{\\text{N} ijkl}`).
        """
        return (self.full_isothermal_compliance_tensor
                - np.einsum('ij, kl->ijkl',
                            self.thermal_expansivity_tensor,
                            self.thermal_expansivity_tensor)
                * self.V * self.temperature / self.C_p)

    @material_property
    def isentropic_compliance_tensor(self):
        """
        Returns
        -------
        isentropic_compliance_tensor : 2D numpy array
            The isentropic compliance tensor [1/Pa]
            in Voigt form (:math:`\\mathbb{S}_{\\text{N} pq}`).
        """
        S_full = self.full_isentropic_compliance_tensor
        return self._contract_compliances(S_full)

    @material_property
    def isentropic_stiffness_tensor(self):
        """
        Returns
        -------
        isentropic_stiffness_tensor : 2D numpy array
            The isentropic stiffness tensor [Pa]
            in Voigt form (:math:`\\mathbb{C}_{\\text{N} pq}`).
        """
        return np.linalg.inv(self.isentropic_compliance_tensor)

    @material_property
    def full_isentropic_stiffness_tensor(self):
        """
        Returns
        -------
        full_isentropic_stiffness_tensor : 4D numpy array
            The isentropic stiffness tensor [Pa]
            in standard form (:math:`\\mathbb{C}_{\\text{N} ijkl}`).
        """
        C_Voigt = self.isentropic_stiffness_tensor
        return self._voigt_notation_to_stiffness_tensor(C_Voigt)

    @material_property
    def grueneisen_tensor(self):
        """
        Returns
        -------
        grueneisen_tensor : 2D numpy array
            The grueneisen tensor.
            This is defined by :cite:`BarronMunn1967` as
            :math:`\\mathbb{C}_{\\text{N} ijkl} \\alpha_{kl} V/C_{P}`.
        """
        return (np.einsum('ijkl, kl->ij',
                          self.full_isentropic_stiffness_tensor,
                          self.thermal_expansivity_tensor)
                * self.molar_volume / self.molar_heat_capacity_p)

    @material_property
    def grueneisen_parameter(self):
        """
        Anisotropic minerals do not (in general) have a single grueneisen
        parameter. This function returns a NotImplementedError.
        Users should instead consider directly querying the elements in the
        grueneisen_tensor.
        """
        raise NotImplementedError("Anisotropic minerals do not have a single "
                                  "grueneisen parameter. Query "
                                  "the grueneisen_tensor instead.")

    @material_property
    def isothermal_compressibility_tensor(self):
        """
        Returns
        -------
        isothermal_compressibility_tensor : 2D numpy array
            The isothermal compressibility tensor.
        """
        return np.einsum('ijkl, kl->ij',
                         self.full_isothermal_compliance_tensor,
                         np.eye(3))

    @material_property
    def isentropic_compressibility_tensor(self):
        """
        Returns
        -------
        isentropic_compressibility_tensor : 2D numpy array
            The isentropic compressibility tensor.
        """
        return np.einsum('ijkl, kl->ij',
                         self.full_isentropic_compliance_tensor,
                         np.eye(3))

    @material_property
    def thermal_stress_tensor(self):
        """
        Returns
        -------
        thermal stress : 2D numpy array
            The change in stress with temperature at constant strain.
        """
        pi = np.einsum('ijkl, kl', self.full_isothermal_stiffness_tensor,
                       self.thermal_expansivity_tensor)
        return pi

    @material_property
    def molar_isometric_heat_capacity(self):
        """
        Returns
        -------
        molar_isometric_heat_capacity : float
            The molar heat capacity at constant strain.
        """

        pi = self.thermal_stress_tensor
        pipiV = np.einsum('ij, kl -> ijkl', pi, pi)*self.V
        indices = np.where(np.abs(pipiV) > 1.e-5)
        values = ((self.full_isentropic_stiffness_tensor
                   - self.full_isothermal_stiffness_tensor)[indices]
                  / pipiV[indices])
        if not np.allclose(values, np.ones_like(values)*values[0],
                           rtol=1.e-5):
            raise Exception('Could not calculate the molar heat '
                            'capacity at constant strain. '
                            'There is an inconsistency in the '
                            'equation of state.')
        C_isometric = self.temperature/values[0]

        return C_isometric

    def check_standard_parameters(self, anisotropic_parameters):

        if not np.all(anisotropic_parameters[:, :, 0, 0] == 0):
            raise Exception("anisotropic_parameters_pqmn should be set to "
                            "zero for all m = n = 0")

        sum_ijij_block = np.sum(anisotropic_parameters[:3, :3, :, :],
                                axis=(0, 1))

        if np.abs(sum_ijij_block[1, 0] - 1.) > 1.e-5:
            raise Exception('The sum of the upper 3x3 pq-block of '
                            'anisotropic_parameters_pqmn must equal '
                            '1 for m=1, n=0 for consistency with the volume. '
                            f'Value is {sum_ijij_block[1, 0]}')

        for m in range(2, len(sum_ijij_block)):
            if np.abs(sum_ijij_block[m, 0]) > 1.e-10:
                raise Exception('The sum of the upper 3x3 pq-block of '
                                'anisotropic_parameters_pqmn must equal 0 for'
                                f'm={m}, n=0 for consistency with the volume. '
                                f'Value is {sum_ijij_block[m, 0]}')

        for m in range(len(sum_ijij_block)):
            for n in range(1, len(sum_ijij_block[0])):
                if np.abs(sum_ijij_block[m, n]) > 1.e-10:
                    raise Exception('The sum of the upper 3x3 pq-block of '
                                    'anisotropic_parameters_pqmn must equal '
                                    f'0 for m={m}, n={n} for '
                                    'consistency with the volume. '
                                    f'Value is {sum_ijij_block[m, n]}')

        if cond(anisotropic_parameters[:, :, 1, 0]) > 1/np.finfo(float).eps:
            raise Exception('anisotropic_parameters[:, :, 1, 0] is singular')

        sum_lower_left_block = np.sum(anisotropic_parameters[3:, :3, :, :],
                                      axis=1)

        self.orthotropic = True
        for i, s in enumerate(sum_lower_left_block):
            if not np.all(np.abs(s) < 1.e-10):
                self.orthotropic = False
