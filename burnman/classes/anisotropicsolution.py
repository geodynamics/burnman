# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.
import numpy as np
import copy
from scipy.linalg import expm, logm
from scipy.optimize import minimize
from .solution import Solution
from .anisotropicmineral import (
    AnisotropicMineral,
    convert_f_Pth_to_f_T_derivatives,
    deformation_gradient_alpha_and_compliance,
)
from .material import material_property, cached_property, Material
from ..utils.anisotropy import (
    contract_stresses,
    expand_stresses,
    voigt_notation_to_compliance_tensor,
)
from ..utils.misc import copy_documentation


class AnisotropicSolution(Solution, AnisotropicMineral):
    """
    A class implementing the anisotropic solution model described
    in :cite:`Myhill2024a`.
    This class is derived from Solution and AnisotropicMineral,
    and inherits most of the methods from those classes.

    Instantiation of an AnisotropicSolution is similar to that of a scalar
    Solution, except that each of the endmembers must be an instance of
    an AnisotropicMineral, and an additional function and parameter dictionary
    are passed to the constructor of AnisotropicSolution that describe
    excess contributions to the anisotropic state tensor (Psi_xs)
    and its derivatives with respect to volume and temperature.
    The function arguments should be ln(V), Pth,
    p (a vector of proportions) and the parameter dictionary,
    in that order.
    The output variables Psi_xs_Voigt, dPsidf_Pth_Voigt_xs and
    dPsidPth_f_Voigt_xs (all 6x6 matrices) and dPsiIdp_xs
    (a 3 x 3 x n_endmember matrix) must be returned in that
    order in a tuple.
    States of the mineral can only be queried after setting the
    pressure and temperature using set_state() and the composition using
    set_composition().

    This class is available as ``burnman.AnisotropicSolution``.
    """

    def __init__(
        self,
        name,
        solution_model,
        psi_excess_function,
        anisotropic_parameters,
        molar_fractions=None,
    ):
        # Always set orthotropic as false to ensure correct
        # derivatives are taken.
        # To do: relax this to allow faster calculations
        self.orthotropic = False

        self.T_0 = 298.15

        # Initialise as both Material and Solution object
        Material.__init__(self)
        Solution.__init__(self, name, solution_model, molar_fractions)

        # Store a scalar copy of the solution model to speed up set_state
        scalar_solution_model = copy.copy(solution_model)
        scalar_solution_model.endmembers = [
            [mbr.isotropic_mineral, formula] for mbr, formula in self.endmembers
        ]
        self._scalar_solution = Solution(name, scalar_solution_model, molar_fractions)

        self._logm_M0_mbr = np.einsum(
            "kij->ijk", np.array([logm(m[0].cell_vectors_0.T) for m in self.endmembers])
        )
        self.frame_convention = self.endmembers[0][0].frame_convention
        self.anisotropic_params = anisotropic_parameters
        self.psi_excess_function = psi_excess_function

        # Finally, set the composition
        if molar_fractions is not None:
            self.set_composition(molar_fractions)

    def set_state(self, pressure, temperature):
        # Set solution conditions
        ss = self._scalar_solution
        ss.set_state(pressure, temperature)

        try:
            self.compute_base_properties()
        except AttributeError:
            raise Exception("You must use set_composition() before set_state().")

        Material.set_state(self, pressure, temperature)

    def compute_base_properties(self):
        pressure = self._scalar_solution.pressure
        temperature = self._scalar_solution.temperature

        ss = self._scalar_solution
        V = ss.molar_volume
        KT_at_T = ss.isothermal_bulk_modulus_reuss
        f = np.log(V)
        self._f = f
        # Evaluate endmember properties at V, T
        # Here we explicitly manipulate each of the anisotropic endmembers
        pressure_guesses = [np.max([0.0e9, pressure - 2.0e9]), pressure + 2.0e9]
        mbrs = self.endmembers
        for mbr in mbrs:
            mbr[0].set_state_with_volume(V, temperature, pressure_guesses)

        # Endmember cell vectors and other functions of Psi (all unrotated)
        self._PsiI_mbr = np.einsum("kij->ijk", np.array([m[0]._PsiI for m in mbrs]))
        dPsidf_T_Voigt_mbr = np.array([m[0]._dPsidf_T_Voigt for m in mbrs])
        dPsiIdf_T_mbr = np.array([m[0]._dPsiIdf_T for m in mbrs])
        dPsiIdT_f_mbr = np.array([m[0]._dPsiIdT_f for m in mbrs])

        fr = self.molar_fractions
        PsiI_ideal = np.einsum("ijp, p->ij", self._PsiI_mbr, fr)
        dPsidf_T_Voigt_ideal = np.einsum("p,pij->ij", fr, dPsidf_T_Voigt_mbr)
        dPsiIdf_T_ideal = np.einsum("p,pij->ij", fr, dPsiIdf_T_mbr)
        dPsiIdT_f_ideal = np.einsum("p,pij->ij", fr, dPsiIdT_f_mbr)

        # Now calculate the thermal pressure terms by
        # evaluating the scalar solution at V, T_0
        ss.set_state_with_volume(V, self.T_0)
        P_at_T0 = ss.pressure
        KT_at_T0 = ss.isothermal_bulk_modulus_reuss
        self.dPthdf = KT_at_T0 - KT_at_T
        self.Pth = pressure - P_at_T0

        # And we're done with calculating endmember properties at V
        # Now we can set state of this object and the scalar one
        ss.set_state(pressure, temperature)
        Solution.set_state(self, pressure, temperature)

        # Calculate excess properties at the given f and Pth
        out = self.psi_excess_function(
            f, self.Pth, self.molar_fractions, self.anisotropic_params
        )
        Psi_xs_Voigt, dPsidf_Pth_Voigt_xs, dPsidPth_f_Voigt_xs = out[:3]
        self._dPsiIdp_xs = out[3]
        Psi_xs_full = voigt_notation_to_compliance_tensor(Psi_xs_Voigt)
        PsiI_xs = np.einsum("ijkl, kl", Psi_xs_full, np.eye(3))

        # Change of variables: (f, Pth) -> Psi(f, T)
        aK_T = ss.alpha * ss.isothermal_bulk_modulus_reuss
        out = convert_f_Pth_to_f_T_derivatives(
            dPsidf_Pth_Voigt_xs, dPsidPth_f_Voigt_xs, aK_T, self.dPthdf
        )
        dPsidf_T_Voigt_xs, dPsiIdf_T_xs, dPsiIdT_f_xs = out

        self._PsiI = PsiI_ideal + PsiI_xs

        out = deformation_gradient_alpha_and_compliance(
            ss.alpha,
            ss.isothermal_compressibility_reuss,
            self._PsiI,
            dPsidf_T_Voigt_ideal + dPsidf_T_Voigt_xs,
            dPsiIdf_T_ideal + dPsiIdf_T_xs,
            dPsiIdT_f_ideal + dPsiIdT_f_xs,
        )
        (
            self._unrotated_F,
            self._unrotated_dFdf_T,
            self._unrotated_alpha,
            self._unrotated_S_T_Voigt,
        ) = out

    def set_composition(self, molar_fractions):
        self._scalar_solution.set_composition(molar_fractions)
        self.cell_vectors_0 = expm(
            np.einsum("ijk, k ->ji", self._logm_M0_mbr, molar_fractions)
        )
        # Calculate all other required properties
        Solution.set_composition(self, molar_fractions)
        if self._scalar_solution.pressure is not None:
            self.compute_base_properties()

    @cached_property
    def ones(self):
        return np.ones(self.n_endmembers)

    @cached_property
    def eye(self):
        return np.eye(self.n_endmembers)

    @material_property
    def _dM0dp_fixed_VT(self):
        """
        Gradient in the standard state cell tensor
        with respect to molar proportions at constant
        volume and temperature under hydrostatic conditions.
        """
        lnM0_ones = np.einsum("ij, k->kij", logm(self.cell_vectors_0.T), self.ones)
        dp = 1.0e-5

        dlnM0 = np.einsum("ijk->kij", self._logm_M0_mbr) * dp / 2.0
        logmM0_0 = lnM0_ones - dlnM0
        logmM0_1 = lnM0_ones + dlnM0
        return np.einsum("kij->ijk", (expm(logmM0_1) - expm(logmM0_0)) / dp)

    @material_property
    def _unrotated_dFdp_fixed_VT(self):
        """
        Gradient in the unrotated deformation gradient tensor
        with respect to molar proportions at constant
        volume and temperature under hydrostatic conditions.
        """
        dp = 1.0e-5
        PsiI_ones = np.einsum("ij, k->kij", self._PsiI, self.ones)
        dPsiI = np.einsum("ijk->kij", self._PsiI_mbr + self._dPsiIdp_xs) * dp / 2.0
        PsiI_0 = PsiI_ones - dPsiI
        PsiI_1 = PsiI_ones + dPsiI

        return np.einsum("kij->ijk", (expm(PsiI_1) - expm(PsiI_0)) / dp)

    @material_property
    def _unrotated_dMdn_fixed_VT(self):
        """
        Gradient in unrotated cell tensor with respect to composition
        at fixed volume and temperature under hydrostatic conditions.
        """
        F = self._unrotated_F
        M0 = self.cell_vectors_0.T
        fr = self.molar_fractions
        n_dpdn = self.eye - np.einsum("l, m->lm", self.ones, fr)
        n = np.sum(self.molar_fractions)
        n23 = np.power(n, 2.0 / 3.0)

        dM0dn = (
            np.einsum("kj, l->kjl", M0, self.ones / 3.0)
            + np.einsum("kjm,lm->kjl", self._dM0dp_fixed_VT, n_dpdn)
        ) / n23
        n_dVmoldn = -self.molar_volume * self.ones

        dFdV = self._unrotated_dFdf_T / self.molar_volume
        dFdn = (
            np.einsum("ij, l->ijl", dFdV, n_dVmoldn)
            + np.einsum("ikm,lm->ikl", self._unrotated_dFdp_fixed_VT, n_dpdn)
        ) / n

        dFdnM0 = np.einsum("ikl,kj->ijl", dFdn, M0)
        FdM0dn = np.einsum("ik,kjl->ijl", F, dM0dn)
        return dFdnM0 + FdM0dn

    @material_property
    def _dMdn_fixed_VT(self):
        """
        Gradient in cell tensor with respect to composition
        at fixed volume and temperature under hydrostatic conditions.
        """
        R = self.rotation_matrix
        return np.einsum("mi, nj, ijk->mnk", R, R, self._unrotated_dMdn_fixed_VT)

    @material_property
    def depsdn_fixed_VT(self):
        """
        Gradient in strain with respect to composition
        at fixed volume and temperature under hydrostatic conditions.
        """
        invM = np.linalg.inv(self.cell_vectors.T)
        Ln = np.einsum("ijm,kj->ikm", self._dMdn_fixed_VT, invM)
        LnT = np.einsum("ikm->kim", Ln)
        return 0.5 * (Ln + LnT)


depsdeps = np.einsum("pm, qn->pqmn", np.eye(3), np.eye(3))


class RelaxedAnisotropicSolution(AnisotropicSolution):
    """
    A class implementing the relaxed anisotropic solution
    model described in :cite:`Myhill2024b`.
    This class is derived from AnisotropicSolution,
    and inherits most of the methods from that class.

    Instantiation of a RelaxedAnisotropicSolution involves passing
    an AnisotropicSolution, plus a set of vectors that represent rapid
    deformation modes. For example, a solution of MgO, FeHSO and FeLSO
    (high and low spin wuestite) can rapidly change proportion of
    high spin and low spin iron, and so a single vector should be passed:
    np.array([[0., -1., 1.]]) or some multiple thereof.

    States of the mineral can only be queried after setting the
    pressure and temperature using set_state() and the composition using
    set_composition().

    This class is available as ``burnman.RelaxedAnisotropicSolution``.
    """

    def __init__(
        self,
        anisotropic_solution,
        relaxation_vectors,
        unrelaxed_vectors,
    ):
        # Make an attribute with the unrelaxed solution
        self.unrelaxed = anisotropic_solution

        # The relaxation vectors
        self.n_relaxation_vectors = len(relaxation_vectors)
        self.n_unrelaxed_vectors = len(unrelaxed_vectors)

        self.dndq = np.array(relaxation_vectors).T
        assert len(self.dndq.shape) == 2
        assert len(self.dndq) == self.unrelaxed.n_endmembers

        self.dndx = np.array(unrelaxed_vectors).T
        assert len(self.dndx.shape) == 2
        assert len(self.dndx) == self.unrelaxed.n_endmembers
        assert (
            len(unrelaxed_vectors) + len(relaxation_vectors)
            == self.unrelaxed.n_endmembers
        )

        self.q_initial = np.zeros(len(relaxation_vectors)) + 0.001

        try:
            molar_fractions = anisotropic_solution.molar_fractions
        except AttributeError:
            molar_fractions = None

        # Give the relaxed solution the same base properties as the
        # unrelaxed solution
        AnisotropicSolution.__init__(
            self,
            name=anisotropic_solution.name,
            solution_model=anisotropic_solution.solution_model,
            psi_excess_function=anisotropic_solution.psi_excess_function,
            anisotropic_parameters=anisotropic_solution.anisotropic_params,
            molar_fractions=molar_fractions,
        )

    def set_state(self, pressure, temperature, relaxed=True):
        """
        Sets the state of the solution. Also relaxes the
        structure parameters if set_composition() has already
        been used and if the relaxed argument has been set to
        True.

        :param pressure: The pressure of the solution [Pa]
        :type pressure: float
        :param temperature: The temperature of the solution [K]
        :type temperature: float
        :param relaxed: Whether to minimize the Gibbs energy
            of the material by changing the values of the structure
            parameters. Defaults to True.
        :type relaxed: bool, optional
        """
        self.unrelaxed.set_state(pressure, temperature)

        if hasattr(self.unrelaxed, "molar_fractions") and relaxed:
            self._relax_at_PTX()
            AnisotropicSolution.set_state(self, pressure, temperature)

    def set_composition(self, molar_fractions, q_initial=None, relaxed=True):
        """
        Sets the composition of the model. Also relaxes the
        structure parameters if set_state() has already
        been used and if the relaxed argument has been set to
        True.

        :param molar_fractions: Molar fractions of the
            independent endmembers corresponding to the
            unrelaxed vectors specified during initialisation.
        :type molar_fractions: 1D numpy array
        :param q_initial: Initial values of the structure parameters.
            Defaults to None, in which case the preexisting
            initial values are used
            (first set to 0.001 during initialisation).
        :type q_initial: 1D numpy array, optional
        :param relaxed: Whether to minimize the Gibbs energy
            of the material by changing the values of the structure
            parameters. Defaults to True.
        :type relaxed: bool, optional
        """
        self.unrelaxed_vectors = molar_fractions
        if q_initial is not None:
            self.q_initial = np.array(q_initial)
        n = np.einsum("ij, j", self.dndq, self.q_initial) + np.einsum(
            "ij, j", self.dndx, molar_fractions
        )

        self.unrelaxed.set_composition(n)

        if self.unrelaxed.pressure is not None and relaxed:
            self._relax_at_PTX()
            AnisotropicSolution.set_composition(self, n)
            self.unrelaxed._scalar_solution.set_composition(self.molar_fractions)
            self.unrelaxed.set_composition(self.molar_fractions)

    def _relax_at_PTX(self):
        """
        Minimizes the Gibbs energy at constant pressure and
        temperature by changing the structural parameters.

        Run during set_state() and set_composition(), as long as both
        state and composition have already been set. This function
        should not generally be needed by the user.
        """
        n0 = copy.copy(self.unrelaxed._scalar_solution.molar_fractions)

        def G_func(dq):
            n = n0 + np.einsum("ij, j", self.dndq, dq)
            self.unrelaxed._scalar_solution.set_composition(n)
            return self.unrelaxed._scalar_solution.molar_gibbs

        sol = minimize(G_func, np.zeros(len(self.dndq[0])), method="Nelder-Mead")
        assert sol.success
        n = n0 + np.einsum("ij, j", self.dndq, sol.x)
        self.unrelaxed.set_composition(n)
        AnisotropicSolution.set_composition(self, n)

    def set_state_with_volume(
        self, volume, temperature, pressure_guesses=[0.0e9, 10.0e9]
    ):
        """
        This function acts similarly to set_state, but takes volume and
        temperature as input to find the pressure. In order to ensure
        self-consistency, this function does not use any pressure functions
        from the material classes, but instead finds the pressure using the
        brentq root-finding and Nelder-Mead minimization methods.

        Composition should have been set before this function is used.

        :param volume: The desired molar volume of the mineral [m^3].
        :type volume: float

        :param temperature: The desired temperature of the mineral [K].
        :type temperature: float

        :param pressure_guesses: A list of floats denoting the initial
            low and high guesses for bracketing of the pressure [Pa].
            These guesses should preferably bound the correct pressure,
            but do not need to do so. More importantly,
            they should not lie outside the valid region of
            the equation of state. Defaults to [0.e9, 10.e9].
        :type pressure_guesses: list
        """
        ss = self.unrelaxed._scalar_solution
        n0 = copy.copy(ss.molar_fractions)

        # Initial set_state without changing structural parameters
        # to estimate pressure
        ss.set_state_with_volume(volume, temperature, pressure_guesses)
        # Store in a mutable so that it can be updated in the loop
        pressure = [ss.pressure]

        def F_func(dq):
            n = n0 + np.einsum("ij, j", self.dndq, dq)
            ss.set_composition(n)
            ss.set_state_with_volume(
                volume, temperature, [pressure[0] - 1.0e9, pressure[0] + 1.0e9]
            )
            pressure[0] = ss.pressure
            return ss.molar_helmholtz

        sol = minimize(F_func, np.zeros(len(self.dndq[0])), method="Nelder-Mead")
        assert sol.success

        # Run one more time with root
        F_func(sol.x)
        self.unrelaxed.set_composition(ss.molar_fractions)
        AnisotropicSolution.set_composition(self, ss.molar_fractions)
        self.set_state(pressure[0], temperature)

    @material_property
    def _depsdq_fixed_VT(self):
        """
        Gradient in strain with respect to structural state
        at fixed volume and temperature under hydrostatic conditions.
        """
        return np.einsum("mnq, qu->mnu", self.depsdn_fixed_VT, self.dndq)

    @material_property
    def _depsdT_fixed_volume(self):
        """
        Gradient in strain with respect to temperature
        at fixed volume and structural state under hydrostatic conditions.
        """
        beta_T = self.unrelaxed.isothermal_compressibility_tensor
        beta_TR = self.unrelaxed.isothermal_compressibility_reuss
        alpha = self.unrelaxed.thermal_expansivity_tensor
        alpha_V = self.unrelaxed.alpha

        return (alpha / alpha_V - beta_T / beta_TR) * alpha_V

    @material_property
    def _dbarepsdeps(self):
        """
        Gradient in non-hydrostatic strain with
        respect to total strain at fixed temperature
        and structural state under hydrostatic conditions.
        """
        beta_T = self.unrelaxed.isothermal_compressibility_tensor
        beta_TR = self.unrelaxed.isothermal_compressibility_reuss
        return depsdeps - np.einsum("pq, mn->pqmn", beta_T / beta_TR, np.eye(3))

    @material_property
    def _dSdq(self):
        """
        Gradient in strain with respect to structural state
        at constant volume and temperature under hydrostatic
        conditions.

        This function converts the partial entropies
        from Gibbs natural variables (P, T)
        to Helmholtz natural variables (V, T).
        """
        dSdn = (
            self._scalar_solution.partial_entropies
            - self._scalar_solution.alpha
            * self._scalar_solution.isothermal_bulk_modulus_reuss
            * self._scalar_solution.partial_volumes
        )
        return np.einsum("i, ij->j", dSdn, self.dndq)

    @material_property
    def _dPdq(self):
        """
        Pressure gradient with respect to structural state
        calculated at constant volume and temperature
        under hydrostatic conditions.

        This function converts from Gibbs natural variables
        (partial volumes as a function of P and T)
        to Helmholtz natural variables
        (partial pressures as a function of V and T).
        """
        dPdn = (
            self._scalar_solution.isothermal_bulk_modulus_reuss
            / self._scalar_solution.V
            * self._scalar_solution.partial_volumes
        )
        return np.einsum("i, ij->j", dPdn, self.dndq)

    @material_property
    def _d2Fdqdq_fixed_VT(self):
        """
        Helmholtz structural hessian
        calculated at constant volume and temperature
        under hydrostatic conditions.

        This function converts from Gibbs natural variables
        (the Gibbs compositional hessian at constant P and T)
        to Helmholtz natural variables
        (the Helmholtz compositional hessian
        at constant V and T).
        """
        d2Fdndn = (
            self._scalar_solution.gibbs_hessian
            + np.einsum(
                "i, j->ij",
                self._scalar_solution.partial_volumes,
                self._scalar_solution.partial_volumes,
            )
            * self._scalar_solution.isothermal_bulk_modulus_reuss
            / self._scalar_solution.V
        )
        return np.einsum("ij, ik, jl->kl", d2Fdndn, self.dndq, self.dndq)

    @material_property
    def _d2Fdqdq_fixed_epsT(self):
        """
        Second structure parameter derivative of the Helmholtz energy
        at constant strain and temperature.
        """
        C_T = self.unrelaxed.full_isothermal_stiffness_tensor
        return self._d2Fdqdq_fixed_VT + self.V * np.einsum(
            "kli, klmn, mnj->ij", self._depsdq_fixed_VT, C_T, self._depsdq_fixed_VT
        )

    @material_property
    def _d2Fdqdz(self):
        """
        Second derivatives of the Helmholtz energy with
        respect to composition and z (strain or temperature).
        """
        C_T = self.unrelaxed.full_isothermal_stiffness_tensor
        d2FdqdT = -self._dSdq + self.V * np.einsum(
            "kli, klmn, mn->i", self._depsdq_fixed_VT, C_T, self._depsdT_fixed_volume
        )

        d2Fdqdeps = -self.V * (
            np.einsum("i, mn->imn", self._dPdq, np.eye(3))
            + np.einsum(
                "kli, klpq, pqmn->imn", self._depsdq_fixed_VT, C_T, self._dbarepsdeps
            )
        )
        d2Fdqdeps = np.array([contract_stresses(arr) for arr in d2Fdqdeps])

        return np.concatenate((d2Fdqdeps, d2FdqdT[:, np.newaxis]), axis=1)

    @material_property
    def _d2Fdqdq_fixed_epsT_pinv(self):
        """
        The inverse of the second derivative of the Helmholtz energy
        with respect to the structure parameters at constant
        strain and temperature. Often referred to as the
        susceptibility matrix.
        """
        return np.linalg.pinv(self._d2Fdqdq_fixed_epsT)

    @material_property
    def dqdz_relaxed(self):
        """
        The change of the structure parameters with respect to
        strain and temperature that minimizes the Helmholtz
        energy.
        """
        return -np.einsum("kl, lj->kj", self._d2Fdqdq_fixed_epsT_pinv, self._d2Fdqdz)

    @material_property
    def _d2Fdzdz_Q(self):
        """
        Block matrix of V*C_T, V*pi, -c_eps/T
        at fixed Q
        """
        CT = self.unrelaxed.isothermal_stiffness_tensor
        pi = contract_stresses(self.unrelaxed.thermal_stress_tensor)
        c_eps = self.unrelaxed.molar_isometric_heat_capacity
        V = self.molar_volume
        T = self.temperature
        return np.block([[V * CT, V * pi[:, np.newaxis]], [V * pi, -c_eps / T]])

    @material_property
    def _d2Fdzdz(self):
        """
        Block matrix of V*C_T, V*pi, -c_eps/T
        under Helmholtz-minimizing Q
        """
        return self._d2Fdzdz_Q + np.einsum(
            "ki, kj->ij", self._d2Fdqdz, self.dqdz_relaxed
        )

    # The next three functions provide the three relaxed
    # second derivative properties from which other properties
    # may be derived.
    @material_property
    @copy_documentation(AnisotropicMineral.isothermal_stiffness_tensor)
    def isothermal_stiffness_tensor(self):
        return self._d2Fdzdz[:6, :6] / self.V

    @material_property
    @copy_documentation(AnisotropicMineral.thermal_stress_tensor)
    def thermal_stress_tensor(self):
        return expand_stresses(self._d2Fdzdz[6, :6] / self.V)

    @material_property
    @copy_documentation(AnisotropicMineral.molar_isometric_heat_capacity)
    def molar_isometric_heat_capacity(self):
        return -self._d2Fdzdz[6, 6] * self.temperature

    # The following tensor properties are all
    # second derivatives that are calculated in AnisotropicMineral
    # (rather than being derived from other properties),
    # and must therefore be redefined in relaxed materials

    # The full and Voigt isothermal stiffness tensors,
    # full and contracted isentropic stiffness and compliance tensors
    # Grueneisen tensor are derived properties
    # in AnisotropicMineral
    # and so do not need to be redefined.

    @material_property
    @copy_documentation(AnisotropicMineral.isothermal_compliance_tensor)
    def isothermal_compliance_tensor(self):
        return np.linalg.pinv(self.isothermal_stiffness_tensor)

    @material_property
    @copy_documentation(AnisotropicMineral.full_isothermal_compliance_tensor)
    def full_isothermal_compliance_tensor(self):
        S_Voigt = self.isothermal_compliance_tensor
        return voigt_notation_to_compliance_tensor(S_Voigt)

    @material_property
    @copy_documentation(AnisotropicMineral.thermal_expansivity_tensor)
    def thermal_expansivity_tensor(self):
        alpha = -np.einsum(
            "ijkl, kl",
            self.full_isothermal_compliance_tensor,
            self.thermal_stress_tensor,
        )
        return alpha

    # The following scalar properties are all
    # second derivatives that are calculated in AnisotropicMineral
    # or Solution (rather than being derived from other properties),
    # and must therefore be redefined in relaxed materials

    # The volumetric molar heat capacity and grueneisen parameter are
    # derived properties in AnisotropicMineral
    # and so do not need to be redefined.

    @material_property
    @copy_documentation(AnisotropicMineral.isothermal_compressibility_reuss)
    def isothermal_compressibility_reuss(self):
        return np.trace(self.isothermal_compressibility_tensor)

    @material_property
    @copy_documentation(AnisotropicMineral.isothermal_bulk_modulus_reuss)
    def isothermal_bulk_modulus_reuss(self):
        return 1.0 / self.isothermal_compressibility_reuss

    @material_property
    @copy_documentation(AnisotropicMineral.isentropic_compressibility_reuss)
    def isentropic_compressibility_reuss(self):
        return np.trace(self.isentropic_compressibility_tensor)

    @material_property
    @copy_documentation(AnisotropicMineral.isentropic_bulk_modulus_reuss)
    def isentropic_bulk_modulus_reuss(self):
        return 1.0 / self.isentropic_compressibility_reuss

    @material_property
    @copy_documentation(AnisotropicMineral.thermal_expansivity)
    def thermal_expansivity(self):
        return np.trace(self.thermal_expansivity_tensor)

    @material_property
    @copy_documentation(AnisotropicMineral.molar_heat_capacity_p)
    def molar_heat_capacity_p(self):
        alpha = self.thermal_expansivity_tensor
        pi = self.thermal_stress_tensor
        C_p = (
            self.molar_isometric_heat_capacity
            - self.V * self.temperature * np.einsum("ij, ij", alpha, pi)
        )
        return C_p
