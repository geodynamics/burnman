# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.
import numpy as np
import copy
from scipy.linalg import expm, logm
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
    voigt_notation_to_stiffness_tensor,
    contract_compliances,
    contract_stiffnesses,
)


class AnisotropicSolution(Solution, AnisotropicMineral):
    """
    A class implementing the anisotropic solution model described
    in :cite:`Myhill2024`.
    This class is derived from Solution and AnisotropicMineral,
    and inherits most of the methods from those classes.

    Instantiation of an AnisotropicSolution is similar to that of a scalar
    Solution, except that each of the endmembers must be an instance of
    an AnisotropicMineral, and an additional function and parameter dictionary
    are passed to the constructor of AnisotropicSolution that describe
    excess contributions to the anisotropic state tensor (Psi_xs)
    and its derivatives with respect to volume and temperature.
    The function arguments should be ln(V), Pth,
    X (a vector) and the parameter dictionary, in that order.
    The output variables Psi_xs_Voigt, dPsidf_Pth_Voigt_xs and
    dPsidPth_f_Voigt_xs (all 6x6 matrices) must be returned in that
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
        Solution.__init__(self, name, solution_model)

        # Store a scalar copy of the solution model to speed up set_state
        scalar_solution_model = copy.copy(solution_model)
        scalar_solution_model.endmembers = [
            [mbr.isotropic_mineral, formula] for mbr, formula in self.endmembers
        ]
        self._scalar_solution = Solution(name, scalar_solution_model, molar_fractions)

        self._logm_M_T_0_mbr = np.array(
            [logm(m[0].cell_vectors_0) for m in self.endmembers]
        )
        self.anisotropic_params = anisotropic_parameters
        self.psi_excess_function = psi_excess_function

        if molar_fractions is not None:
            self.set_composition(molar_fractions)

    def set_state(self, pressure, temperature):
        # Set solution conditions
        ss = self._scalar_solution
        if not hasattr(ss, "molar_fractions"):
            raise Exception("To use this EoS, you must first set the composition")

        ss.set_state(pressure, temperature)
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
        PsiI_mbr = np.array([m[0]._PsiI for m in mbrs])
        dPsidf_T_Voigt_mbr = np.array([m[0]._dPsidf_T_Voigt for m in mbrs])
        dPsiIdf_T_mbr = np.array([m[0]._dPsiIdf_T for m in mbrs])
        dPsiIdT_f_mbr = np.array([m[0]._dPsiIdT_f for m in mbrs])

        fr = self.molar_fractions
        PsiI_ideal = np.einsum("p,pij->ij", fr, PsiI_mbr)
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
        Psi_xs_Voigt, dPsidf_Pth_Voigt_xs, dPsidPth_f_Voigt_xs = out
        Psi_xs_full = voigt_notation_to_compliance_tensor(Psi_xs_Voigt)
        PsiI_xs = np.einsum("ijkl, kl", Psi_xs_full, np.eye(3))

        # Change of variables: (f, Pth) -> Psi(f, T)
        aK_T = ss.alpha * ss.isothermal_bulk_modulus_reuss
        out = convert_f_Pth_to_f_T_derivatives(
            dPsidf_Pth_Voigt_xs, dPsidPth_f_Voigt_xs, aK_T, self.dPthdf
        )
        dPsidf_T_Voigt_xs, dPsiIdf_T_xs, dPsiIdT_f_xs = out

        out = deformation_gradient_alpha_and_compliance(
            self.alpha,
            self.isothermal_compressibility_reuss,
            PsiI_ideal + PsiI_xs,
            dPsidf_T_Voigt_ideal + dPsidf_T_Voigt_xs,
            dPsiIdf_T_ideal + dPsiIdf_T_xs,
            dPsiIdT_f_ideal + dPsiIdT_f_xs,
        )
        self._unrotated_F, self._unrotated_alpha, self._unrotated_S_T_Voigt = out

    def set_composition(self, molar_fractions):
        self._scalar_solution.set_composition(molar_fractions)
        Solution.set_composition(self, molar_fractions)
        self.cell_vectors_0 = expm(
            np.einsum("p, pij->ij", molar_fractions, self._logm_M_T_0_mbr)
        )
        # Calculate all other required properties
        if self.pressure is not None:
            self.set_state(self.pressure, self.temperature)
