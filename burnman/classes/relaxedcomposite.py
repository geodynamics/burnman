# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
from scipy.linalg import block_diag
from collections import Counter

from .material import material_property
from .composite import Composite
from ..utils.misc import copy_documentation
from ..tools.equilibration import equilibrate


class RelaxedComposite(Composite):
    """
    A class implementing a relaxed Composite mineral assemblage.
    Unlike a standard Composite, the set_state() method of
    this class uses the `burnman.equilibrate` function to find
    the equilibrium phase fractions and compositions at given
    pressure and temperature. Once equilibrated, the
    thermodynamic properties of the relaxed composite can be
    queried in the same way as for a standard Composite.
    One additional feature of this class is that properties
    related to the second derivatives of the Gibbs energy
    (e.g. compressibility, thermal expansivity, heat capacity)
    can be calculated taking into account relaxation of one or more
    user-defined isochemical reaction vectors via the relaxation_vectors
    argument to the constructor. If complete relaxation is desired,
    the user can pass relaxation_vectors=composite.reaction_basis.
    More limited relaxation can be specified
    by passing a linear subset of the reaction basis vectors.
    For example, on short timescales,a composite of olivine and wadsleyite
    might be able to exchange Mg and Fe, even if no grain growth or
    reduction can occur. In this case, a single vector representing the
    exchange of Mg and Fe between the two minerals can be used to relax
    the composition at given P and T:
    relaxation_vectors=[[1, -1, -1, 1]], where the first two
    entries correspond to olivine (fo, fa) and the second two to
    wadsleyite (mwd, fwd). This is a subset of the full reaction basis,
    which would be: [[1, 0, -1, 0], [0, 1, 0, -1]].

    This class is available as ``burnman.RelaxedComposite``.
    """

    def __init__(self, composite, relaxation_vectors):
        # Make an attribute with the unrelaxed composite
        self.unrelaxed = composite

        # The relaxation vectors
        self.n_relaxation_vectors = len(relaxation_vectors)

        # Check the shape of the relaxation vectors:
        # should be a 2D array with n_endmembers rows and
        # n_relaxation_vectors columns
        self.dndq = np.array(relaxation_vectors).T
        assert len(self.dndq.shape) == 2
        assert len(self.dndq) == self.unrelaxed.n_endmembers

        # Check that the relaxation vectors are isochemical
        delta_c = self.unrelaxed.stoichiometric_array.T @ self.dndq
        if not np.allclose(delta_c, 0.0):
            raise ValueError("Relaxation vectors must be isochemical.")

        try:
            molar_fractions = composite.molar_fractions
        except AttributeError:
            molar_fractions = None

        # Give the relaxed composite the same base properties as the
        # unrelaxed composite
        Composite.__init__(
            self,
            phases=composite.phases,
            fractions=molar_fractions,
            fraction_type="molar",
            name=composite.name,
        )

        self.number_of_moles = composite.number_of_moles

    def set_state(self, pressure, temperature, bulk_composition=None):
        """
        Sets the state of the composite. Also relaxes the
        structure parameters if set_composition() has already
        been used and if the relaxed argument has been set to
        True.

        :param pressure: The pressure of the solution [Pa]
        :type pressure: float
        :param temperature: The temperature of the solution [K]
        :type temperature: float
        :param bulk_composition: The bulk composition of the
            composite in terms of elemental amounts. If None,
            the bulk composition stored in the unrelaxed
            composite is used.
        :type bulk_composition: dict
        """
        if bulk_composition is None:
            bulk_composition = self.unrelaxed.formula

        self.unrelaxed.set_state(pressure, temperature)

        sol, _ = equilibrate(
            bulk_composition, self.unrelaxed, [["P", pressure], ["T", temperature]]
        )
        assert sol.success, "Equilibration failed in RelaxedComposite.set_state()"

        # Set the relaxed composite to the equilibrated state
        Composite.set_state(self, pressure, temperature)
        Composite.set_fractions(self, self.unrelaxed.molar_fractions)
        self.number_of_moles = self.unrelaxed.number_of_moles

    @material_property
    def _d2Gdqdq_fixed_PT(self):
        """
        Gibbs structural hessian
        calculated at constant pressure and temperature.
        """
        phases = self.unrelaxed.phases
        molar_amounts = self.unrelaxed.molar_fractions * self.unrelaxed.number_of_moles
        hessians = [
            (
                phase.gibbs_hessian / molar_amounts[i]
                if hasattr(phase, "gibbs_hessian")
                else np.array([[0.0]])
            )
            for i, phase in enumerate(phases)
        ]
        H = block_diag(*hessians)
        return self.dndq.T @ H @ self.dndq

    @material_property
    def _d2Gdqdz(self):
        """
        Second derivatives of the Gibbs energy with
        respect to composition and z (pressure and temperature).
        """
        dVdq = self.unrelaxed.endmember_partial_volumes @ self.dndq
        dSdq = self.unrelaxed.endmember_partial_entropies @ self.dndq

        return np.array([dVdq, -dSdq]).T

    @material_property
    def _d2Gdqdq_fixed_PT_pinv(self):
        """
        The second derivative of the Gibbs energy
        with respect to the structure parameters at constant
        pressure and temperature. Often referred to as the
        susceptibility matrix.
        """
        return np.linalg.pinv(self._d2Gdqdq_fixed_PT, rcond=1.0e-22)

    @material_property
    def dqdz_relaxed(self):
        """
        The change of the structure parameters with respect to
        pressure and temperature that minimizes the Gibbs
        energy.
        """
        return -self._d2Gdqdq_fixed_PT_pinv @ self._d2Gdqdz

    @material_property
    def _d2Gdzdz_Q(self):
        """
        Block matrix of -V*beta_TR, V*alpha, -c_p/T
        at fixed Q
        """
        beta_TR = self.unrelaxed.isothermal_compressibility_reuss
        alpha = self.unrelaxed.thermal_expansivity
        c_p = self.unrelaxed.heat_capacity_p
        V = self.unrelaxed.volume
        T = self.unrelaxed.temperature
        return np.array([[-V * beta_TR, V * alpha], [V * alpha, -c_p / T]])

    @material_property
    def _d2Gdzdz(self):
        """
        Block matrix of -V*beta_TR, V*alpha, -c_p/T
        under Gibbs-minimizing Q
        """
        return self._d2Gdzdz_Q + self._d2Gdqdz.T @ self.dqdz_relaxed

    # The following scalar properties are second derivatives that are
    # calculated in Composite (rather than being derived from other
    # properties), and must therefore be redefined in relaxed materials.

    # The volumetric molar heat capacity and Grueneisen parameter
    # are derived properties in Composite and so do not need to be redefined.

    @material_property
    @copy_documentation(Composite.isothermal_compressibility_reuss)
    def isothermal_compressibility_reuss(self):
        return -self._d2Gdzdz[0, 0] / self.volume

    @material_property
    @copy_documentation(Composite.thermal_expansivity)
    def thermal_expansivity(self):
        return self._d2Gdzdz[0, 1] / self.volume

    @material_property
    @copy_documentation(Composite.molar_heat_capacity_p)
    def molar_heat_capacity_p(self):
        return -self._d2Gdzdz[1, 1] * self.T / self.number_of_moles

    @material_property
    @copy_documentation(Composite.isothermal_bulk_modulus_reuss)
    def isothermal_bulk_modulus_reuss(self):
        return 1.0 / self.isothermal_compressibility_reuss
