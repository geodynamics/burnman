# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
from sympy import Matrix, nsimplify
from collections import OrderedDict

from .material import material_property, cached_property
from .mineral import Mineral
from .solutionmodel import MechanicalSolution
from .solutionmodel import PolynomialSolution
from .averaging_schemes import reuss_average_function

from ..utils.reductions import independent_row_indices
from ..utils.chemistry import sum_formulae, sort_element_list_to_IUPAC_order
from ..utils.misc import copy_documentation
from ..utils.math import complete_basis

import copy
from scipy.optimize import minimize


class Solution(Mineral):
    """
    This is the base class for all solutions.
    Site occupancies, endmember activities and the constant
    and pressure and temperature dependencies of the excess
    properties can be queried after using set_composition().
    States of the solution can only be queried after setting
    the pressure, temperature and composition using set_state().

    This class is available as :class:`burnman.Solution`.
    It uses an instance of :class:`burnman.SolutionModel` to
    calculate interaction terms between endmembers.

    All the solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T
    and P derivatives in J/K/mol and m^3/mol.

    The parameters are relevant to all solution models. Please
    see the documentation for individual models for details about
    other parameters.

    :param name: Name of the solution.
    :type name: string
    :param solution_model: The SolutionModel object defining the properties
        of the solution.
    :type solution_model: :class:`burnman.SolutionModel`
    :param molar_fractions: The molar fractions of each endmember
        in the solution. Can be reset using the set_composition() method.
    :type molar_fractions: numpy.array
    """

    def __init__(self, name=None, solution_model=None, molar_fractions=None):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.
        """
        Mineral.__init__(self)

        # Solution needs a method attribute to call Mineral.set_state().
        # Note that set_method() below will not change self.method
        self.method = "SolutionMethod"

        if name is not None:
            self.name = name
        if solution_model is not None:
            self.solution_model = solution_model

        # Equation of state
        for mbr in self.solution_model.endmembers:
            mbr[0].set_method(mbr[0].params["equation_of_state"])

        # Molar fractions
        if molar_fractions is not None:
            self.set_composition(molar_fractions)

    @cached_property
    def endmembers(self):
        return self.solution_model.endmembers

    def set_composition(self, molar_fractions):
        """
        Set the composition for this solution.
        Resets cached properties.

        :param molar_fractions: Molar abundance for each endmember,
            needs to sum to one.
        :type molar_fractions: list of float
        """
        assert len(self.solution_model.endmembers) == len(molar_fractions)

        if type(self.solution_model) is not MechanicalSolution:
            if np.abs(sum(molar_fractions) - 1.0) > 1.0e-4:
                raise ValueError(
                    "Molar fractions do not sum to one for "
                    "the current instance of "
                    f"<{self.name}>: {molar_fractions}"
                )

        if type(self.solution_model) is PolynomialSolution:
            self.solution_model.set_composition(molar_fractions)

        self.reset()
        self.molar_fractions = np.array(molar_fractions)

    def set_method(self, method):
        for i in range(self.n_endmembers):
            self.solution_model.endmembers[i][0].set_method(method)
        # note: do not set self.method here!
        self.reset()

    def set_state(self, pressure, temperature):
        if type(self.solution_model) is PolynomialSolution:
            self.solution_model.set_state(pressure, temperature)

        Mineral.set_state(self, pressure, temperature)
        for mbr in self.solution_model.endmembers:
            mbr[0].set_state(pressure, temperature)

    @material_property
    def formula(self):
        """
        Returns molar chemical formula of the solution.
        :rtype: Counter
        """
        return sum_formulae(self.endmember_formulae, self.molar_fractions)

    @material_property
    def site_occupancies(self):
        """
        :returns: The fractional occupancies of species on each site.
        :rtype: list of OrderedDicts
        """
        f = self.molar_fractions
        occs = np.einsum("ij, i", self.solution_model.endmember_occupancies, f)
        site_occs = []
        k = 0
        for site in self.solution_model.sites:
            site_occs.append(OrderedDict())
            for species in site:
                site_occs[-1][species] = occs[k]
                k += 1

        return site_occs

    def site_formula(self, precision=2):
        """
        Returns the molar chemical formula of the solution with
            site occupancies. For example, [Mg0.4Fe0.6]2SiO4.

        :param precision: Precision with which to print the site occupancies
        :type precision: int

        :returns: Molar chemical formula of the solution with site occupancies
        :rtype: str
        """
        split_empty = self.solution_model.empty_formula.split("[")
        formula = split_empty[0]
        for i, site_occs in enumerate(self.site_occupancies):
            formula += "["
            for species, occ in site_occs.items():
                if np.abs(occ) < 1.0e-12:
                    occ = np.abs(occ)
                formula += f"{species}{occ:0.{precision}f}"
            formula += split_empty[i + 1]
        return formula

    @material_property
    def activities(self):
        """
        Returns a list of endmember activities [unitless].
        """
        return self.solution_model.activities(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def activity_coefficients(self):
        """
        Returns a list of endmember activity coefficients
        (gamma = activity / ideal activity) [unitless].
        """
        return self.solution_model.activity_coefficients(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_internal_energy(self):
        """
        Returns molar internal energy of the mineral [J/mol].
        Aliased with self.energy
        """
        return self.molar_helmholtz + self.temperature * self.molar_entropy

    @material_property
    def excess_partial_gibbs(self):
        """
        Returns excess partial molar gibbs free energy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_gibbs_free_energies(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def excess_partial_volumes(self):
        """
        Returns excess partial volumes [m^3].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_volumes(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def excess_partial_entropies(self):
        """
        Returns excess partial entropies [J/K].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_entropies(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def partial_gibbs(self):
        """
        Returns endmember partial molar gibbs free energy [J/mol].
        Property specific to solutions.
        """
        return (
            np.array([mbr[0].gibbs for mbr in self.solution_model.endmembers])
            + self.excess_partial_gibbs
        )

    @material_property
    def partial_volumes(self):
        """
        Returns endmember partial volumes [m^3].
        Property specific to solutions.
        """
        return (
            np.array([mbr[0].molar_volume for mbr in self.solution_model.endmembers])
            + self.excess_partial_volumes
        )

    @material_property
    def partial_entropies(self):
        """
        Returns endmember partial entropies [J/K].
        Property specific to solutions.
        """
        return (
            np.array([mbr[0].molar_entropy for mbr in self.solution_model.endmembers])
            + self.excess_partial_entropies
        )

    @material_property
    def excess_gibbs(self):
        """
        Returns molar excess gibbs free energy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_gibbs_free_energy(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def gibbs_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the Gibbs free energy [J]. Property specific to solutions.
        """
        return self.solution_model.gibbs_hessian(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def entropy_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the entropy [J/K]. Property specific to solutions.
        """
        return self.solution_model.entropy_hessian(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def volume_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the volume [m^3]. Property specific to solutions.
        """
        return self.solution_model.volume_hessian(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_gibbs(self):
        """
        Returns molar Gibbs free energy of the solution [J/mol].
        Aliased with self.gibbs.
        """
        f = self.molar_fractions
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].gibbs * f[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_gibbs
        )

    @material_property
    def molar_helmholtz(self):
        """
        Returns molar Helmholtz free energy of the solution [J/mol].
        Aliased with self.helmholtz.
        """
        return self.molar_gibbs - self.pressure * self.molar_volume

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the solution [kg/mol].
        """
        return sum(
            [
                self.solution_model.endmembers[i][0].molar_mass
                * self.molar_fractions[i]
                for i in range(self.n_endmembers)
            ]
        )

    @material_property
    def excess_volume(self):
        """
        Returns excess molar volume of the solution [m^3/mol].
        Specific property for solutions.
        """
        return self.solution_model.excess_volume(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the solution [m^3/mol].
        Aliased with self.V.
        """
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].molar_volume
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_volume
        )

    @material_property
    def density(self):
        """
        Returns density of the solution [kg/m^3].
        Aliased with self.rho.
        """
        return self.molar_mass / self.molar_volume

    @material_property
    def excess_entropy(self):
        """
        Returns excess molar entropy [J/K/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_entropy(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropy of the solution [J/K/mol].
        Aliased with self.S.
        """
        f = self.molar_fractions
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].S * f[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_entropy
        )

    @material_property
    def excess_enthalpy(self):
        """
        Returns excess molar enthalpy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_enthalpy(
            self.pressure, self.temperature, self.molar_fractions
        )

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the solution [J/mol].
        Aliased with self.H.
        """
        f = self.molar_fractions
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].H * f[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.excess_enthalpy
        )

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        Returns isothermal bulk modulus of the solution [Pa].
        Aliased with self.K_T.
        """
        return (
            self.V
            * 1.0
            / (
                sum(
                    [
                        self.solution_model.endmembers[i][0].V
                        / (
                            self.solution_model.endmembers[i][
                                0
                            ].isothermal_bulk_modulus_reuss
                        )
                        * self.molar_fractions[i]
                        for i in range(self.n_endmembers)
                    ]
                )
                + self.solution_model.VoverKT_excess()
            )
        )

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Returns adiabatic bulk modulus of the solution [Pa].
        Aliased with self.K_S.
        """
        if self.temperature < 1e-10:
            return self.isothermal_bulk_modulus_reuss
        else:
            return (
                self.isothermal_bulk_modulus_reuss
                * self.molar_heat_capacity_p
                / self.molar_heat_capacity_v
            )

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns isothermal compressibility of the solution.
        (or inverse isothermal bulk modulus) [1/Pa].
        Aliased with self.K_T.
        """
        return 1.0 / self.isothermal_bulk_modulus_reuss

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        Returns adiabatic compressibility of the solution.
        (or inverse adiabatic bulk modulus) [1/Pa].
        Aliased with self.K_S.
        """
        return 1.0 / self.isentropic_bulk_modulus_reuss

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the solution [Pa].
        Aliased with self.G.
        """
        G_list = np.fromiter(
            (e[0].G for e in self.solution_model.endmembers),
            dtype=float,
            count=self.n_endmembers,
        )
        return reuss_average_function(self.molar_fractions, G_list)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the solution [m/s].
        Aliased with self.v_p.
        """
        return np.sqrt(
            (self.isentropic_bulk_modulus_reuss + 4.0 / 3.0 * self.shear_modulus)
            / self.density
        )

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the solution [m/s].
        Aliased with self.v_phi.
        """
        return np.sqrt(self.isentropic_bulk_modulus_reuss / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the solution [m/s].
        Aliased with self.v_s.
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the solution [unitless].
        Aliased with self.gr.
        """
        if self.temperature < 1e-10:
            return float("nan")
        else:
            return (
                self.thermal_expansivity
                * self.isothermal_bulk_modulus_reuss
                * self.molar_volume
                / self.molar_heat_capacity_v
            )

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha)
        of the solution [1/K].
        Aliased with self.alpha.
        """
        return (1.0 / self.V) * (
            sum(
                [
                    self.solution_model.endmembers[i][0].alpha
                    * self.solution_model.endmembers[i][0].V
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.solution_model.alphaV_excess()
        )

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the
        solution [J/K/mol].
        Aliased with self.C_v.
        """
        return (
            self.molar_heat_capacity_p
            - self.molar_volume
            * self.temperature
            * self.thermal_expansivity
            * self.thermal_expansivity
            * self.isothermal_bulk_modulus_reuss
        )

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure
        of the solution [J/K/mol].
        Aliased with self.C_p.
        """
        return (
            sum(
                [
                    self.solution_model.endmembers[i][0].molar_heat_capacity_p
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)
                ]
            )
            + self.solution_model.Cp_excess()
        )

    @cached_property
    def stoichiometric_matrix(self):
        """
        A sympy Matrix where each element M[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """

        def f(i, j):
            e = self.elements[j]
            if e in self.endmember_formulae[i]:
                return nsimplify(self.endmember_formulae[i][e])
            else:
                return 0

        return Matrix(len(self.endmember_formulae), len(self.elements), f)

    @cached_property
    def stoichiometric_array(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of atoms of element[j] in endmember[i].
        """
        return np.array(self.stoichiometric_matrix)

    @cached_property
    def reaction_basis(self):
        """
        An array where each element arr[i,j] corresponds
        to the number of moles of endmember[j] involved in reaction[i].
        """
        reaction_basis = np.array(
            [v[:] for v in self.stoichiometric_matrix.T.nullspace()], dtype=float
        )

        if len(reaction_basis) == 0:
            reaction_basis = np.empty((0, len(self.endmember_names)))

        return reaction_basis

    @cached_property
    def n_reactions(self):
        """
        The number of reactions in reaction_basis.
        """
        return len(self.reaction_basis[:, 0])

    @cached_property
    def compositional_basis(self):
        """_summary_

        :return: _description_
        :rtype: _type_
        """
        return complete_basis(self.reaction_basis)[self.n_reactions :]

    @cached_property
    def independent_element_indices(self):
        """
        A list of an independent set of element indices. If the amounts of
        these elements are known (element_amounts),
        the amounts of the other elements can be inferred by
        -compositional_null_basis[independent_element_indices].dot(element_amounts).
        """
        return sorted(independent_row_indices(self.stoichiometric_matrix.T))

    @cached_property
    def dependent_element_indices(self):
        """
        The element indices not included in the independent list.
        """
        return [
            i
            for i in range(len(self.elements))
            if i not in self.independent_element_indices
        ]

    @cached_property
    def compositional_null_basis(self):
        """
        An array N such that N.b = 0 for all bulk compositions that can
        be produced with a linear sum of the endmembers in the solution.
        """
        null = self.stoichiometric_matrix.nullspace()
        null_basis = np.array([v[:] for v in null])

        M = null_basis[:, self.dependent_element_indices]
        assert M.shape[0] == M.shape[1]
        assert (np.abs(M - np.eye(M.shape[0])) < 1.0e-20).all()

        return null_basis

    @cached_property
    def endmember_formulae(self):
        """
        A list of formulae for all the endmember in the solution.
        """
        mbrs = self.solution_model.endmembers
        return [mbr[0].params["formula"] for mbr in mbrs]

    @cached_property
    def endmember_names(self):
        """
        A list of names for all the endmember in the solution.
        """
        return [mbr[0].name for mbr in self.solution_model.endmembers]

    @cached_property
    def n_endmembers(self):
        """
        The number of endmembers in the solution.
        """
        return len(self.solution_model.endmembers)

    @cached_property
    def elements(self):
        """
        A list of the elements which could be contained in the solution,
        returned in the IUPAC element order.
        """
        keys = []
        for f in self.endmember_formulae:
            keys.extend(f.keys())

        return sort_element_list_to_IUPAC_order(set(keys))


SolidSolution = Solution


class RelaxedSolution(Solution):
    """
    A class implementing a relaxed solution
    model.
    This class is derived from Solution,
    and inherits most of the methods from that class.

    Instantiation of a RelaxedSolution involves passing
    an Solution, plus a set of vectors that represent rapid
    deformation modes. For example, a solution of MgO, FeHSO and FeLSO
    (high and low spin wuestite) can rapidly change proportion of
    high spin and low spin iron, and so a single vector should be passed:
    np.array([[0., -1., 1.]]) or some multiple thereof.

    States of the mineral can only be queried after setting the
    pressure and temperature using set_state() and the composition using
    set_composition().

    This class is available as ``burnman.RelaxedSolution``.
    """

    def __init__(
        self,
        solution,
        relaxation_vectors,
        unrelaxed_vectors,
    ):
        # Make an attribute with the unrelaxed solution
        self.unrelaxed = solution

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
            molar_fractions = solution.molar_fractions
        except AttributeError:
            molar_fractions = None

        # Give the relaxed solution the same base properties as the
        # unrelaxed solution
        Solution.__init__(
            self,
            name=solution.name,
            solution_model=solution.solution_model,
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
            Solution.set_state(self, pressure, temperature)

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
            Solution.set_composition(self, n)
            self.unrelaxed.set_composition(self.molar_fractions)

    def _relax_at_PTX(self):
        """
        Minimizes the Gibbs energy at constant pressure and
        temperature by changing the structural parameters.

        Run during set_state() and set_composition(), as long as both
        state and composition have already been set. This function
        should not generally be needed by the user.
        """
        n0 = copy.copy(self.unrelaxed.molar_fractions)

        def G_func(dq):
            n = n0 + np.einsum("ij, j", self.dndq, dq)
            self.unrelaxed.set_composition(n)
            return self.unrelaxed.molar_gibbs

        sol = minimize(G_func, np.zeros(len(self.dndq[0])), method="Nelder-Mead")
        assert sol.success
        n = n0 + np.einsum("ij, j", self.dndq, sol.x)
        self.unrelaxed.set_composition(n)
        Solution.set_composition(self, n)

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
        ss = self.unrelaxed
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
        Solution.set_composition(self, ss.molar_fractions)
        self.set_state(pressure[0], temperature)

    @material_property
    def _d2Gdqdq_fixed_PT(self):
        """
        Gibbs structural hessian
        calculated at constant pressure and temperature.
        """
        return np.einsum(
            "ij, ik, jl->kl", self.unrelaxed.gibbs_hessian, self.dndq, self.dndq
        )

    @material_property
    def _d2Gdqdz(self):
        """
        Second derivatives of the Helmholtz energy with
        respect to composition and z (pressure or temperature).
        """
        dVdq = np.einsum("i, ij->j", self.unrelaxed.partial_volumes, self.dndq)
        dSdq = np.einsum("i, ij->j", self.unrelaxed.partial_entropies, self.dndq)

        return np.array([dVdq, -dSdq]).T

    @material_property
    def _d2Gdqdq_fixed_PT_pinv(self):
        """
        The second derivative of the Helmholtz energy
        with respect to the structure parameters at constant
        strain and temperature. Often referred to as the
        susceptibility matrix.
        """
        return np.linalg.pinv(self._d2Gdqdq_fixed_PT)

    @material_property
    def dqdz_relaxed(self):
        """
        The change of the structure parameters with respect to
        strain and temperature that minimizes the Helmholtz
        energy.
        """
        return -np.einsum("kl, lj->kj", self._d2Gdqdq_fixed_PT_pinv, self._d2Gdqdz)

    @material_property
    def _d2Gdzdz_Q(self):
        """
        Block matrix of -V*beta_TR, V*alpha, -c_p/T
        at fixed Q
        """
        beta_TR = self.unrelaxed.isothermal_compressibility_reuss
        alpha = self.unrelaxed.thermal_expansivity
        c_p = self.unrelaxed.molar_heat_capacity_p
        V = self.molar_volume
        T = self.temperature
        return np.array([[-V * beta_TR, V * alpha], [V * alpha, -c_p / T]])

    @material_property
    def _d2Gdzdz(self):
        """
        Block matrix of -V*beta_TR, V*alpha, -c_p/T
        under Gibbs-minimizing Q
        """
        return self._d2Gdzdz_Q + np.einsum(
            "ki, kj->ij", self._d2Gdqdz, self.dqdz_relaxed
        )

    # The following scalar properties are all
    # second derivatives that are calculated in AnisotropicMineral
    # or Solution (rather than being derived from other properties),
    # and must therefore be redefined in relaxed materials

    # The volumetric molar heat capacity and grueneisen parameter are
    # derived properties in AnisotropicMineral
    # and so do not need to be redefined.

    @material_property
    @copy_documentation(Solution.isothermal_compressibility_reuss)
    def isothermal_compressibility_reuss(self):
        return -self._d2Gdzdz[0, 0] / self.V

    @material_property
    @copy_documentation(Solution.thermal_expansivity)
    def thermal_expansivity(self):
        return self._d2Gdzdz[0, 1] / self.V

    @material_property
    @copy_documentation(Solution.molar_heat_capacity_p)
    def molar_heat_capacity_p(self):
        return -self._d2Gdzdz[1, 1] * self.T

    @material_property
    @copy_documentation(Solution.isothermal_bulk_modulus_reuss)
    def isothermal_bulk_modulus_reuss(self):
        return 1.0 / self.isothermal_compressibility_reuss

    @material_property
    @copy_documentation(Solution.molar_heat_capacity_v)
    def molar_heat_capacity_v(self):
        return (
            self.molar_heat_capacity_p
            - self.molar_volume
            * self.temperature
            * self.thermal_expansivity
            * self.thermal_expansivity
            * self.isothermal_bulk_modulus_reuss
        )

    @material_property
    @copy_documentation(Solution.isentropic_bulk_modulus_reuss)
    def isentropic_bulk_modulus_reuss(self):
        if self.temperature < 1.0e-10:
            return self.isothermal_bulk_modulus_reuss
        else:
            return (
                self.isothermal_bulk_modulus_reuss
                * self.molar_heat_capacity_p
                / self.molar_heat_capacity_v
            )

    @material_property
    @copy_documentation(Solution.isentropic_compressibility_reuss)
    def isentropic_compressibility_reuss(self):
        return 1.0 / self.isentropic_bulk_modulus_reuss
