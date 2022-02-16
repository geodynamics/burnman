# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
from sympy import Matrix, nsimplify
from .material import material_property, cached_property
from .mineral import Mineral
from .solutionmodel import SolutionModel
from .solutionmodel import MechanicalSolution, IdealSolution
from .solutionmodel import SymmetricRegularSolution, AsymmetricRegularSolution
from .solutionmodel import SubregularSolution
from .averaging_schemes import reuss_average_function

from ..utils.reductions import independent_row_indices
from ..utils.chemistry import sum_formulae, sort_element_list_to_IUPAC_order


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

    Parameters
    ----------
    name : string
        Name of the solution
    solution_type : string
        String determining which SolutionModel to use. One of 'mechanical',
        'ideal', 'symmetric', 'asymmetric' or 'subregular'.
    endmembers : list of lists
        List of endmembers in this solution. The first item of each
        list should be a :class:`burnman.Mineral` object. The second item
        should be a string with the site formula of the endmember.
    molar_fractions : numpy array (optional)
        The molar fractions of each endmember in the solution.
        Can be reset using the set_composition() method.
    """

    def __init__(self,
                 name=None,
                 solution_type=None,
                 endmembers=None,
                 energy_interaction=None,
                 volume_interaction=None,
                 entropy_interaction=None,
                 energy_ternary_terms=None,
                 volume_ternary_terms=None,
                 entropy_ternary_terms=None,
                 alphas=None,
                 molar_fractions=None):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.
        """
        Mineral.__init__(self)

        # Solution needs a method attribute to call Mineral.set_state().
        # Note that set_method() below will not change self.method
        self.method = 'SolutionMethod'

        if name is not None:
            self.name = name
        if solution_type is not None:
            self.solution_type = solution_type
        if endmembers is not None:
            self.endmembers = endmembers
        if energy_interaction is not None:
            self.energy_interaction = energy_interaction
        if volume_interaction is not None:
            self.volume_interaction = volume_interaction
        if entropy_interaction is not None:
            self.entropy_interaction = entropy_interaction
        if energy_ternary_terms is not None:
            self.energy_ternary_terms = energy_ternary_terms
        if volume_ternary_terms is not None:
            self.volume_ternary_terms = volume_ternary_terms
        if entropy_ternary_terms is not None:
            self.entropy_ternary_terms = entropy_ternary_terms
        if alphas is not None:
            self.alphas = alphas
        if endmembers is not None:
            self.endmembers = endmembers

        if hasattr(self, 'endmembers') is False:
            raise Exception("'endmembers' attribute missing "
                            "from solution")

        # Set default solution model type
        if hasattr(self, 'solution_type'):
            if self.solution_type == 'mechanical':
                self.solution_model = MechanicalSolution(self.endmembers)
            elif self.solution_type == 'ideal':
                self.solution_model = IdealSolution(self.endmembers)
            else:
                if hasattr(self, 'energy_interaction') is False:
                    self.energy_interaction = None
                if hasattr(self, 'volume_interaction') is False:
                    self.volume_interaction = None
                if hasattr(self, 'entropy_interaction') is False:
                    self.entropy_interaction = None

                if self.solution_type == 'symmetric':
                    self.solution_model = SymmetricRegularSolution(
                        self.endmembers, self.energy_interaction,
                        self.volume_interaction, self.entropy_interaction)
                elif self.solution_type == 'asymmetric':
                    if hasattr(self, 'alphas') is False:
                        raise Exception(
                            "'alphas' attribute missing from solution")
                    self.solution_model = AsymmetricRegularSolution(
                        self.endmembers, self.alphas, self.energy_interaction,
                        self.volume_interaction, self.entropy_interaction)
                elif self.solution_type == 'subregular':
                    if hasattr(self, 'energy_ternary_terms') is False:
                        self.energy_ternary_terms = None
                    if hasattr(self, 'volume_ternary_terms') is False:
                        self.volume_ternary_terms = None
                    if hasattr(self, 'entropy_ternary_terms') is False:
                        self.entropy_ternary_terms = None

                    self.solution_model = SubregularSolution(
                        self.endmembers,
                        self.energy_interaction,  self.volume_interaction,
                        self.entropy_interaction,
                        self.energy_ternary_terms, self.volume_ternary_terms,
                        self.entropy_ternary_terms)
                else:
                    raise Exception("Solution model type "
                                    + self.solution_type + "not recognised.")
        else:
            self.solution_model = SolutionModel()

        # Equation of state
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(
                self.endmembers[i][0].params['equation_of_state'])

        # Molar fractions
        if molar_fractions is not None:
            self.set_composition(molar_fractions)

    def get_endmembers(self):
        return self.endmembers

    def set_composition(self, molar_fractions):
        """
        Set the composition for this solution.
        Resets cached properties.

        Parameters
        ----------
        molar_fractions: list of float
            molar abundance for each endmember, needs to sum to one.
        """
        assert(len(self.endmembers) == len(molar_fractions))

        if self.solution_type != 'mechanical':
            assert(sum(molar_fractions) > 0.9999)
            assert(sum(molar_fractions) < 1.0001)

        self.reset()
        self.molar_fractions = np.array(molar_fractions)

    def set_method(self, method):
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_method(method)
        # note: do not set self.method here!
        self.reset()

    def set_state(self, pressure, temperature):

        Mineral.set_state(self, pressure, temperature)
        for i in range(self.n_endmembers):
            self.endmembers[i][0].set_state(pressure, temperature)

    @material_property
    def formula(self):
        """
        Returns molar chemical formula of the solution.
        """
        return sum_formulae(self.endmember_formulae, self.molar_fractions)

    @material_property
    def activities(self):
        """
        Returns a list of endmember activities [unitless].
        """
        return self.solution_model.activities(self.pressure, self.temperature,
                                              self.molar_fractions)

    @material_property
    def activity_coefficients(self):
        """
        Returns a list of endmember activity coefficients
        (gamma = activity / ideal activity) [unitless].
        """
        return self.solution_model.activity_coefficients(self.pressure,
                                                         self.temperature,
                                                         self.molar_fractions)

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
        return self.solution_model.excess_partial_gibbs_free_energies(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def excess_partial_volumes(self):
        """
        Returns excess partial volumes [m^3].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_volumes(self.pressure,
                                                          self.temperature,
                                                          self.molar_fractions)

    @material_property
    def excess_partial_entropies(self):
        """
        Returns excess partial entropies [J/K].
        Property specific to solutions.
        """
        return self.solution_model.excess_partial_entropies(self.pressure,
                                                            self.temperature,
                                                            self.molar_fractions)

    @material_property
    def partial_gibbs(self):
        """
        Returns endmember partial molar gibbs free energy [J/mol].
        Property specific to solutions.
        """
        return (np.array([self.endmembers[i][0].gibbs
                          for i in range(self.n_endmembers)])
                + self.excess_partial_gibbs)

    @material_property
    def partial_volumes(self):
        """
        Returns endmember partial volumes [m^3].
        Property specific to solutions.
        """
        return (np.array([self.endmembers[i][0].molar_volume
                          for i in range(self.n_endmembers)])
                + self.excess_partial_volumes)

    @material_property
    def partial_entropies(self):
        """
        Returns endmember partial entropies [J/K].
        Property specific to solutions.
        """
        return (np.array([self.endmembers[i][0].molar_entropy
                          for i in range(self.n_endmembers)])
                + self.excess_partial_entropies)

    @material_property
    def excess_gibbs(self):
        """
        Returns molar excess gibbs free energy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_gibbs_free_energy(self.pressure,
                                                            self.temperature,
                                                            self.molar_fractions)

    @material_property
    def gibbs_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the Gibbs free energy [J]. Property specific to solutions.
        """
        return self.solution_model.gibbs_hessian(self.pressure,
                                                 self.temperature,
                                                 self.molar_fractions)

    @material_property
    def entropy_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the entropy [J/K]. Property specific to solutions.
        """
        return self.solution_model.entropy_hessian(self.pressure,
                                                   self.temperature,
                                                   self.molar_fractions)

    @material_property
    def volume_hessian(self):
        """
        Returns an array containing the second compositional derivative
        of the volume [m^3]. Property specific to solutions.
        """
        return self.solution_model.volume_hessian(self.pressure,
                                                  self.temperature,
                                                  self.molar_fractions)

    @material_property
    def molar_gibbs(self):
        """
        Returns molar Gibbs free energy of the solution [J/mol].
        Aliased with self.gibbs.
        """
        return sum([self.endmembers[i][0].gibbs * self.molar_fractions[i]
                    for i in range(self.n_endmembers)]) + self.excess_gibbs

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
        return sum([self.endmembers[i][0].molar_mass
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)])

    @material_property
    def excess_volume(self):
        """
        Returns excess molar volume of the solution [m^3/mol].
        Specific property for solutions.
        """
        return self.solution_model.excess_volume(self.pressure,
                                                 self.temperature,
                                                 self.molar_fractions)

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the solution [m^3/mol].
        Aliased with self.V.
        """
        return sum([self.endmembers[i][0].molar_volume
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)]) + self.excess_volume

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
        return self.solution_model.excess_entropy(self.pressure,
                                                  self.temperature,
                                                  self.molar_fractions)

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropy of the solution [J/K/mol].
        Aliased with self.S.
        """
        return sum([self.endmembers[i][0].S * self.molar_fractions[i]
                    for i in range(self.n_endmembers)]) + self.excess_entropy

    @material_property
    def excess_enthalpy(self):
        """
        Returns excess molar enthalpy [J/mol].
        Property specific to solutions.
        """
        return self.solution_model.excess_enthalpy(self.pressure,
                                                   self.temperature,
                                                   self.molar_fractions)

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the solution [J/mol].
        Aliased with self.H.
        """
        return sum([self.endmembers[i][0].H * self.molar_fractions[i]
                    for i in range(self.n_endmembers)]) + self.excess_enthalpy

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the solution [Pa].
        Aliased with self.K_T.
        """
        return self.V * 1. / (sum([self.endmembers[i][0].V
                                   / (self.endmembers[i][0].K_T)
                                   * self.molar_fractions[i]
                                   for i in range(self.n_endmembers)]))

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the solution [Pa].
        Aliased with self.K_S.
        """
        if self.temperature < 1e-10:
            return self.isothermal_bulk_modulus
        else:
            return (self.isothermal_bulk_modulus
                    * self.molar_heat_capacity_p / self.molar_heat_capacity_v)

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the solution.
        (or inverse isothermal bulk modulus) [1/Pa].
        Aliased with self.K_T.
        """
        return 1. / self.isothermal_bulk_modulus

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the solution.
        (or inverse adiabatic bulk modulus) [1/Pa].
        Aliased with self.K_S.
        """
        return 1. / self.adiabatic_bulk_modulus

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the solution [Pa].
        Aliased with self.G.
        """
        G_list = np.fromiter((e[0].G for e in self.endmembers), dtype=float,
                             count=self.n_endmembers)
        return reuss_average_function(self.molar_fractions, G_list)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the solution [m/s].
        Aliased with self.v_p.
        """
        return np.sqrt((self.adiabatic_bulk_modulus
                        + 4. / 3. * self.shear_modulus) / self.density)

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the solution [m/s].
        Aliased with self.v_phi.
        """
        return np.sqrt(self.adiabatic_bulk_modulus / self.density)

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
            return float('nan')
        else:
            return (self.thermal_expansivity * self.isothermal_bulk_modulus
                    * self.molar_volume / self.molar_heat_capacity_v)

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha)
        of the solution [1/K].
        Aliased with self.alpha.
        """
        return (1. / self.V) * sum([self.endmembers[i][0].alpha
                                    * self.endmembers[i][0].V
                                    * self.molar_fractions[i]
                                    for i in range(self.n_endmembers)])

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the
        solution [J/K/mol].
        Aliased with self.C_v.
        """
        return (self.molar_heat_capacity_p
                - self.molar_volume * self.temperature
                * self.thermal_expansivity * self.thermal_expansivity
                * self.isothermal_bulk_modulus)

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure
        of the solution [J/K/mol].
        Aliased with self.C_p.
        """
        return sum([self.endmembers[i][0].molar_heat_capacity_p
                    * self.molar_fractions[i]
                    for i in range(self.n_endmembers)])

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
        reaction_basis = np.array([v[:] for v in
                                   self.stoichiometric_matrix.T.nullspace()])

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
        return [i for i in range(len(self.elements))
                if i not in self.independent_element_indices]

    @cached_property
    def compositional_null_basis(self):
        """
        An array N such that N.b = 0 for all bulk compositions that can
        be produced with a linear sum of the endmembers in the solution.
        """
        null_basis = np.array([v[:] for v in
                               self.stoichiometric_matrix.nullspace()])

        M = null_basis[:, self.dependent_element_indices]
        assert (M.shape[0] == M.shape[1]) and (M == np.eye(M.shape[0])).all()

        return null_basis

    @cached_property
    def endmember_formulae(self):
        """
        A list of formulae for all the endmember in the solution.
        """
        return [mbr[0].params['formula'] for mbr in self.endmembers]

    @cached_property
    def endmember_names(self):
        """
        A list of names for all the endmember in the solution.
        """
        return [mbr[0].name for mbr in self.endmembers]

    @cached_property
    def n_endmembers(self):
        """
        The number of endmembers in the solution.
        """
        return len(self.endmembers)

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
