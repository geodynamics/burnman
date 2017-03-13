# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import warnings

from .mineral import Mineral, material_property
from .solutionmodel import *
from .processchemistry import sum_formulae
from .averaging_schemes import reuss_average_function
from . import constants


class SolidSolution(Mineral):

    """
    This is the base class for all solid solutions.
    Site occupancies, endmember activities and the constant
    and pressure and temperature dependencies of the excess
    properties can be queried after using set_composition()
    States of the solid solution can only be queried after setting
    the pressure, temperature and composition using set_state().

    This class is available as :class:`burnman.SolidSolution`.
    It uses an instance of :class:`burnman.SolutionModel` to
    calculate interaction terms between endmembers.

    All the solid solution parameters are expected to be in SI units.  This
    means that the interaction parameters should be in J/mol, with the T
    and P derivatives in J/K/mol and m^3/mol.
    """

    def __init__(self,
                 name=None,
                 solution_type=None,
                 endmembers=None,
                 energy_interaction=None,
                 volume_interaction=None,
                 entropy_interaction=None,
                 alphas=None,
                 molar_fractions=None):
        """
        Set up matrices to speed up calculations for when P, T, X is defined.

        Parameters
        ----------
        endmembers: list of :class:`burnman.Mineral`
            List of endmembers in this solid solution.
        solution_model: :class:`burnman.SolutionModel`
            SolutionModel to use.
        """
        Mineral.__init__(self)


        class SolidSolutionMethod(object):

            """Dummy class because SolidSolution needs a method to call
            Mineral.set_state(), but should never have a method that
            is used for minerals. Note that set_method() below will
            not change self.method"""
            pass
        self.method = SolidSolutionMethod()

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
        if alphas is not None:
            self.alphas = alphas
        if endmembers is not None:
            self.endmembers = endmembers

        if hasattr(self, 'endmembers') == False:
            raise Exception(
                "'endmembers' attribute missing from solid solution")

        
        # Set default solution model type
        if hasattr(self, 'solution_type'):
            if self.solution_type == 'mechanical':
                self.solution_model = MechanicalSolution(self.endmembers)
            elif self.solution_type == 'ideal':
                self.solution_model = IdealSolution(self.endmembers)
            else:
                if hasattr(self, 'energy_interaction') == False:
                    self.energy_interaction = None
                if hasattr(self, 'volume_interaction') == False:
                    self.volume_interaction = None
                if hasattr(self, 'entropy_interaction') == False:
                    self.entropy_interaction = None

                if self.solution_type == 'symmetric':
                    self.solution_model = SymmetricRegularSolution(
                        self.endmembers, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
                elif self.solution_type == 'asymmetric':
                    try:
                        self.solution_model = AsymmetricRegularSolution(
                            self.endmembers, self.alphas, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
                    except:
                        raise Exception(
                            "'alphas' attribute missing from solid solution")
                elif self.solution_type == 'subregular':
                    self.solution_model = SubregularSolution(
                        self.endmembers, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
                else:
                    raise Exception(
                        "Solution model type " + self.solution_type + "not recognised.")
        else:
            if hasattr(self, 'energy_interaction') == False:
                self.energy_interaction = None
            if hasattr(self, 'volume_interaction') == False:
                self.volume_interaction = None
            if hasattr(self, 'entropy_interaction') == False:
                self.entropy_interaction = None

            if self.solution_type == 'symmetric':
                self.solution_model = SymmetricRegularSolution(
                    self.endmembers, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
            elif self.solution_type == 'asymmetric':
                try:
                    self.solution_model = AsymmetricRegularSolution(
                        self.endmembers, self.alphas, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
                except:
                    raise Exception(
                        "'alphas' attribute missing from solid solution")
            elif self.solution_type == 'subregular':
                self.solution_model = SubregularSolution(
                    self.endmembers, self.energy_interaction, self.volume_interaction, self.entropy_interaction)
            else:
                raise Exception(
                    "Solution model type " + self.solution_type + "not recognised.")
            self.solution_model = SolutionModel()

        # Number of endmembers in the solid solution
        self.n_endmembers = len(self.endmembers)

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
        Set the composition for this solid solution.

        Parameters
        ----------
        molar_fractions: list of float
            molar abundance for each endmember, needs to sum to one.
        """
        assert(len(self.endmembers) == len(molar_fractions))

        if self.solution_type != 'mechanical':
            assert(sum(molar_fractions) > 0.9999)
            assert(sum(molar_fractions) < 1.0001)
            
        self.molar_fractions = molar_fractions
        
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
        Returns chemical formula of the solid solution
        """
        return sum_formulae([self.endmembers[i][0].params['formula'] for i in range(self.n_endmembers)],
                            self.molar_fractions)
    
    @material_property
    def activities(self):
        """
        Returns a list of endmember activities [unitless]
        """
        return self.solution_model.activities(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def activity_coefficients(self):
        """
        Returns a list of endmember activity coefficients (gamma = activity / ideal activity) [unitless]
        """
        return self.solution_model.activity_coefficients(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def internal_energy(self):
        """
        Returns internal energy of the mineral [J]
        Aliased with self.energy
        """
        return self.molar_helmholtz - self.pressure * self.molar_volume

    @material_property
    def excess_partial_gibbs(self):
        """
        Returns excess partial gibbs free energy [J]
        Property specific to solid solutions.
        """
        return self.solution_model.excess_partial_gibbs_free_energies(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def partial_gibbs(self):
        """
        Returns excess partial gibbs free energy [J]
        Property specific to solid solutions.
        """
        return np.array([self.endmembers[i][0].gibbs for i in range(self.n_endmembers)]) + self.excess_partial_gibbs

    @material_property
    def excess_gibbs(self):
        """
        Returns excess gibbs free energy [J]
        Property specific to solid solutions.
        """
        return self.solution_model.excess_gibbs_free_energy(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the solid solution [J]
        Aliased with self.gibbs
        """
        return sum([self.endmembers[i][0].gibbs * self.molar_fractions[i] for i in range(self.n_endmembers)]) + self.excess_gibbs

    @material_property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the solid solution [J]
        Aliased with self.helmholtz
        """
        return self.molar_gibbs + self.temperature * self.molar_entropy

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the solid solution [kg/mol]
        """
        return sum([self.endmembers[i][0].molar_mass * self.molar_fractions[i] for i in range(self.n_endmembers)])
    
    @material_property
    def formula(self):
        """
        Returns chemical formula of the solid solution
        """
        return sum_formulae([self.endmembers[i][0].params['formula'] for i in range(self.n_endmembers)],
                            self.molar_fractions)

    @material_property
    def excess_volume(self):
        """
        Returns excess volume of the solid solution [m^3/mol]
        Specific property for solid solutions
        """
        return self.solution_model.excess_volume(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the solid solution [m^3/mol]
        Aliased with self.V
        """
        return sum([self.endmembers[i][0].molar_volume * self.molar_fractions[i] for i in range(self.n_endmembers)]) + self.excess_volume

    @material_property
    def density(self):
        """
        Returns density of the solid solution [kg/m^3]
        Aliased with self.rho
        """
        return self.molar_mass / self.molar_volume

    @material_property
    def excess_entropy(self):
        """
        Returns excess entropy [J]
        Property specific to solid solutions.
        """
        return self.solution_model.excess_entropy(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the solid solution [J]
        Aliased with self.S
        """
        return sum([self.endmembers[i][0].S * self.molar_fractions[i] for i in range(self.n_endmembers)]) + self.excess_entropy

    @material_property
    def excess_enthalpy(self):
        """
        Returns excess enthalpy [J]
        Property specific to solid solutions.
        """
        return self.solution_model.excess_enthalpy(self.pressure, self.temperature, self.molar_fractions)

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the solid solution [J]
        Aliased with self.H
        """
        return sum([self.endmembers[i][0].H * self.molar_fractions[i] for i in range(self.n_endmembers)]) + self.excess_enthalpy

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the solid solution [Pa]
        Aliased with self.K_T
        """
        return self.V * 1. / (sum([self.endmembers[i][0].V / (self.endmembers[i][0].K_T) * self.molar_fractions[i] for i in range(self.n_endmembers)]))

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the solid solution [Pa]
        Aliased with self.K_S
        """
        if self.temperature < 1e-10:
            return self.isothermal_bulk_modulus
        else:
            return self.isothermal_bulk_modulus * self.heat_capacity_p / self.heat_capacity_v

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the solid solution (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.K_T
        """
        return 1. / self.isothermal_bulk_modulus

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the solid solution (or inverse adiabatic bulk modulus) [1/Pa]
        Aliased with self.K_S
        """
        return 1. / self.adiabatic_bulk_modulus

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the solid solution [Pa]
        Aliased with self.G
        """
        G_list = np.fromiter(
            (e[0].G for e in self.endmembers), dtype=np.float, count=self.n_endmembers)
        return reuss_average_function(self.molar_fractions, G_list)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the solid solution [m/s]
        Aliased with self.v_p
        """
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. *
                        self.shear_modulus) / self.density)

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the solid solution [m/s]
        Aliased with self.v_phi
        """
        return np.sqrt(self.adiabatic_bulk_modulus / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the solid solution [m/s]
        Aliased with self.v_s
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the solid solution [unitless]
        Aliased with self.gr
        """
        if self.temperature < 1e-10:
            return float('nan')
        else:
            return self.thermal_expansivity * self.isothermal_bulk_modulus * self.molar_volume / self.heat_capacity_v

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the solid solution [1/K]
        Aliased with self.alpha
        """
        return (1. / self.V) * sum([self.endmembers[i][0].alpha * self.endmembers[i][0].V * self.molar_fractions[i] for i in range(self.n_endmembers)])

    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the solid solution [J/K/mol]
        Aliased with self.C_v
        """
        return self.heat_capacity_p - self.molar_volume * self.temperature * self.thermal_expansivity * self.thermal_expansivity * self.isothermal_bulk_modulus

    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the solid solution [J/K/mol]
        Aliased with self.C_p
        """
        return sum([self.endmembers[i][0].heat_capacity_p * self.molar_fractions[i] for i in range(self.n_endmembers)])
