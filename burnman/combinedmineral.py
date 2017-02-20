# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import warnings

from .mineral import Mineral, material_property
from .solidsolution import SolidSolution
from .solutionmodel import *
from .averaging_schemes import reuss_average_function
from .processchemistry import sum_formulae
from . import constants


class CombinedMineral(Mineral):

    """
    This is the base class for endmembers constructed from a 
    linear combination of other minerals.

    Instances of this class should be initialised with a 
    list of Mineral instances, a second list containing the 
    number of moles of each mineral, and (optionally) a 
    third list containing three floats describing a 
    free energy adjustment which is linear in pressure and temperature 
    (i.e. a constant energy, entropy and volume adjustment).

    For example, a crude approximation to a bridgmanite model might be
    bdg = CombinedMineral([per, stv], [1.0, 1.0], [-15.e3, 0., 0.])

    This class is available as :class:`burnman.CombinedMineral`.
    """

    def __init__(self, mineral_list, molar_amounts, free_energy_adjustment=[]):
        self.mixture = SolidSolution(solution_type='mechanical',
                                     endmembers=[[m, ''] for m in mineral_list],
                                     molar_fractions = molar_amounts)

        self.params = {
            'name': 'User-created endmember',
            'formula': self.mixture.formula,
            'molar_mass': self.mixture.molar_mass,
            'n': sum(self.mixture.formula.values())
        }

        if free_energy_adjustment != []:
            assert(len(free_energy_adjustment) == 3)
            self.property_modifiers = [['linear', {'delta_E': free_energy_adjustment[0],
                                                   'delta_S': free_energy_adjustment[1],
                                                   'delta_V': free_energy_adjustment[2]}]]

        Mineral.__init__(self)
        
        class MadeMemberMethod(object):

            """Dummy class because SolidSolution needs a method to call
            Mineral.set_state(), but should never have a method that
            is used for minerals. Note that set_method() below will
            not change self.method"""
            pass
        self.method = MadeMemberMethod()


    def set_state(self, pressure, temperature):
        self.mixture.set_state(pressure, temperature)
        Mineral.set_state(self, pressure, temperature)


    @material_property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the solid solution [J]
        Aliased with self.gibbs
        """
        return self.mixture.molar_gibbs + self._property_modifiers['G']

    @material_property
    def _molar_volume_unmodified(self):
        return self.mixture.molar_volume

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the solid solution [m^3/mol]
        Aliased with self.V
        """
        return self.mixture.molar_volume + self._property_modifiers['dGdP']

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the solid solution [J]
        Aliased with self.S
        """
        return self.mixture.molar_entropy - self._property_modifiers['dGdT']

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the solid solution [Pa]
        Aliased with self.K_T
        """
        K_T_orig = self.mixture.isothermal_bulk_modulus

        return self.molar_volume \
            / ((self._molar_volume_unmodified / K_T_orig) - self._property_modifiers['d2GdP2'])

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the solid solution [Pa]
        Aliased with self.G
        """
        return self.mixture.shear_modulus

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the solid solution [1/K]
        Aliased with self.alpha
        """
        return ((self.mixture.thermal_expansivity * self._molar_volume_unmodified)
                + self._property_modifiers['d2GdPdT']) / self.molar_volume
    
    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the solid solution [J/K/mol]
        Aliased with self.C_p
        """
        return self.mixture.heat_capacity_p - self.temperature * self._property_modifiers['d2GdT2']


    """
    Properties from mineral parameters,
    Legendre transformations
    or Maxwell relations
    """

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the solid solution [kg/mol]
        """
        return self.mixture.molar_mass
    
    @material_property
    def formula(self):
        """
        Returns chemical formula of the solid solution
        """
        return self.mixture.formula

    @material_property
    def density(self):
        """
        Returns density of the solid solution [kg/m^3]
        Aliased with self.rho
        """
        return self.molar_mass / self.molar_volume

    @material_property
    def internal_energy(self):
        """
        Returns internal energy of the mineral [J]
        Aliased with self.energy
        """
        return self.molar_gibbs - self.pressure * self.molar_volume + self.temperature * self.molar_entropy
    
    @material_property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the solid solution [J]
        Aliased with self.helmholtz
        """
        return self.molar_gibbs - self.pressure * self.molar_volume


    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the solid solution [J]
        Aliased with self.H
        """
        return self.molar_gibbs + self.temperature * self.molar_entropy


    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the solid solution [Pa]
        Aliased with self.K_S
        """
        if self.temperature < 1.e-10:
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
        if self.temperature < 1.e-12:
            return 0.
        else:
            return self.thermal_expansivity * self.isothermal_bulk_modulus \
                * self.molar_volume / self.heat_capacity_v

    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the solid solution [J/K/mol]
        Aliased with self.C_v
        """
        return self.heat_capacity_p - self.molar_volume * self.temperature \
            * self.thermal_expansivity * self.thermal_expansivity \
            * self.isothermal_bulk_modulus

