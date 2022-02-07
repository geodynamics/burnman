# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np

from .mineral import Mineral, material_property
from .solution import Solution


class CombinedMineral(Mineral):

    """
    This is the base class for endmembers constructed from a
    linear combination of other minerals.

    Instances of this class should be initialised with a
    list of Mineral instances, a second list containing the
    number of moles of each mineral, and (optionally) a
    third list containing three floats describing a
    free energy adjustment which is linear in pressure and temperature
    (i.e. a constant energy [J/mol], entropy [J/K/mol]
    and volume adjustment [J/Pa/mol or m^3/mol]).

    For example, a crude approximation to a bridgmanite model might be
    bdg = CombinedMineral([per, stv], [1.0, 1.0], [-15.e3, 0., 0.])

    This class is available as :class:`burnman.CombinedMineral`.
    """

    def __init__(self, mineral_list, molar_amounts,
                 free_energy_adjustment=[],
                 name='User-created endmember'):
        self.mixture = Solution(solution_type='mechanical',
                                endmembers=[[m, ''] for m in mineral_list],
                                molar_fractions=molar_amounts)

        # Remove elements from the chemical formula if they have
        # negligible concentrations
        for key, value in list(self.mixture.formula.items()):
            if np.abs(value) < 1.e-10:
                self.mixture.formula.pop(key)

        self.params = {'name': name,
                       'formula': self.mixture.formula,
                       'equation_of_state': 'combined',
                       'molar_mass': self.mixture.molar_mass,
                       'n': sum(self.mixture.formula.values())}

        if free_energy_adjustment != []:
            assert(len(free_energy_adjustment) == 3)
            dE, dS, dV = free_energy_adjustment
            self.property_modifiers = [['linear', {'delta_E': dE,
                                                   'delta_S': dS,
                                                   'delta_V': dV}]]

        Mineral.__init__(self)

    def set_state(self, pressure, temperature):
        self.mixture.set_state(pressure, temperature)
        Mineral.set_state(self, pressure, temperature)

    @material_property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        Aliased with self.gibbs
        """
        return self.mixture.molar_gibbs + self._property_modifiers['G']

    @material_property
    def _molar_volume_unmodified(self):
        return self.mixture.molar_volume

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        Aliased with self.V
        """
        return self.mixture.molar_volume + self._property_modifiers['dGdP']

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the mineral [J]
        Aliased with self.S
        """
        return self.mixture.molar_entropy - self._property_modifiers['dGdT']

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        Aliased with self.K_T
        """
        K_T_orig = self.mixture.isothermal_bulk_modulus

        return (self.molar_volume
                / ((self._molar_volume_unmodified / K_T_orig)
                   - self._property_modifiers['d2GdP2']))

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        """
        return self.mixture.shear_modulus

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the mineral [1/K]
        Aliased with self.alpha
        """
        return ((self.mixture.thermal_expansivity
                 * self._molar_volume_unmodified)
                + self._property_modifiers['d2GdPdT']) / self.molar_volume

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        Aliased with self.C_p
        """
        return (self.mixture.molar_heat_capacity_p
                - self.temperature * self._property_modifiers['d2GdT2'])

    """
    Properties from mineral parameters,
    Legendre transformations
    or Maxwell relations
    """

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        return self.mixture.molar_mass

    @material_property
    def formula(self):
        """
        Returns molar chemical formula of the mineral
        """
        return self.mixture.formula

    @material_property
    def density(self):
        """
        Returns density of the mineral [kg/m^3]
        Aliased with self.rho
        """
        return self.molar_mass / self.molar_volume

    @material_property
    def molar_internal_energy(self):
        """
        Returns molar internal energy of the mineral [J/mol]
        Aliased with self.energy
        """
        return (self.molar_gibbs
                - self.pressure * self.molar_volume
                + self.temperature * self.molar_entropy)

    @material_property
    def molar_helmholtz(self):
        """
        Returns molar Helmholtz free energy of the mineral [J/mol]
        Aliased with self.helmholtz
        """
        return self.molar_gibbs - self.pressure * self.molar_volume

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the mineral [J/mol]
        Aliased with self.H
        """
        return self.molar_gibbs + self.temperature * self.molar_entropy

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        """
        if self.temperature < 1.e-10:
            return self.isothermal_bulk_modulus
        else:
            return (self.isothermal_bulk_modulus
                    * self.molar_heat_capacity_p
                    / self.molar_heat_capacity_v)

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the mineral
        (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.K_T
        """
        return 1. / self.isothermal_bulk_modulus

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the mineral
        (or inverse adiabatic bulk modulus) [1/Pa]
        Aliased with self.K_S
        """
        return 1. / self.adiabatic_bulk_modulus

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the mineral [m/s]
        Aliased with self.v_p
        """
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. *
                        self.shear_modulus) / self.density)

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the mineral [m/s]
        Aliased with self.v_phi
        """
        return np.sqrt(self.adiabatic_bulk_modulus / self.density)

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the mineral [m/s]
        Aliased with self.v_s
        """
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
        Aliased with self.gr
        """
        if self.temperature < 1.e-12:
            return 0.
        else:
            return self.thermal_expansivity * self.isothermal_bulk_modulus \
                * self.molar_volume / self.molar_heat_capacity_v

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the mineral [J/K/mol]
        Aliased with self.C_v
        """
        return (self.molar_heat_capacity_p
                - self.molar_volume * self.temperature
                * self.thermal_expansivity * self.thermal_expansivity
                * self.isothermal_bulk_modulus)
