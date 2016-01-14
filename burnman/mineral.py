# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function
import warnings

import numpy as np

from .material import Material, material_property
from . import eos
from .tools import copy_documentation



class Mineral(Material):
    """
    This is the base class for all minerals. States of the mineral
    can only be queried after setting the pressure and temperature
    using set_state(). The method for computing properties of
    the material is set using set_method(). This is done during 
    initialisation if the param 'equation_of_state' has been defined.
    The method can be overridden later by the user.

    This class is available as ``burnman.Mineral``.

    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the material parameters are expected to be in plain SI units.  This
    means that the elastic moduli should be in Pascals and NOT Gigapascals,
    and the Debye temperature should be in K not C.  Additionally, the
    reference volume should be in m^3/(mol molecule) and not in unit cell
    volume and 'n' should be the number of atoms per molecule.  Frequently in
    the literature the reference volume is given in Angstrom^3 per unit cell.
    To convert this to m^3/(mol of molecule) you should multiply by 10^(-30) *
    N_a / Z, where N_a is Avogadro's number and Z is the number of formula units per
    unit cell. You can look up Z in many places, including www.mindat.org
    """

    def __init__(self):
        Material.__init__(self)
        if 'params' not in self.__dict__:
            self.params = {}
        if 'property_modifiers' not in self.__dict__:
            self.property_modifiers = []
        self.method = None
        if 'equation_of_state' in self.params:
            self.set_method(self.params['equation_of_state'])
        if 'name' in self.params:
            self.name = self.params['name']




    def set_method(self, equation_of_state):
        """
        Set the equation of state to be used for this mineral.
        Takes a string corresponding to any of the predefined
        equations of state:  'bm2', 'bm3', 'mgd2', 'mgd3', 'slb2', 'slb3',
        'mt', 'hp_tmt', or 'cork'.  Alternatively, you can pass a user defined
        class which derives from the equation_of_state base class.
        After calling set_method(), any existing derived properties
        (e.g., elastic parameters or thermodynamic potentials) will be out
        of date, so set_state() will need to be called again.
        """

        if equation_of_state is None:
            self.method = None
            return

        new_method = eos.create(equation_of_state)
        if self.method is not None and 'equation_of_state' in self.params:
            self.method = eos.create(self.params['equation_of_state'])

        if type(new_method).__name__ == 'instance':
            raise Exception("Please derive your method from object (see python old style classes)")

        if self.method is not None and type(new_method) is not type(self.method):

            # Warn user that they are changing the EoS
            warnings.warn('Warning, you are changing the method to ' + new_method.__class__.__name__ +
                          ' even though the material is designed to be used with the method ' +
                          self.method.__class__.__name__ +
                          '.  This does not overwrite any mineral attributes', stacklevel=2)
            self.reset()

        self.method = new_method

        #Validate the params object on the requested EOS. 
        try:
            self.method.validate_parameters(self.params)
        except Exception as e:
            print('Mineral ' + self.to_string() + ' failed to validate parameters with message : \" ' + e.message + '\"')
            raise

        # Invalidate the cache upon resetting the method
        self.reset()

    def to_string(self):
        """
        Returns the name of the mineral class
        """
        return "'" + self.__class__.__module__.replace(".minlib_",".") + "." + self.__class__.__name__ + "'"

    def debug_print(self, indent=""):
        print("%s%s" % (indent, self.to_string()))

    def unroll(self):
        return ([self], [1.0])

    @copy_documentation(Material.set_state)
    def set_state(self, pressure, temperature):
        Material.set_state(self, pressure, temperature)
        self._property_modifiers = eos.property_modifiers.calculate_property_modifications(self)

        if self.method is None:
             raise AttributeError("no method set for mineral, or equation_of_state given in mineral.params")

    
    """
    Properties from equations of state
    We choose the P, T properties (e.g. Gibbs(P, T) rather than Helmholtz(V, T)), 
    as it allows us to more easily apply corrections to the free energy
    """
    @material_property
    @copy_documentation(Material.molar_gibbs)
    def molar_gibbs(self):
        return self.method.gibbs_free_energy(self.pressure, self.temperature, self.molar_volume, self.params) \
            + self._property_modifiers['G']

    
    @material_property
    @copy_documentation(Material.molar_volume)
    def molar_volume(self):
        return self.method.volume(self.pressure, self.temperature, self.params) \
            + self._property_modifiers['dGdP']


    @material_property
    @copy_documentation(Material.molar_entropy)
    def molar_entropy(self):
        return self.method.entropy(self.pressure, self.temperature, self.molar_volume, self.params) \
            - self._property_modifiers['dGdT']

    @material_property
    @copy_documentation(Material.isothermal_bulk_modulus)
    def isothermal_bulk_modulus(self):
        V_orig = self.method.volume(self.pressure, self.temperature, self.params)
        K_T_orig = self.method.isothermal_bulk_modulus(self.pressure, self.temperature,
                                                  self.molar_volume, self.params)

        return self.molar_volume \
            / ((V_orig/K_T_orig) - self._property_modifiers['d2GdP2'])
    
    @material_property
    @copy_documentation(Material.heat_capacity_p)
    def heat_capacity_p(self):
        return self.method.heat_capacity_p(self.pressure, self.temperature,
                                           self.molar_volume , self.params) \
            - self.temperature*self._property_modifiers['d2GdT2']

    
    @material_property
    @copy_documentation(Material.thermal_expansivity)
    def thermal_expansivity(self):
        return ((self.method.thermal_expansivity(self.pressure, self.temperature,
                                                 self.molar_volume , self.params) \
                 * self.method.volume(self.pressure, self.temperature, self.params)) \
                + self._property_modifiers['d2GdPdT']) / self.molar_volume

    
    @material_property
    @copy_documentation(Material.shear_modulus)
    def shear_modulus(self):
        return self.method.shear_modulus(self.pressure, self.temperature, self.molar_volume , self.params)
    
    
    """
    Properties from mineral parameters,
    Legendre transformations 
    or Maxwell relations
    """
    
    @material_property
    @copy_documentation(Material.molar_mass)
    def molar_mass(self):
        if 'molar_mass' in self.params:
            return self.params['molar_mass']
        else:
            raise ValueError("No molar_mass parameter for mineral "+self.to_string+".")

    @material_property
    @copy_documentation(Material.density)
    def density(self):
        return self.molar_mass/self.molar_volume

    
    @material_property
    @copy_documentation(Material.internal_energy)
    def internal_energy(self):
        return self.molar_gibbs - self.pressure*self.molar_volume + self.temperature*self.molar_entropy

    
    @material_property
    @copy_documentation(Material.molar_helmholtz)
    def molar_helmholtz(self):
        return self.molar_gibbs - self.pressure*self.molar_volume

    
    @material_property
    @copy_documentation(Material.molar_enthalpy)
    def molar_enthalpy(self):
        return self.molar_gibbs + self.temperature*self.molar_entropy

    @material_property
    @copy_documentation(Material.adiabatic_bulk_modulus)
    def adiabatic_bulk_modulus(self):
        if self.temperature < 1.e-10:
            return self.isothermal_bulk_modulus
        else:
            return self.isothermal_bulk_modulus * self.heat_capacity_p / self.heat_capacity_v

    @material_property
    @copy_documentation(Material.isothermal_compressibility)
    def isothermal_compressibility(self):
        return 1./self.isothermal_bulk_modulus

    @material_property
    @copy_documentation(Material.adiabatic_compressibility)
    def adiabatic_compressibility(self):
        return 1./self.adiabatic_bulk_modulus

    
    @material_property
    @copy_documentation(Material.p_wave_velocity)
    def p_wave_velocity(self):
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
                             self.shear_modulus) / self.density)

    @material_property
    @copy_documentation(Material.bulk_sound_velocity)
    def bulk_sound_velocity(self):
        return np.sqrt(self.adiabatic_bulk_modulus / self.density)

    @material_property
    @copy_documentation(Material.shear_wave_velocity)
    def shear_wave_velocity(self):
        return np.sqrt(self.shear_modulus / self.density)

    @material_property
    @copy_documentation(Material.grueneisen_parameter)
    def grueneisen_parameter(self):
        if self.temperature < 1.e-12:
            return 0.
        else:
            return self.thermal_expansivity * self.isothermal_bulk_modulus \
                * self.molar_volume / self.heat_capacity_v

    @material_property
    @copy_documentation(Material.heat_capacity_v)
    def heat_capacity_v(self):
        return self.heat_capacity_p - self.molar_volume * self.temperature \
            * self.thermal_expansivity * self.thermal_expansivity \
            * self.isothermal_bulk_modulus
