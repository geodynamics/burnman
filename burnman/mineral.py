# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function
import warnings

import numpy as np

from .material import Material, material_property
from . import eos



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
            warnings.warn('Warning, you are changing the method to ' + new_method.__class__.__name__ + ' even though the material is designed to be used with the method ' + self.method.__class__.__name__ + '.  This does not overwrite any mineral attributes other than temperature and pressure, which are set to nan. Stale attributes will be preserved unless they are refreshed (for example, by set_state).', stacklevel=2)

            try:
                self.pressure=self.temperature=float('nan')
            except AttributeError:
                pass

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

    def eos_pressure(self, temperature, volume):
        return self.method.pressure(temperature, volume, self.params)

    def set_state(self, pressure, temperature):
        """
        Update the material to the given pressure [Pa] and temperature [K].

        This updates the other properties of this class (v_s, v_p, ...).
        """
        Material.set_state(self, pressure, temperature)

        if self.method is None:
             raise AttributeError("no method set for mineral, or equation_of_state given in mineral.params")

    
    @material_property
    def internal_energy(self):
        """
        Returns internal energy of the mineral [J]
        Aliased with self.energy
        """
        return self.method.internal_energy(self.pressure, self.temperature, self.molar_volume, self.params )

    
    @material_property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        Aliased with self.gibbs
        """
        return self.method.gibbs_free_energy(self.pressure, self.temperature, self.molar_volume, self.params)


    @material_property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        Aliased with self.helmholtz
        """
        return self.method.helmholtz_free_energy(self.pressure, self.temperature, self.molar_volume, self.params)
    
    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        if 'molar_mass' in self.params:
            return self.params['molar_mass']
        else:
            raise ValueError("No molar_mass parameter for mineral "+self.to_string+".")



    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        Aliased with self.V
        """
        return self.method.volume(self.pressure, self.temperature, self.params)
  

    @material_property
    def density(self):
        """
        Returns density of the mineral [kg/m^3]
        Aliased with self.rho
        """
        return self.molar_mass/self.molar_volume


    @material_property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.S
        """
        return self.method.entropy(self.pressure, self.temperature, self.molar_volume, self.params)


    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.H
        """
        return self.method.enthalpy(self.pressure, self.temperature, self.molar_volume, self.params)



    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        Aliased with self.K_T
        """
        return  self.method.isothermal_bulk_modulus(self.pressure, self.temperature, self.molar_volume , self.params)



    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        """
        return self.method.adiabatic_bulk_modulus(self.pressure, self.temperature, self.molar_volume , self.params)

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the mineral (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_T
        """
        return 1./self.isothermal_bulk_modulus


    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the mineral (or inverse adiabatic bulk modulus) [1/Pa]
        Aliased with self.beta_S
        """
        return 1./self.adiabatic_bulk_modulus


    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        """
        return self.method.shear_modulus(self.pressure, self.temperature, self.molar_volume , self.params)

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the mineral [m/s]
        Aliased with self.v_p
        """
        return np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
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
        return self.method.grueneisen_parameter(self.pressure, self.temperature, self.molar_volume, self.params)



    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the mineral [1/K]
        Aliased with self.alpha
        """
        return self.method.thermal_expansivity(self.pressure, self.temperature, self.molar_volume , self.params)


    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        Aliased with self.C_v
        """
        return self.method.heat_capacity_v(self.pressure, self.temperature, self.molar_volume , self.params)


    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        Aliased with self.C_p
        """
        return self.method.heat_capacity_p(self.pressure, self.temperature, self.molar_volume , self.params)





