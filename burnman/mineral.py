# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function
import warnings

import numpy as np

from .material import Material
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
        if 'params' not in self.__dict__:
            self.params = {}
        self.method = None
        if 'equation_of_state' in self.params:
            self.set_method(self.params['equation_of_state'])
        if 'name' in self.params:
            self.name=self.params['name']
        class_items = Mineral.__dict__.iteritems()
        for k, v in class_items:
            if isinstance(v, property):
                setattr(self,'_'+k,None)

            
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
        self.pressure = pressure
        self.temperature = temperature

        if self.method is None:
            raise AttributeError("no method set for mineral, or equation_of_state given in mineral.params")



    def composition(self):
        return self.params['formula']
    
    @property
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        """
        if self._molar_mass is None:
            self._molar_mass = self.params['molar_mass']
        return self._molar_mass

    @property
    def density(self):
        """
        Returns density of the mineral [kg/m^3]
        """
        if self._density is None:
            self._density = self.molar_mass/self.molar_volume
        return self._density

    @property
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        """
        if self._molar_volume is None:
            self._molar_volume=self.method.volume(self.pressure, self.temperature, self.params)
        return self._molar_volume

    @property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
        """
        if self._grueneisen_parameter is None:
            self._grueneisen_parameter = self.method.grueneisen_parameter(self.pressure, self.temperature, self.molar_volume, self.params)
        return self._grueneisen_parameter

    @property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        """
        if self._isothermal_bulk_modulus is None:
            self._isothermal_bulk_modulus = self.method.isothermal_bulk_modulus(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._isothermal_bulk_modulus

    @property
    def compressibility(self):
        """
        Returns compressibility of the mineral (or inverse isothermal bulk modulus) [1/Pa]
        """
        if self._compressibility is None:
            self._compressibility = 1./self.isothermal_bulk_modulus
        return self._compressibility

    @property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        """
        if self._adiabatic_bulk_modulus is None:
            self._adiabatic_bulk_modulus = self.method.adiabatic_bulk_modulus(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._adiabatic_bulk_modulus

    @property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        """
        if self._shear_modulus is None:
            self._shear_modulus = self.method.shear_modulus(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._shear_modulus

    @property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient (alpha) of the mineral [1/K]
        """
        if self._thermal_expansivity is None:
            self._thermal_expansivity = self.method.thermal_expansivity(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._thermal_expansivity

    @property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        """
        if self._heat_capacity_v is None :
            self._heat_capacity_v = self.method.heat_capacity_v(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._heat_capacity_v

    @property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        """
        if self._heat_capacity_p is None :
            self._heat_capacity_p = self.method.heat_capacity_p(self.pressure, self.temperature, self.molar_volume , self.params)
        return self._heat_capacity_p

    @property
    def v_s(self):
        """
        Returns shear wave speed of the mineral [m/s]
        """
        if self._v_s is None:
            self._v_s = np.sqrt(self.shear_modulus / self.density)
        return self._v_s

    @property
    def v_p(self):
        """
        Returns P wave speed of the mineral [m/s]
        """
        if self._v_p is None:
            self._v_p = np.sqrt((self.adiabatic_bulk_modulus + 4. / 3. * \
                             self.shear_modulus) / self.density)
        return self._v_p

    @property
    def v_phi(self):
        """
        Returns bulk sound speed of the mineral [m/s]
        """
        if self._v_phi is None:
            self._v_phi = np.sqrt(self.adiabatic_bulk_modulus / self.density)
        return self._v_phi
    
    @property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        """
        if self._molar_gibbs is None:
            self._molar_gibbs = self.method.gibbs_free_energy(self.pressure, self.temperature, self.molar_volume, self.params)
        return self._molar_gibbs

    @property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        """
        if self._molar_helmholtz is None:
            self._molar_helmholtz = self.method.helmholtz_free_energy(self.pressure, self.temperature, self.molar_volume, self.params)
        return self._molar_helmholtz

    @property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        if self._molar_enthalpy is None:
            self._molar_enthalpy = self.method.enthalpy(self.pressure, self.temperature, self.molar_volume, self.params)
        return self._molar_enthalpy

    @property
    def molar_entropy(self):
        """
        Returns enthalpy of the mineral [J]
        """
        if self._molar_entropy is None :
            self._molar_entropy = self.method.entropy(self.pressure, self.temperature, self.molar_volume, self.params)
        return self._molar_entropy

    gibbs = molar_gibbs
    V = molar_volume
    H = molar_enthalpy
    S  = molar_entropy
    C_p = heat_capacity_p
    C_v = heat_capacity_v
    alpha = thermal_expansivity
    K_T = isothermal_bulk_modulus
    K_S = adiabatic_bulk_modulus
    gr = grueneisen_parameter
    G = shear_modulus
