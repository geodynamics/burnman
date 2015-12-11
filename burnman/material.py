from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

import numpy as np
import copy

class Material(object):
    """
    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type :class:`~burnman.mineral.Mineral` and their molar
    fractions. This class is available as ``burnman.Material``.

    The user needs to call set_method() (once in the beginning) and set_state()
    before querying the material with unroll() or density().

    Attributes
    ----------
    pressure : float
        The current pressure as set by :func:`~burnman.Material.set_state`. [Pa]
    temperature : float
        The current temperature as set by :func:`~burnman.Material.set_state`. [K]
    """

    def __init__(self):
        self._pressure = None
        self._temperature = None

    def set_method(self, method):
        """
        Set the averaging method. See :doc:`averaging` for details.

        Notes
        -----
        Needs to be implemented in derived classes.
        """
        raise NotImplementedError("need to implement set_method() in derived class!")

    def to_string(self):
        """
        Returns a human-readable name of this material. The default implementation will return the name of the class,
        which is a reasonable default.

        Returns
        -------
        name : string
            Name of this material.
        """
        return "'" + self.__class__.__name__ + "'"

    def debug_print(self, indent=""):
        """
        Print a human-readable representation of this Material.
        """
        raise NotImplementedError("Derived classes need to implement debug_print(). This is '" + self.__class__.__name__ + "'")

    def print_minerals_of_current_state(self):
        """
        Print a human-readable representation of this Material at the current P, T as a list of minerals.
        This requires set_state() has been called before.
        """
        (minerals, fractions) = self.unroll()
        if len(minerals)==1:
            print minerals[0].to_string()
        else:
            print "Material %s:" % self.to_string()
            for (mineral, fraction) in zip(minerals, fractions):
                print "  %g of phase %s" % (fraction, mineral.to_string())


    def set_state(self, pressure, temperature):
        """
        Set the material to the given pressure and temperature.

        Parameters
        ----------
        pressure : float
            The desired pressure in [Pa].
        temperature : float
            The desired temperature in [K].
        """
        self._pressure = pressure
        self._temperature = temperature

    def unroll(self):
        """
        Unroll this material into a list of :class:`burnman.Mineral` and their molar fractions. All averaging schemes
        then operate on this list of minerals. Note that the return value of this function may depend on the current
        state (temperature, pressure).

        Notes
        -----
        Needs to be implemented in derived classes.

        Returns
        -------
        fractions : list of float
            List of molar fractions, should sum to 1.0.
        minerals : list of :class:`burnman.Mineral`
            List of minerals.
        """
        raise NotImplementedError("need to implement unroll() in derived class!")
        return ([], [])


    def evaluate(self,vars_list,pressures, temperatures):
        """
        Returns an array of material properties requested through a list of strings at given pressure and temperature conditions. At the end it resets the set_state to the original values.
        The user needs to call set_method() before.
        
            
        Parameters
        ----------
        vars_list : list of strings
        Variables to be returned for given conditions
        pressure : array of float
        Array of pressures in [Pa].
        temperature : float
        Array of temperaturesin [K].
        
        Returns
        -------
        
        """
        old_pressure = self.pressure
        old_temperature = self.temperature
        output = np.empty((len(vars_list),len(pressures)))
        for i in range(len(pressures)):
            self.set_state(pressures[i],temperatures[i])
            for j in range(len(vars_list)):
                output[j,i]=getattr(self,vars_list[j])
        self.set_state(old_pressure,old_temperature)
        return output


    @property
    def pressure (self):
        """
        Returns pressure set in state [Pa]
        Aliased with self.P
        
        Returns
        -------
        pressure : float
        
        """
        return self._pressure
    
    @property
    def temperature (self):
        """
        Returns pressure set in state [K]
        Aliased with self.T
            
        Returns
        -------
        pressure : float
        
        """
        return self._temperature
 
    @property
    def internal_energy(self):
        """
        Returns internal energy of the mineral [J]
        Aliased with self.energy
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        internal energies : float
        """
        raise NotImplementedError("need to implement internal_energy() in derived class!")
        return None
    
    @property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        Aliased with self.gibbs
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        Gibbs free energies : float
        """
        raise NotImplementedError("need to implement molar_gibbs() in derived class!")
        return None
    
    @property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        Aliased with self.helholtz
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        Helmholtz free energies : float
        """
        raise NotImplementedError("need to implement molar_helmholtz() in derived class!")
        return None

    @property
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        molar_mass : float
        
        """
        raise NotImplementedError("need to implement molar_mass() in derived class!")
        return None
    
    
    @property
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        Aliased with self.V
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        molar_volume : float
        """
        raise NotImplementedError("need to implement molar_volume() in derived class!")
        return None

    @property
    def density(self):
        """
        Returns the density of this material.
        Aliased with self.rho
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        density : float
        The density of this material in [kg/m^3]
        
        """
        raise NotImplementedError("need to implement density() in derived class!")
        return None


    @property
    def molar_entropy(self):
        """
        Returns entropy of the mineral [J]
        Aliased with self.S
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        entropies : float
        """
        raise NotImplementedError("need to implement molar_entropy() in derived class!")
        return None


    @property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        Aliased with self.H
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        enthalpies : float
        """
        raise NotImplementedError("need to implement molar_enthalpy() in derived class!")
        return None
    
    
    @property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        Aliased with self.K_T
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        isothermal_bulk_modulus : float
        """
        raise NotImplementedError("need to implement isothermal_bulk_moduls() in derived class!")
        return None

    @property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        Aliased with self.K_S
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        adiabatic_bulk_modulus : array of floats
        """
        raise NotImplementedError("need to implement adiabatic_bulk_modulus() in derived class!")
        return None
    
    
    @property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the mineral (or inverse isothermal bulk modulus) [1/Pa]
        Aliased with self.beta_T
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        (K_T)^-1 : float
        
        """
        raise NotImplementedError("need to implement compressibility() in derived class!")
        return None
    
    
    @property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the mineral (or inverse adiabatic bulk modulus) [1/Pa]
        Aliased with self.beta_S
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        (K_S)^-1 : float
        
        """
        raise NotImplementedError("need to implement compressibility() in derived class!")
        return None

    @property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        Aliased with self.G
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        shear_modulus : float
        """
        raise NotImplementedError("need to implement shear_modulus() in derived class!")
        return None
    
    @property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the mineral [m/s]
        Aliased with self.v_p
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        p_wave_velocity : float
        """
        raise NotImplementedError("need to implement p_wave_velocity() in derived class!")
        return None
    
    @property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the mineral [m/s]
        Aliased with self.v_phi
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        bulk sound velocity: float
        """
        raise NotImplementedError("need to implement bulk_sound_velocity() in derived class!")
        return None
    
    @property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the mineral [m/s]
        Aliased with self.v_s
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        shear wave velociy : float
        """
        raise NotImplementedError("need to implement shear_wave_velocity() in derived class!")
        return None


    @property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
        Aliased with self.gr
                
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        gr : float
        """
        raise NotImplementedError("need to implement grueneisen_parameter() in derived class!")
        return None
    
    
    @property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the mineral [1/K]
        Aliased with self.alpha
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        alpha : float
        """
        raise NotImplementedError("need to implement thermal_expansivity() in derived class!")
        return None


    @property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        Aliased with self.C_v
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        heat_capacity_v : float
        """
        raise NotImplementedError("need to implement heat_capacity_v() in derived class!")
        return None

    @property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        Aliased with self.C_p
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        heat_capacity_p : float
        """
        raise NotImplementedError("need to implement heat_capacity_p() in derived class!")
        return None


    def copy(self):
        return copy.deepcopy(self)








    ########################################################################
    # Aliased properties
    
    @property
    def P(self):
        return self.pressure
    @property
    def T(self):
        return self.temperature
    @property
    def energy(self):
        return self.internal_energy
    @property
    def helmholtz(self):
        return self.molar_helmholtz
    @property
    def gibbs(self):
        return self.molar_gibbs
    @property
    def V(self):
        return self.molar_volume
    @property
    def rho(self):
        return self.density
    @property
    def S(self):
        return self.molar_entropy
    @property
    def H(self):
        return self.molar_enthalpy
    @property
    def K_T(self):
        return self.isothermal_bulk_modulus
    @property
    def K_S(self):
        return self.adiabatic_bulk_modulus
    @property
    def beta_T(self):
        return self.isothermal_compressibility
    @property
    def beta_S(self):
        return self.adiabatic_compressibility
    @property
    def G(self):
        return self.shear_modulus
    @property
    def v_p(self):
        return self.p_wave_velocity
    @property
    def v_phi(self):
        return self.bulk_sound_velocity
    @property
    def v_s(self):
        return self.shear_wave_velocity
    @property
    def gr(self):
        return self.grueneisen_parameter
    @property
    def alpha(self):
        return self.thermal_expansivity
    @property
    def C_v(self):
        return self.heat_capacity_v
    @property
    def C_p(self):
        return self.heat_capacity_p
    
















