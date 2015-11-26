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
        self.pressure = None
        self.temperature = None

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
        self.pressure = pressure
        self.temperature = temperature

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
    def molar_mass(self):
        """
        Returns molar mass of the mineral [kg/mol]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        molar_mass : array of floats
        
        """
        raise NotImplementedError("need to implement molar_mass() in derived class!")
        return None

    @property
    def density(self):
        """
        Returns the density of this material.
        
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
    def molar_volume(self):
        """
        Returns molar volume of the mineral [m^3/mol]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        V : array of floats
        """
        raise NotImplementedError("need to implement molar_volume() in derived class!")
        return None

    @property
    def grueneisen_parameter(self):
        """
        Returns grueneisen parameter of the mineral [unitless]
                
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        gr : array of floats
        """
        raise NotImplementedError("need to implement grueneisen_parameter() in derived class!")
        return None

    @property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the mineral [Pa]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        K_T : array of floats
        """
        raise NotImplementedError("need to implement isothermal_bulk_moduls() in derived class!")
        return None

    @property
    def compressibility(self):
        """
        Returns compressibility of the mineral (or inverse isothermal bulk modulus) [1/Pa]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        (K_T)^-1 : array of floats
        
        """
        raise NotImplementedError("need to implement compressibility() in derived class!")
        return None

    @property
    def adiabatic_bulk_modulus(self):
        """
        Returns adiabatic bulk modulus of the mineral [Pa]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        K_S : array of floats
        """
        raise NotImplementedError("need to implement adiabatic_bulk_modulus() in derived class!")
        return None

    @property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral [Pa]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        G : array of floats
        """
        raise NotImplementedError("need to implement shear_modulus() in derived class!")
        return None

    @property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the mineral [1/K]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        alpha : array of floats
        """
        raise NotImplementedError("need to implement thermal_expansivity() in derived class!")
        return None

    @property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral [J/K/mol]
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        C_v : array of floats
        """
        raise NotImplementedError("need to implement heat_capacity_v() in derived class!")
        return None

    @property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral [J/K/mol]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        C_p : array of floats
        """
        raise NotImplementedError("need to implement heat_capacity_p() in derived class!")
        return None

    @property
    def v_s(self):
        """
        Returns shear wave speed of the mineral [m/s]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        v_s : array of floats
        """
        raise NotImplementedError("need to implement v_s() in derived class!")
        return None

    @property
    def v_p(self):
        """
        Returns P wave speed of the mineral [m/s]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        v_p : array of floats
        """
        raise NotImplementedError("need to implement v_p() in derived class!")
        return None

    @property
    def v_phi(self):
        """
        Returns bulk sound speed of the mineral [m/s]
            
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        v_phi: array of floats
        """
        raise NotImplementedError("need to implement v_phi() in derived class!")
        return None

    @property
    def molar_gibbs(self):
        """
        Returns Gibbs free energy of the mineral [J]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        Gibbs free energies : array of floats
        """
        raise NotImplementedError("need to implement molar_gibbs() in derived class!")
        return None
 
    @property
    def molar_helmholtz(self):
        """
        Returns Helmholtz free energy of the mineral [J]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        Helmholtz free energies : array of floats
        """
        raise NotImplementedError("need to implement molar_helmholtz() in derived class!")
        return None

    @property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral [J]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        enthalpies : array of floats
        """
        raise NotImplementedError("need to implement molar_enthalpy() in derived class!")
        return None

    @property
    def molar_entropy(self):
        """
        Returns entropy of the mineral [J]
        
        Notes
        -----
        Needs to be implemented in derived classes.
        
        Returns
        -------
        entropies : array of floats
        """
        raise NotImplementedError("need to implement molar_entropy() in derived class!")
        return None

    def copy(self):
        return copy.deepcopy(self)




