from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np


def material_property(func):
    """
    Decorator @material_property to be used for cached properties of materials.

    To be used on function in Material or derived classes that should be exposed
    as read-only properties that are cached. The function Material.reset() will
    reset the cached values.

    Internally, the values are stored in a dictionary member called _cached, which
    is emptied by .reset().
    """
    class mat_obj():

        def __init__(self, func):
            self.func = func
            self.varname = self.func.__name__

        def get(self, obj):
            if not hasattr(obj, "_cached"):
                raise Exception("The material_property decorator could not find class member _cached. "
                                "Did you forget to call Material.__init__(self) in __init___?")
            cache_array = getattr(obj, "_cached")
            if self.varname not in cache_array:
                cache_array[self.varname] = self.func(obj)
            return cache_array[self.varname]

    return property(mat_obj(func).get, doc=func.__doc__)


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
        if not hasattr(self, "name"):
            # if a derived class decides to set .name before calling this
            # constructor (I am looking at you, SLB_2011.py!), do not
            # overwrite the name here.
            self._name = self.__class__.__name__
        self._cached = {}

    @property
    def name(self):
        """ Human-readable name of this material.

        By default this will return the name of the class, but it can be set
        to an arbitrary string. Overriden in Mineral.
        """
        return self._name

    @name.setter
    def name(self, value):
        self._name = value

    def set_method(self, method):
        """
        Set the averaging method. See :doc:`averaging` for details.

        Notes
        -----
        Needs to be implemented in derived classes.
        """
        raise NotImplementedError(
            "need to implement set_method() in derived class!")

    def to_string(self):
        """
        Returns a human-readable name of this material. The default implementation will return the name of the class,
        which is a reasonable default.

        Returns
        -------
        name : string
            Name of this material.
        """
        return "'" + self.name + "'"

    def debug_print(self, indent=""):
        """
        Print a human-readable representation of this Material.
        """
        raise NotImplementedError(
            "Derived classes need to implement debug_print(). This is '" + self.__class__.__name__ + "'")

    def print_minerals_of_current_state(self):
        """
        Print a human-readable representation of this Material at the current P, T as a list of minerals.
        This requires set_state() has been called before.
        """
        (minerals, fractions) = self.unroll()
        if len(minerals) == 1:
            print(minerals[0].to_string())
        else:
            print("Material %s:" % self.to_string())
            for (mineral, fraction) in zip(minerals, fractions):
                print("  %g of phase %s" % (fraction, mineral.to_string()))

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
        if not hasattr(self, "_pressure"):
            raise Exception("Material.set_state() could not find class member _pressure. "
                            "Did you forget to call Material.__init__(self) in __init___?")
        self.reset()

        self._pressure = pressure
        self._temperature = temperature

    def reset(self):
        """
        Resets all cached material properties.

        It is typically not required for the user to call this function.
        """
        self._cached = {}

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
        raise NotImplementedError(
            "need to implement unroll() in derived class!")

    def evaluate(self, vars_list, pressures, temperatures):
        """
        Returns an array of material properties requested through a list of strings at given pressure and temperature
        conditions. At the end it resets the set_state to the original values.
        The user needs to call set_method() before.

        Parameters
        ----------
        vars_list : list of strings
            Variables to be returned for given conditions
        pressures : ndlist or ndarray of float
            n-dimensional array of pressures in [Pa].
        temperatures : ndlist or ndarray of float
            n-dimensional array of temperatures in [K].

        Returns
        -------
        output : array of array of float
            Array returning all variables at given pressure/temperature values. output[i][j] is property vars_list[j]
            and temperatures[i] and pressures[i].

        """
        old_pressure = self.pressure
        old_temperature = self.temperature
        pressures = np.array(pressures)
        temperatures = np.array(temperatures)

        assert(pressures.shape == temperatures.shape)
        
        output = np.empty((len(vars_list),) + pressures.shape)
        for i, p in np.ndenumerate(pressures):
            self.set_state(p, temperatures[i])
            for j in range(len(vars_list)):
                output[(j,) + i] = getattr(self, vars_list[j])
        if old_pressure is None or old_temperature is None:
            # do not set_state if old values were None. Just reset to None
            # manually
            self._pressure = self._temperature = None
            self.reset()
        else:
            self.set_state(old_pressure, old_temperature)

        return output

    @property
    def pressure(self):
        """
        Returns current pressure that was set with :func:`~burnman.material.Material.set_state`.


        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.P`.

        Returns
        -------
        pressure : float
            Pressure in [Pa].
        """
        return self._pressure

    @property
    def temperature(self):
        """
        Returns current temperature that was set with :func:`~burnman.material.Material.set_state`.

        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.T`.

        Returns
        -------
        temperature : float
            Temperature in [K].
        """
        return self._temperature

    @material_property
    def internal_energy(self):
        """
        Returns the internal energy of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.energy`.

        Returns
        -------
        internal_energy : float
            The internal energy in [J].
        """
        raise NotImplementedError(
            "need to implement internal_energy() in derived class!")

    @material_property
    def molar_gibbs(self):
        """
        Returns the Gibbs free energy of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gibbs`.

        Returns
        -------
        molar_gibbs : float
            Gibbs free energy in [J].
        """
        raise NotImplementedError(
            "need to implement molar_gibbs() in derived class!")

    @material_property
    def molar_helmholtz(self):
        """
        Returns the Helmholtz free energy of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.helmholtz`.

        Returns
        -------
        molar_helmholtz : float
            Helmholtz free energy in [J].
        """
        raise NotImplementedError(
            "need to implement molar_helmholtz() in derived class!")

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.

        Returns
        -------
        molar_mass : float
            Molar mass in [kg/mol].
        """
        raise NotImplementedError(
            "need to implement molar_mass() in derived class!")

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.V`.

        Returns
        -------
        molar_volume : float
            Molar volume in [m^3/mol].
        """
        raise NotImplementedError(
            "need to implement molar_volume() in derived class!")

    @material_property
    def density(self):
        """
        Returns the density of this material.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.rho`.

        Returns
        -------
        density : float
            The density of this material in [kg/m^3].
        """
        raise NotImplementedError(
            "need to implement density() in derived class!")

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.S`.

        Returns
        -------
        entropy : float
            Entropy in [J].
        """
        raise NotImplementedError(
            "need to implement molar_entropy() in derived class!")

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.H`.

        Returns
        -------
        enthalpy : float
            Enthalpy in [J].
        """
        raise NotImplementedError(
            "need to implement molar_enthalpy() in derived class!")

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the material.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_T`.

        Returns
        -------
        isothermal_bulk_modulus : float
            Bulk modulus in [Pa].
        """
        raise NotImplementedError(
            "need to implement isothermal_bulk_moduls() in derived class!")

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns the adiabatic bulk modulus of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_S`.

        Returns
        -------
        adiabatic_bulk_modulus : float
            Adiabatic bulk modulus in [Pa].
        """
        raise NotImplementedError(
            "need to implement adiabatic_bulk_modulus() in derived class!")

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the mineral (or inverse isothermal bulk modulus).

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_T`.

        Returns
        -------
        (K_T)^-1 : float
            Compressibility in [1/Pa].
        """
        raise NotImplementedError(
            "need to implement compressibility() in derived class!")

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the mineral (or inverse adiabatic bulk modulus).


        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_S`.

        Returns
        -------
        adiabatic_compressibility : float
            adiabatic compressibility in [1/Pa].
        """
        raise NotImplementedError(
            "need to implement compressibility() in derived class!")

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_G`.

        Returns
        -------
        shear_modulus : float
            Shear modulus in [Pa].
        """
        raise NotImplementedError(
            "need to implement shear_modulus() in derived class!")

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_p`.

        Returns
        -------
        p_wave_velocity : float
            P wave speed in [m/s].
        """
        raise NotImplementedError(
            "need to implement p_wave_velocity() in derived class!")

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_phi`.

        Returns
        -------
        bulk sound velocity: float
            Sound velocity in [m/s].
        """
        raise NotImplementedError(
            "need to implement bulk_sound_velocity() in derived class!")

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_s`.

        Returns
        -------
        shear_wave_velocity : float
            Wave speed in [m/s].
        """
        raise NotImplementedError(
            "need to implement shear_wave_velocity() in derived class!")

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameter of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gr`.

        Returns
        -------
        gr : float
            Grueneisen parameters [unitless].
        """
        raise NotImplementedError(
            "need to implement grueneisen_parameter() in derived class!")

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.alpha`.

        Returns
        -------
        alpha : float
            Thermal expansivity in [1/K].
        """
        raise NotImplementedError(
            "need to implement thermal_expansivity() in derived class!")

    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_v`.

        Returns
        -------
        heat_capacity_v : float
            Heat capacity in [J/K/mol].
        """
        raise NotImplementedError(
            "need to implement heat_capacity_v() in derived class!")

    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the mineral.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_p`.

        Returns
        -------
        heat_capacity_p : float
            Heat capacity in [J/K/mol].
        """
        raise NotImplementedError(
            "need to implement heat_capacity_p() in derived class!")

    #
    # Aliased properties
    @property
    def P(self):
        """Alias for :func:`~burnman.material.Material.pressure`"""
        return self.pressure

    @property
    def T(self):
        """Alias for :func:`~burnman.material.Material.temperature`"""
        return self.temperature

    @property
    def energy(self):
        """Alias for :func:`~burnman.material.Material.internal_energy`"""
        return self.internal_energy

    @property
    def helmholtz(self):
        """Alias for :func:`~burnman.material.Material.molar_helmholtz`"""
        return self.molar_helmholtz

    @property
    def gibbs(self):
        """Alias for :func:`~burnman.material.Material.molar_gibbs`"""
        return self.molar_gibbs

    @property
    def V(self):
        """Alias for :func:`~burnman.material.Material.molar_volume`"""
        return self.molar_volume

    @property
    def rho(self):
        """Alias for :func:`~burnman.material.Material.density`"""
        return self.density

    @property
    def S(self):
        """Alias for :func:`~burnman.material.Material.molar_entropy`"""
        return self.molar_entropy

    @property
    def H(self):
        """Alias for :func:`~burnman.material.Material.molar_enthalpy`"""
        return self.molar_enthalpy

    @property
    def K_T(self):
        """Alias for :func:`~burnman.material.Material.isothermal_bulk_modulus`"""
        return self.isothermal_bulk_modulus

    @property
    def K_S(self):
        """Alias for :func:`~burnman.material.Material.adiabatic_bulk_modulus`"""
        return self.adiabatic_bulk_modulus

    @property
    def beta_T(self):
        """Alias for :func:`~burnman.material.Material.isothermal_compressibility`"""
        return self.isothermal_compressibility

    @property
    def beta_S(self):
        """Alias for :func:`~burnman.material.Material.adiabatic_compressibility`"""
        return self.adiabatic_compressibility

    @property
    def G(self):
        """Alias for :func:`~burnman.material.Material.shear_modulus`"""
        return self.shear_modulus

    @property
    def v_p(self):
        """Alias for :func:`~burnman.material.Material.p_wave_velocity`"""
        return self.p_wave_velocity

    @property
    def v_phi(self):
        """Alias for :func:`~burnman.material.Material.bulk_sound_velocity`"""
        return self.bulk_sound_velocity

    @property
    def v_s(self):
        """Alias for :func:`~burnman.material.Material.shear_wave_velocity`"""
        return self.shear_wave_velocity

    @property
    def gr(self):
        """Alias for :func:`~burnman.material.Material.grueneisen_parameter`"""
        return self.grueneisen_parameter

    @property
    def alpha(self):
        """Alias for :func:`~burnman.material.Material.thermal_expansivity`"""
        return self.thermal_expansivity

    @property
    def C_v(self):
        """Alias for :func:`~burnman.material.Material.heat_capacity_v`"""
        return self.heat_capacity_v

    @property
    def C_p(self):
        """Alias for :func:`~burnman.material.Material.heat_capacity_p`"""
        return self.heat_capacity_p
