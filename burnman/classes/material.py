# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from copy import deepcopy
from ..utils.math import bracket
from scipy.optimize import brentq
from functools import cached_property

cached_property = cached_property  # for easier access


def material_property(func):
    """
    Decorator @material_property to be used for cached properties of materials.

    To be used on function in Material or derived classes that should be exposed
    as read-only properties that are cached. The function Material.reset() will
    reset the cached values.

    Internally, the values are stored in a dictionary member called _cached, which
    is emptied by .reset().
    """

    class mat_obj:
        def __init__(self, func):
            self.func = func
            self.varname = self.func.__name__

        def get(self, obj):
            if not hasattr(obj, "_cached"):
                raise Exception(
                    "The material_property decorator could not find "
                    "class member _cached. "
                    "Did you forget to call Material.__init__(self) in __init___?"
                )
            cache_array = getattr(obj, "_cached")
            if self.varname not in cache_array:
                cache_array[self.varname] = self.func(obj)
            return cache_array[self.varname]

    return property(mat_obj(func).get, doc=func.__doc__)


class Material(object):
    """
    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type :class:`~burnman.Mineral` and their molar
    fractions. This class is available as ``burnman.Material``.

    The user needs to call set_method() (once in the beginning) and set_state()
    before querying the material with unroll() or density().
    """

    def __init__(self):
        self._pressure = None
        self._temperature = None
        self._number_of_moles = 1.0
        if not hasattr(self, "name"):
            # if a derived class decides to set .name before calling this
            # constructor (I am looking at you, SLB_2011.py!), do not
            # overwrite the name here.
            self._name = self.__class__.__name__
        self._cached = {}

    @property
    def name(self):
        """Human-readable name of this material.

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

        .. note:: Needs to be implemented in derived classes.
        """
        raise NotImplementedError("need to implement set_method() in derived class!")

    def to_string(self):
        """
        Returns a human-readable name of this material.
        The default implementation will return the name of the class,
        which is a reasonable default.

        :returns: A human-readable name of the material.
        :rtype: str
        """
        return "'" + self.name + "'"

    def debug_print(self, indent=""):
        """
        Print a human-readable representation of this Material.
        """
        raise NotImplementedError(
            "Derived classes need to implement debug_print(). This is '"
            + self.__class__.__name__
            + "'"
        )

    def print_minerals_of_current_state(self):
        """
        Print a human-readable representation of this Material at the current
        P, T as a list of minerals.
        This requires set_state() has been called before.
        """
        (minerals, fractions) = self.unroll()
        if len(minerals) == 1:
            print(minerals[0].to_string())
        else:
            print("Material %s:" % self.to_string())
            for mineral, fraction in zip(minerals, fractions):
                print("  %g of phase %s" % (fraction, mineral.to_string()))

    def set_state(self, pressure, temperature):
        """
        Set the material to the given pressure and temperature.

        :param pressure: The desired pressure in [Pa].
        :type pressure: float

        :param temperature: The desired temperature in [K].
        :type temperature: float
        """
        if not hasattr(self, "_pressure"):
            raise Exception(
                "Material.set_state() could not find class member _pressure. "
                "Did you forget to call Material.__init__(self) in __init___?"
            )
        self.reset()

        self._pressure = pressure
        self._temperature = temperature

    def set_state_with_volume(
        self, volume, temperature, pressure_guesses=[0.0e9, 10.0e9]
    ):
        """
        This function acts similarly to set_state, but takes volume and
        temperature as input to find the pressure.

        In order to ensure self-consistency, this function does not
        use any pressure functions from the material classes,
        but instead finds the pressure using the brentq root-finding
        method. To provide more context, even if a mineral is being
        evaluated with a P(V, T) equation of state, there might be a
        property modifier G_mod(P, T) added on top - which might then
        introduce a pressure dependent V_mod(P, T).
        We must therefore find the pressure P(V, T) iteratively:
        V = V_base(P, T) + V_mod(P, T), even if there is a function
        for P_base(V_base, T).

        This function is overloaded by the Mineral class.

        :param volume: The desired molar volume of the mineral [m^3].
        :type volume: float

        :param temperature: The desired temperature of the mineral [K].
        :type temperature: float

        :param pressure_guesses: A list of floats denoting the initial
            low and high guesses for bracketing of the pressure [Pa].
            These guesses should preferably bound the correct pressure,
            but do not need to do so. More importantly,
            they should not lie outside the valid region of
            the equation of state. Defaults to [0.e9, 10.e9].
        :type pressure_guesses: list
        """

        def _delta_volume(pressure, volume, temperature):
            # Try to set the state with this pressure,
            # but if the pressure is too low most equations of state
            # fail. In this case, treat the molar_volume as infinite
            # and brentq will try a larger pressure.
            try:
                self.set_state(pressure, temperature)
                return volume - self.molar_volume
            except Exception:
                return -np.inf

        # we need to have a sign change in [a,b] to find a zero.
        args = (volume, temperature)
        try:
            sol = bracket(_delta_volume, pressure_guesses[0], pressure_guesses[1], args)
        except ValueError:
            try:  # Try again with 0 Pa lower bound on the pressure
                sol = bracket(_delta_volume, 0.0e9, pressure_guesses[1], args)
            except ValueError:
                raise Exception(
                    "Cannot find a pressure, perhaps the volume or starting pressures "
                    "are outside the range of validity for the equation of state?"
                )
        pressure = brentq(_delta_volume, sol[0], sol[1], args=args)
        self.set_state(pressure, temperature)

    @property
    def number_of_moles(self):
        """
        Returns the number of moles in this material.

        :returns: Number of moles.
        :rtype: float
        """
        return self._number_of_moles

    @number_of_moles.setter
    def number_of_moles(self, value):
        """
        Sets the number of moles in this material.

        :param value: Number of moles.
        :type value: float
        """
        if value < 0:
            raise ValueError("Moles cannot be negative")
        self._number_of_moles = value

    def reset(self):
        """
        Resets all cached material properties.

        It is typically not required for the user to call this function.
        """
        self._cached = {}

    def copy(self):
        return deepcopy(self)

    def unroll(self):
        """
        Unroll this material into a list of :class:`burnman.Mineral` and their molar
        fractions. All averaging schemes then operate on this list of minerals.
        Note that the return value of this function may depend on the current
        state (temperature, pressure).

        .. note:: Needs to be implemented in derived classes.

        :returns: A list of molar fractions which should sum to 1.0,
            and a list of :class:`burnman.Mineral` objects
            containing the minerals in the material.
        :rtype: tuple
        """
        raise NotImplementedError("need to implement unroll() in derived class!")

    def evaluate(self, vars_list, pressures, temperatures, molar_fractions=None):
        """
        Returns an array of material properties requested through a list of strings
        at given pressure and temperature conditions.
        At the end it resets the set_state to the original values.
        The user needs to call set_method() before.

        :param vars_list: Variables to be returned for given conditions
        :type vars_list: list of strings

        :param pressures: ndlist or ndarray of float of pressures in [Pa].
        :type pressures: :class:`numpy.array`, n-dimensional

        :param temperatures: ndlist or ndarray of float of temperatures in [K].
        :type temperatures: :class:`numpy.array`, n-dimensional

        :returns: List or array returning all variables at given pressure/temperature values.
            output[i][j] is property vars_list[j] for temperatures[i] and pressures[i].
            Attempts to return an array, falls back to a list if the returned properties
            have different shapes.
        :rtype: list or :class:`numpy.array`, n-dimensional
        """
        old_pressure = self.pressure
        old_temperature = self.temperature
        try:
            old_molar_fractions = self.molar_fractions
        except AttributeError:
            old_molar_fractions = None

        pressures = np.array(pressures)
        temperatures = np.array(temperatures)

        assert pressures.shape == temperatures.shape

        first_index = list(np.ndenumerate(temperatures))[0][0]
        if molar_fractions is not None:
            molar_fractions = np.array(molar_fractions)
            self.set_composition(molar_fractions[first_index])
            assert temperatures.shape == molar_fractions.shape[:-1]

        # First, check the output types of all the requested variables:
        self.set_state(pressures.flat[0], temperatures.flat[0])

        output = []
        for j in range(len(vars_list)):
            try:
                var_shape = getattr(self, vars_list[j]).shape
            except AttributeError:
                var_shape = ()
            output.append(np.empty(pressures.shape + var_shape))

        for i, p in np.ndenumerate(pressures):
            if molar_fractions is not None:
                self.set_composition(molar_fractions[i])
            self.set_state(p, temperatures[i])
            for j in range(len(vars_list)):
                output[j][i] = getattr(self, vars_list[j])
        if old_pressure is None or old_temperature is None:
            # do not set_state if old values were None. Just reset to None
            # manually
            self._pressure = self._temperature = None
            self.reset()
        else:
            self.set_state(old_pressure, old_temperature)
            if old_molar_fractions is not None:
                try:
                    self.set_composition(old_molar_fractions)
                except AttributeError:
                    pass

        try:
            output = np.array(output)
        except ValueError:  # if the lists are different shapes
            pass
        return output

    def evaluate_with_volumes(
        self, vars_list, volumes, temperatures, molar_fractions=None
    ):
        """
        Returns an array of material properties requested through a list of strings
        at given volume and temperature conditions.
        At the end it resets the set_state to the original values.
        The user needs to call set_method() before.

        :param vars_list: Variables to be returned for given conditions
        :type vars_list: list of strings

        :param volumes: ndlist or ndarray of float of volumes in [m^3].
        :type volumes: :class:`numpy.array`, n-dimensional

        :param temperatures: ndlist or ndarray of float of temperatures in [K].
        :type temperatures: :class:`numpy.array`, n-dimensional

        :returns: List or array returning all variables at given pressure/temperature values.
            output[i][j] is property vars_list[j] for temperatures[i] and pressures[i].
            Attempts to return an array, falls back to a list if the returned properties
            have different shapes.
        :rtype: list or :class:`numpy.array`, n-dimensional
        """
        old_pressure = self.pressure
        old_temperature = self.temperature
        try:
            old_molar_fractions = self.molar_fractions
        except AttributeError:
            old_molar_fractions = None

        volumes = np.array(volumes)
        temperatures = np.array(temperatures)

        assert volumes.shape == temperatures.shape

        if molar_fractions is not None:
            molar_fractions = np.array(molar_fractions)
            self.set_composition(molar_fractions[0])
            assert temperatures.shape == molar_fractions.shape[:-1]

        # First, check the output types of all the requested variables:
        self.set_state_with_volume(volumes.flat[0], temperatures.flat[0])

        output = []
        for j in range(len(vars_list)):
            try:
                var_shape = getattr(self, vars_list[j]).shape
            except AttributeError:
                var_shape = ()
            output.append(np.empty(volumes.shape + var_shape))

        for i, v in np.ndenumerate(volumes):
            if molar_fractions is not None:
                self.set_composition(molar_fractions[i])
            self.set_state_with_volume(v, temperatures[i])
            for j in range(len(vars_list)):
                output[j][i] = getattr(self, vars_list[j])
        if old_pressure is None or old_temperature is None:
            # do not set_state if old values were None. Just reset to None
            # manually
            self._pressure = self._temperature = None
            self.reset()
        else:
            self.set_state(old_pressure, old_temperature)

            if old_molar_fractions is not None:
                try:
                    self.set_composition(old_molar_fractions)
                except AttributeError:
                    pass

        try:
            output = np.array(output)
        except ValueError:  # if the lists are different shapes
            pass
        return output

    @property
    def pressure(self):
        """
        Returns current pressure that was set with :func:`~burnman.Material.set_state`.

        .. note:: Aliased with :func:`~burnman.Material.P`.

        :returns: Pressure in [Pa].
        :rtype: float
        """
        return self._pressure

    @property
    def temperature(self):
        """
        Returns current temperature that was set with
        :func:`~burnman.Material.set_state`.

        .. note:: Aliased with :func:`~burnman.Material.T`.

        :returns: Temperature in [K].
        :rtype: float
        """
        return self._temperature

    @material_property
    def molar_internal_energy(self):
        """
        Returns the molar internal energy of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.energy`.

        :returns: The internal energy in [J/mol].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement molar_internal_energy() in derived class!"
        )

    @material_property
    def molar_gibbs(self):
        """
        Returns the molar Gibbs free energy of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.gibbs`.

        :returns: Gibbs free energy in [J/mol].
        :rtype: float
        """
        raise NotImplementedError("need to implement molar_gibbs() in derived class!")

    @material_property
    def molar_helmholtz(self):
        """
        Returns the molar Helmholtz free energy of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.helmholtz`.

        :returns: Helmholtz free energy in [J/mol].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement molar_helmholtz() in derived class!"
        )

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the mineral.

        .. note:: Needs to be implemented in derived classes.

        :returns: Molar mass in [kg/mol].
        :rtype: float
        """
        raise NotImplementedError("need to implement molar_mass() in derived class!")

    @material_property
    def mass(self):
        """
        Returns the mass of this material.

        :returns: Mass in [kg].
        :rtype: float
        """
        raise self.molar_mass * self.number_of_moles

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.V`.

        :returns: Molar volume in [m^3/mol].
        :rtype: float
        """
        raise NotImplementedError("need to implement molar_volume() in derived class!")

    @material_property
    def density(self):
        """
        Returns the density of this material.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.rho`.

        :returns: The density of this material in [kg/m^3].
        :rtype: float
        """
        raise NotImplementedError("need to implement density() in derived class!")

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropy of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.S`.

        :returns: Entropy in [J/K/mol].
        :rtype: float
        """
        raise NotImplementedError("need to implement molar_entropy() in derived class!")

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.H`.

        :returns: Enthalpy in [J/mol].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement molar_enthalpy() in derived class!"
        )

    @material_property
    def isothermal_bulk_modulus_reuss(self):
        """
        Returns isothermal bulk modulus of the material.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.K_T`.

        :returns: Isothermal bulk modulus in [Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement isothermal_bulk_modulus_reuss() in derived class!"
        )

    @material_property
    def isentropic_bulk_modulus_reuss(self):
        """
        Returns the adiabatic bulk modulus of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.K_S`.

        :returns: Adiabatic bulk modulus in [Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement isentropic_bulk_modulus_reuss() in derived class!"
        )

    @material_property
    def isothermal_compressibility_reuss(self):
        """
        Returns isothermal compressibility of the mineral
        (or inverse isothermal bulk modulus).

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.beta_T`.

        :returns: Isothermal compressibility in [1/Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement isothermal_compressibility_reuss() in derived class!"
        )

    @material_property
    def isentropic_compressibility_reuss(self):
        """
        Returns adiabatic compressibility of the mineral
        (or inverse adiabatic bulk modulus).


        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.beta_S`.

        :returns: Adiabatic compressibility in [1/Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement isentropic_compressibility_reuss() in derived class!"
        )

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.beta_G`.

        :returns: Shear modulus in [Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement shear_modulus() in derived class! Perhaps you meant to call effective_shear_modulus, or G_eff?"
        )

    @material_property
    def effective_isentropic_bulk_modulus(self):
        """
        Returns the effective isentropic bulk modulus of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.K_eff`.

        :returns: Effective isentropic bulk modulus in [Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement effective_isentropic_bulk_modulus() in derived class!"
        )

    @material_property
    def effective_shear_modulus(self):
        """
        Returns the effective shear modulus of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.G_eff`.

        :returns: Effective shear modulus in [Pa].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement effective_shear_modulus() in derived class!"
        )

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.v_p`.

        :returns: P wave speed in [m/s].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement p_wave_velocity() in derived class!"
        )

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.v_phi`.

        :returns: Bulk sound velocity in [m/s].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement bulk_sound_velocity() in derived class!"
        )

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.v_s`.

        :returns: Shear wave speed in [m/s].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement shear_wave_velocity() in derived class!"
        )

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameter of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.gr`.

        :returns: Grueneisen parameter [unitless].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement grueneisen_parameter() in derived class!"
        )

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.alpha`.

        :returns: Thermal expansivity in [1/K].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement thermal_expansivity() in derived class!"
        )

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.C_v`.

        :returns: Isochoric heat capacity in [J/K/mol].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement molar_heat_capacity_v() in derived class!"
        )

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure of the mineral.

        .. note:: Needs to be implemented in derived classes.
            Aliased with :func:`~burnman.Material.C_p`.

        :returns: Isobaric heat capacity in [J/K/mol].
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement molar_heat_capacity_p() in derived class!"
        )

    @material_property
    def isentropic_thermal_gradient(self):
        """
        :returns: dTdP, the change in temperature with pressure at constant entropy [Pa/K]
        :rtype: float
        """
        raise NotImplementedError(
            "need to implement isentropic_thermal_gradient() in derived class!"
        )

    # Extensive properties
    @property
    def internal_energy(self):
        """Extensive internal energy of the material (J)."""
        return self.molar_internal_energy * self.number_of_moles

    @property
    def helmholtz(self):
        """Extensive Helmholtz free energy of the material (J)."""
        return self.molar_helmholtz * self.number_of_moles

    @property
    def gibbs(self):
        """Extensive Gibbs free energy of the material (J)."""
        return self.molar_gibbs * self.number_of_moles

    @property
    def enthalpy(self):
        """Extensive enthalpy of the material (J)."""
        return self.molar_enthalpy * self.number_of_moles

    @property
    def volume(self):
        """Extensive volume of the material (m^3)."""
        return self.molar_volume * self.number_of_moles

    @property
    def entropy(self):
        """Extensive entropy of the material (J/K)."""
        return self.molar_entropy * self.number_of_moles

    @property
    def heat_capacity_v(self):
        """Extensive heat capacity at constant volume of the material (J/K)."""
        return self.molar_heat_capacity_v * self.number_of_moles

    @property
    def heat_capacity_p(self):
        """Extensive heat capacity at constant pressure of the material (J/K)."""
        return self.molar_heat_capacity_p * self.number_of_moles

    # Aliased properties
    @property
    def P(self):
        """Alias for :func:`~burnman.Material.pressure`"""
        return self.pressure

    @property
    def T(self):
        """Alias for :func:`~burnman.Material.temperature`"""
        return self.temperature

    @property
    def V(self):
        """Alias for :func:`~burnman.Material.volume`"""
        return self.volume

    @property
    def rho(self):
        """Alias for :func:`~burnman.Material.density`"""
        return self.density

    @property
    def S(self):
        """Alias for :func:`~burnman.Material.entropy`"""
        return self.entropy

    @property
    def H(self):
        """Alias for :func:`~burnman.Material.enthalpy`"""
        return self.enthalpy

    @property
    def K_T(self):
        """Alias for :func:`~burnman.Material.isothermal_bulk_modulus_reuss`"""
        return self.isothermal_bulk_modulus_reuss

    @property
    def K_S(self):
        """Alias for :func:`~burnman.Material.isentropic_bulk_modulus_reuss`"""
        return self.isentropic_bulk_modulus_reuss

    @property
    def beta_T(self):
        """Alias for :func:`~burnman.Material.isothermal_compressibility_reuss`"""
        return self.isothermal_compressibility_reuss

    @property
    def beta_S(self):
        """Alias for :func:`~burnman.Material.isentropic_compressibility_reuss`"""
        return self.isentropic_compressibility_reuss

    @property
    def G(self):
        """Alias for :func:`~burnman.Material.shear_modulus`"""
        return self.shear_modulus

    @property
    def K_eff(self):
        """Alias for :func:`~burnman.Material.effective_isentropic_bulk_modulus`"""
        return self.effective_isentropic_bulk_modulus

    @property
    def G_eff(self):
        """Alias for :func:`~burnman.Material.effective_shear_modulus`"""
        return self.effective_shear_modulus

    @property
    def v_p(self):
        """Alias for :func:`~burnman.Material.p_wave_velocity`"""
        return self.p_wave_velocity

    @property
    def v_phi(self):
        """Alias for :func:`~burnman.Material.bulk_sound_velocity`"""
        return self.bulk_sound_velocity

    @property
    def v_s(self):
        """Alias for :func:`~burnman.Material.shear_wave_velocity`"""
        return self.shear_wave_velocity

    @property
    def gr(self):
        """Alias for :func:`~burnman.Material.grueneisen_parameter`"""
        return self.grueneisen_parameter

    @property
    def alpha(self):
        """Alias for :func:`~burnman.Material.thermal_expansivity`"""
        return self.thermal_expansivity

    @property
    def C_v(self):
        """Alias for :func:`~burnman.Material.heat_capacity_v`"""
        return self.heat_capacity_v

    @property
    def C_p(self):
        """Alias for :func:`~burnman.Material.heat_capacity_p`"""
        return self.heat_capacity_p
