from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import warnings
from .material import material_property


class Planet(object):
    """
    A class to build (self-consistent) Planets made out of Layers
    (``burnman.Layer``). By default the planet is set to be self-consistent
    (with zero pressure at the surface and zero gravity at the center),
    but this can be overwritte using the set_pressure_mode().
    Pressure_modes defined in the individual layers will be ignored.
    If temperature modes are already set for each of the layers, when the
    planet is initialized, the planet will be built immediately.
    """

    def __init__(self, name, layers, n_max_iterations=50,
                 max_delta=1.e-5, verbose=False):
        """
        Parameters
        ----------
        name : string
            Name of planet
        layers : list of ``burnman.Layer``
            Layers to build the planet out of
            (layers are sorted within the planet)
        n_max_iterations : int
            Maximum number of iterations to reach self-consistent planet
            (default = 50)
        max_delta : float
            Relative update to the center pressure of the planet between
            iterations to stop iterations (default = 1.e-5)
        """
        # sort layers
        self.layers = sorted(layers, key=lambda x: x.inner_radius)
        # assert layers attach to one another
        if len(self.layers) > 1:
            for i in range(1, len(self.layers)):
                assert(self.layers[i].inner_radius
                       == self.layers[i - 1].outer_radius)

        self.name = name

        self.radii = self.evaluate(['radii'])
        self.n_slices = len(self.radii)
        self.radius_planet = max(self.radii)
        self.volume = 4./3.*np.pi*np.power(self.radius_planet, 3.)

        for layer in self.layers:
            layer.n_start = np.where(self.radii == layer.inner_radius)[0][-1]
            layer.n_end = np.where(self.radii == layer.outer_radius)[0][0] + 1
        self._cached = {}
        self.verbose = verbose
        self.set_pressure_mode(
            n_max_iterations=n_max_iterations, max_delta=max_delta)

    def __iter__(self):
        """
        Planet object will iterate over Layers.
        """
        return list(self.layers).__iter__()

    def __str__(self):
        """
        Prints details of the planet
        """
        writing = '{0} consists of {1} layers:\n'.format(
            self.name, len(self.layers))
        for layer in self:
            writing = writing + layer.__str__()
        return writing

    def reset(self):
        """
        Resets all cached material properties.
        It is typically not required for the user to call this function.
        """
        self._cached = {}

    def get_layer(self, name):
        """
        Returns a layer with a given name

        Parameters
        ----------
        name : string
            Given name of a layer
        """
        for layer in self.layers:
            if layer.name == name:
                return layer
        raise LookupError()

    def get_layer_by_radius(self, radius):
        """
        Returns a layer in which this radius lies

        Parameters
        ----------
        radius : float
            radius at which to evaluate the layer

        """
        for layer in self.layers:
            if layer.outer_radius >= radius:
                return layer
        raise LookupError()

    def evaluate(self, properties, radlist=None):
        """
        Function that is generally used to evaluate properties
        of the different layers and stitch them together.
        If asking for different radii than the internal radlist,
        pressure and temperature values are interpolated and the
        layer material evaluated at those pressures and
        temperatures.

        Parameters
        ----------
        properties : list of strings
            List of properties to evaluate
        radlist : array of floats
            Radii to evaluate properties at. If left empty,
            internal radius lists are used.

        Returns
        -------
        properties_array : numpy array
            1D or 2D array of requested properties
            (1D if only one property was requested)
        """
        if radlist is None:
            values = np.empty([len(properties),
                               np.sum([len(layer.radii)
                                       for layer in self.layers])])
            for i, prop in enumerate(properties):
                if prop == 'depth':
                    values[i] = np.array(
                        [self.radius_planet - r for layer in self.layers
                         for r in layer.radii])
                else:
                    j = 0
                    for layer in self.layers:
                        vals = getattr(layer, prop)
                        values[i, j:j + len(vals)] = vals
                        j += len(vals)
        else:
            values = np.empty([len(properties), len(radlist)])
            l_idx = [i for i, layer in enumerate(self.layers) for r in radlist
                     if r >= layer.inner_radius and r <= layer.outer_radius]

            for j, r in enumerate(radlist):
                values[:, j] = self.layers[l_idx[j]].evaluate(
                    properties, [r], self.radius_planet).T[0]

        if values.shape[0] == 1:
            values = values[0]
        return values

    def set_pressure_mode(self, pressure_mode='self-consistent',
                          pressures=None, pressure_top=0.,
                          gravity_bottom=0., n_max_iterations=50,
                          max_delta=1.e-5):
        """
        Sets the pressure mode of the planet by user-defined values are in a
        self-consistent fashion.
        pressure_mode is 'user-defined' or 'self-consistent'.
        The default for the planet is self-consistent, with zero pressure at
        the surface and zero pressure at the center.

        Parameters
        ----------
        pressure_mode : string
            This can be set to 'user-defined' or 'self-consistent'
        pressures : array of floats
            Pressures (Pa) to set layer to ('user-defined'). This should be the
            same length as defined radius array for the layer
        pressure_top : float
            Pressure (Pa) at the top of the layer.
        gravity_bottom : float
            gravity (m/s^2) at the bottom the layer
        n_max_iterations : int
            Maximum number of iterations to reach self-consistent pressures
            (default = 50)
        """
        self.reset()
        assert(pressure_mode == 'user-defined' or pressure_mode == 'self-consistent')

        self.pressure_mode = pressure_mode
        self.gravity_bottom = gravity_bottom

        if pressure_mode == 'user-defined':
            assert(len(pressures) == len(self.radii))
            self._pressures = pressures
            warnings.warn('User-defined pressures mean that the planet is '
                          'unlikely to be self-consistent')

        if pressure_mode == 'self-consistent':
            self.pressure_top = pressure_top
            self.n_max_iterations = n_max_iterations
            self.max_delta = max_delta

    def make(self):
        """
        This routine needs to be called before evaluating any properties. If pressures and temperatures are self-consistent, they
        are computed across the planet here. Also initializes an array of materials in each Layer to compute properties from.
        """

        self.reset()
        for layer in self.layers:
            assert(layer.temperature_mode is not None)

        if self.pressure_mode == 'user-defined':
            self._temperatures = self._evaluate_temperature(self._pressures)

        if self.pressure_mode == 'self-consistent':

            new_press = self.pressure_top + \
                (-self.radii + max(self.radii)) * \
                1.e3  # initial pressure curve guess
            temperatures = self._evaluate_temperature(new_press)

            # Make it self-consistent!!!
            i = 0
            while i < self.n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(new_press,
                                                         temperatures,
                                                         self.gravity_bottom,
                                                         self.pressure_top)
                temperatures = self._evaluate_temperature(new_press)
                rel_err = abs((max(ref_press) - max(new_press))
                              / max(new_press))
                if self.verbose:
                    print(f'Iteration {i:0d} maximum relative pressure error: '
                          f'{rel_err:.1e}')

                if rel_err < self.max_delta:
                    break

            self.pressures = new_press
            self.temperatures = temperatures
            self._gravity = new_grav

        for layer in self.layers:
            layer.sublayers = []
            layer.pressures = self.pressures[layer.n_start: layer.n_end]
            layer.temperatures = self.temperatures[layer.n_start: layer.n_end]
            layer.gravity_bottom = self._gravity[layer.n_start - 1]
            layer.pressure_mode = 'set-in-planet'
            for i in range(len(layer.radii)):
                layer.sublayers.append(layer.material.copy())
                layer.sublayers[i].set_state(
                    layer.pressures[i], layer.temperatures[i])

    def _evaluate_eos(self, pressures, temperatures,
                      gravity_bottom, pressure_top):
        """
        Used to update the pressure profile in set_state()
        """
        density = self._evaluate_density(pressures, temperatures)
        grav = self._compute_gravity(density, gravity_bottom)
        press = self._compute_pressure(density, grav, pressure_top)
        return grav, press

    def _evaluate_density(self, pressures, temperatures):
        """
        Used to update the density profile in _evaluate_eos()
        """
        density = []
        for layer in self.layers:
            density.append(layer.material.evaluate(['density'],
                                                   pressures[layer.n_start:layer.n_end],
                                                   temperatures[layer.n_start:layer.n_end]))
        return np.squeeze(np.hstack(density))

    def _evaluate_temperature(self, pressures):
        """
        Returns the temperatures of different layers for given pressures.
        Used by set_state()
        """
        temps = []
        temperature_top = None
        for layer in self.layers[::-1]:
            if temperature_top is None or layer.temperature_top is not None:
                temperature_top = layer.temperature_top
            temps.extend(layer._evaluate_temperature(
                (pressures[layer.n_start:layer.n_end]), temperature_top)[::-1])
            temperature_top = temps[-1]
        return np.hstack(np.squeeze(temps))[::-1]

    def _compute_gravity(self, density, gravity_bottom):
        """
        Calculate the gravity of the planet, based on a density profile.
        This integrates Poisson's equation in radius, under the assumption
        that the planet is laterally homogeneous.
        Used to update the gravity profile in _evaluate_eos()
        """

        start_gravity = gravity_bottom
        grav = []
        for layer in self.layers:
            grav.extend(layer._compute_gravity(
                density[layer.n_start: layer.n_end], start_gravity))
            start_gravity = grav[-1]

        return np.array(grav)

    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.
        This integrates the equation for hydrostatic equilibrium P = rho g z.
        Used to update the pressure profile in _evaluate_eos()
        """
        start_pressure = pressure_top
        press = []
        for layer in self.layers[::-1]:
            press.extend(layer._compute_pressure(density[layer.n_start:layer.n_end],
                                                 gravity[layer.n_start: layer.n_end],
                                                 start_pressure)[::-1])
            start_pressure = press[-1]
        return np.array(press)[::-1]

    @property
    def mass(self):
        """
        calculates the mass of the entire planet [kg]
        """
        return np.sum([layer.mass for layer in self.layers])

    @property
    def average_density(self):
        """
        calculates the average density of the entire planet [kg/m^3]
        """
        return self.mass/self.volume

    @property
    def moment_of_inertia(self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        return np.sum([layer.moment_of_inertia for layer in self.layers])

    @property
    def moment_of_inertia_factor(self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment_factor = self.moment_of_inertia / self.mass / \
            self.radius_planet / self.radius_planet
        return moment_factor

    @property
    def depth(self):
        """
        Returns depth of the layer [m]
        """
        return self.evaluate(['depth'])

    @property
    def gravity(self):
        """
        Returns gravity of the layer [m s^(-2)]
        """
        return self.evaluate(['gravity'])

    @property
    def bullen(self):
        """
        Returns the Bullen parameter
        """
        return self.evaluate(['bullen'])

    @property
    def brunt_vasala(self):
        return self.evaluate(['brunt_vasala'])

    @property
    def pressure(self):
        """
        Returns current pressure that was set with
        :func:`~burnman.Material.set_state`.

        Aliased with :func:`~burnman.Material.P`.

        Returns
        -------
        pressure : array of floats
            Pressure in [Pa].
        """
        return self.pressures

    @property
    def temperature(self):
        """
        Returns current temperature that was set with
        :func:`~burnman.Material.set_state`.

        Aliased with :func:`~burnman.Material.T`.

        Returns
        -------
        temperature : array of floats
            Temperature in [K].
        """
        return self.temperatures

    @material_property
    def molar_internal_energy(self):
        """
        Returns the molar internal energy of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.energy`.

        Returns
        -------
        molar_internal_energy : array of floats
            The internal energy in [J/mol].
        """
        return self.evaluate(['molar_internal_energy'])

    @material_property
    def molar_gibbs(self):
        """
        Returns the molar Gibbs free energy of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.gibbs`.

        Returns
        -------
        molar_gibbs : array of floats
            Gibbs free energy in [J/mol].
        """
        return self.evaluate(['molar_gibbs'])

    @material_property
    def molar_helmholtz(self):
        """
        Returns the molar Helmholtz free energy of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.helmholtz`.

        Returns
        -------
        molar_helmholtz : array of floats
            Helmholtz free energy in [J/mol].
        """
        return self.evaluate(['molar_helmholtz'])

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the planet.

        Needs to be implemented in derived classes.

        Returns
        -------
        molar_mass : array of floats
            Molar mass in [kg/mol].
        """
        return self.evaluate(['molar_mass'])

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.V`.

        Returns
        -------
        molar_volume : array of floats
            Molar volume in [m^3/mol].
        """
        return self.evaluate(['molar_volume'])

    @material_property
    def density(self):
        """
        Returns the density of this planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.rho`.

        Returns
        -------
        density : array of floats
            The density of this material in [kg/m^3].
        """
        return self.evaluate(['density'])

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropy of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.S`.

        Returns
        -------
        molar_entropy : array of floats
            Entropy in [J/mol].
        """
        return self.evaluate(['molar_entropy'])

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpy of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.H`.

        Returns
        -------
        molar_enthalpy : array of floats
            Enthalpy in [J/mol].
        """
        return self.evaluate(['molar_enthalpy'])

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.K_T`.

        Returns
        -------
        isothermal_bulk_modulus : array of floats
            Bulk modulus in [Pa].
        """
        return self.evaluate(['isothermal_bulk_modulus'])

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns the adiabatic bulk modulus of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.K_S`.

        Returns
        -------
        adiabatic_bulk_modulus : array of floats
            Adiabatic bulk modulus in [Pa].
        """
        return self.evaluate(['adiabatic_bulk_modulus'])

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the planet
        (or inverse isothermal bulk modulus).

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_T`.

        Returns
        -------
        (K_T)^-1 : array of floats
            Compressibility in [1/Pa].
        """
        return self.evaluate(['istothermal_compressibility'])

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the planet
        (or inverse adiabatic bulk modulus).

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_S`.

        Returns
        -------
        adiabatic_compressibility : array of floats
            adiabatic compressibility in [1/Pa].
        """
        return self.evaluate(['adiabatic_compressibility'])

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.beta_G`.

        Returns
        -------
        shear_modulus : array of floats
            Shear modulus in [Pa].
        """
        return self.evaluate(['shear_modulus'])

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_p`.

        Returns
        -------
        p_wave_velocity : array of floats
            P wave speed in [m/s].
        """
        return self.evaluate(['p_wave_velocity'])

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_phi`.

        Returns
        -------
        bulk sound velocity: array of floats
            Sound velocity in [m/s].
        """
        return self.evaluate(['bulk_sound_velocity'])

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.v_s`.

        Returns
        -------
        shear_wave_velocity : array of floats
            Wave speed in [m/s].
        """
        return self.evaluate(['shear_wave_velocity'])

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameter of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.gr`.

        Returns
        -------
        gr : array of floats
            Grueneisen parameters [unitless].
        """
        return self.evaluate(['grueneisen_parameter'])

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.alpha`.

        Returns
        -------
        alpha : array of floats
            Thermal expansivity in [1/K].
        """
        return self.evaluate(['thermal_expansivity'])

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volume of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.C_v`.

        Returns
        -------
        molar_heat_capacity_v : array of floats
            Heat capacity in [J/K/mol].
        """
        return self.evaluate(['molar_heat_capacity_v'])

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar heat capacity at constant pressure of the planet.

        Needs to be implemented in derived classes.
        Aliased with :func:`~burnman.Material.C_p`.

        Returns
        -------
        molar_heat_capacity_p : array of floats
            Heat capacity in [J/K/mol].
        """
        return self.evaluate(['molar_heat_capacity_p'])

#
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
    def energy(self):
        """Alias for :func:`~burnman.Material.molar_internal_energy`"""
        return self.molar_internal_energy

    @property
    def helmholtz(self):
        """Alias for :func:`~burnman.Material.molar_helmholtz`"""
        return self.molar_helmholtz

    @property
    def gibbs(self):
        """Alias for :func:`~burnman.Material.molar_gibbs`"""
        return self.molar_gibbs

    @property
    def V(self):
        """Alias for :func:`~burnman.Material.molar_volume`"""
        return self.molar_volume

    @property
    def rho(self):
        """Alias for :func:`~burnman.Material.density`"""
        return self.density

    @property
    def S(self):
        """Alias for :func:`~burnman.Material.molar_entropy`"""
        return self.molar_entropy

    @property
    def H(self):
        """Alias for :func:`~burnman.Material.molar_enthalpy`"""
        return self.molar_enthalpy

    @property
    def K_T(self):
        """Alias for :func:`~burnman.Material.isothermal_bulk_modulus`"""
        return self.isothermal_bulk_modulus

    @property
    def K_S(self):
        """Alias for :func:`~burnman.Material.adiabatic_bulk_modulus`"""
        return self.adiabatic_bulk_modulus

    @property
    def beta_T(self):
        """Alias for :func:`~burnman.Material.isothermal_compressibility`"""
        return self.isothermal_compressibility

    @property
    def beta_S(self):
        """Alias for :func:`~burnman.Material.adiabatic_compressibility`"""
        return self.adiabatic_compressibility

    @property
    def G(self):
        """Alias for :func:`~burnman.Material.shear_modulus`"""
        return self.shear_modulus

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
        """Alias for :func:`~burnman.Material.molar_heat_capacity_v`"""
        return self.molar_heat_capacity_v

    @property
    def C_p(self):
        """Alias for :func:`~burnman.Material.molar_heat_capacity_p`"""
        return self.molar_heat_capacity_p
