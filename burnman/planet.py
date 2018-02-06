from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from burnman import constants
from .material import Material, material_property

class Planet(object):
    """
    A  class to build (self-consistent) Planets made out of Layers (``burnman.Layer``).
    By default the planet is set to be self-consistent (with zero pressure at the surface and zero 
    gravity at the center), but this can be overwritten in the set_pressure_mode().
    Pressure_modes defined in the individual layers will be ignored. If temperature modes are already
    set, the planet will be built upon initialization.
    """

    def __init__(self, name, layers, n_max_iterations=50, max_delta = 1.e-5, verbose=False):
        """
        Parameters
        ----------
        name : string
        Name of planet
        layers : list of ``burnman.Layer``
        Layers to build the planet out of (layers are sorted within the planet)
        n_max_iterations : int
        Maximum number of iterations to reach self-consistent planet (default = 50)
        max_delta : float
        Relative update to the center pressure of the planet between iterations to stop iterations (default = 1.e-5)
        """
        # sort layers
        self.layers = sorted(layers, key=lambda x: x.inner_radius)
        # assert layers attach to one another
        if len(self.layers) > 1:
            for l in range(1, len(self.layers)):
                assert(self.layers[l].inner_radius ==
                       self.layers[l - 1].outer_radius)


        self.name = name

        self.radii = self.evaluate(['radii'])
        self.n_slices = len(self.radii)
        self.radius_planet = max(self.radii)


        for layer in self.layers:
            layer.n_start = np.where(self.radii == layer.inner_radius)[0][-1]
            layer.n_end = np.where(self.radii == layer.outer_radius)[0][0] + 1
        self._cached = {}
        self.verbose = verbose
        self.set_pressure_mode(n_max_iterations=n_max_iterations, max_delta = max_delta)



    def __iter__(self):
        """
        Planet object will iterate over Layers.
        """
        return list(self.layers).__iter__()

    def __str__(self):
        """
        Prints details of the planet
        """
        writing  = 'Planet ' + self.name + ' consists of ' +str(len(self.layers)) +' layers: \n'
        for layer in self:
            writing  =  writing + ' Layer ' + layer.name + ' made out of ' + layer.material.name + ' with ' + layer.temperature_mode + ' temperatures \n'
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
        
        Parameter 
        ---------
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
        
        Parameter
        ---------
        radius : float
        radius to evaluate layer at
        
        """
        for layer in self.layers:
            if layer.outer_radius >= radius:
                return layer
        raise LookupError()
                



    def evaluate(self, properties, radlist=None):
        """
        Function that is generally used to evaluate properties
        of the different layers and stitch them together. 
        If asking for different radii than the internal radiilist,
        values are linearly interpolated.
        
        Parameter
        ---------
        properties : list of strings
        List of properties to evaluate
        radlist : array of floats
        Radii to evaluate properties at. If left empty,
        internal radius lists are used.
        
        Returns
        -------
        2D array of requested properties
        """
        all = None
        if radlist is None:
            for prop in properties:
                oneprop = None
                for layer in self.layers:
                    if oneprop is None:
                        oneprop = getattr(layer, prop)
                    else:
                        oneprop = np.append(oneprop, getattr(layer, prop))
                if all is None:
                    all = np.array(oneprop)
                else:
                    all = np.vstack((all, np.array(oneprop)))
                        
        if radlist is not None:
            for prop in properties:
                oneprop = []
                for r in radlist:
                    layer = self.get_layer_by_radius(r)
                    oneprop.append(layer.evaluate([prop],r))
                if all is None:
                    all = np.array(oneprop)
                else:
                    all = np.vstack((all,np.array(oneprop)))
        return all

    def set_pressure_mode( self, pressure_mode='self-consistent', pressures=None, pressure_top=0.,
              gravity_bottom=0., n_max_iterations=50, max_delta = 1.e-5):
        """
        Sets the pressure mode of the planet by user-defined values are in a self-consistent fashion.
        pressure_mode is 'user-defined' or 'self-consistent'.
        The default for the planet is self-consistent, with zero pressure at the surface and zero pressure at the center.
        
        Parameters
        ----------
        pressure_mode : string
        This can be set to 'user-defined' or 'self-consistent'
        pressures : array of floats
        Pressures (Pa) to set layer to ('user-defined'). This should be the same length as defined radius array for the layer
        pressure_top : float
        Pressure (Pa) at the top of the layer.
        gravity_bottom : float
        gravity (m/s^2) at the bottom the layer
        n_max_iterations : int
        Maximum number of iterations to reach self-consistent pressures (default = 50)
        """
        self.reset()
        assert(pressure_mode == 'user-defined' or pressure_mode == 'self-consistent')
        
        
        self.pressure_mode = pressure_mode
        self.gravity_bottom = gravity_bottom
 
        if pressure_mode == 'user-defined':
            assert(len(pressures) == len(self.radii))
            self._pressures = pressures
            warnings.warn(
                "By setting the pressures in Planet it is unlikely to be self-consistent")

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
                2.e5  # initial pressure curve guess
            temperatures = self._evaluate_temperature( new_press)
            # Make it self-consistent!!!
            i = 0

            while i < self.n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(
                    new_press, temperatures, self.gravity_bottom, self.pressure_top)
                temperatures = self._evaluate_temperature( new_press)
                rel_err = abs(
                    (max(ref_press) - max(new_press)) / max(new_press))
                if self.verbose:
                    print('Iteration {0:0d} maximum relative pressure error: {1:.1e}'.format(i, rel_err))

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
            for l in range(len(layer.radii)):
                layer.sublayers.append(layer.material.copy())
                layer.sublayers[l].set_state(
                    layer.pressures[l], layer.temperatures[l])

    def _evaluate_eos(self, pressures, temperatures, gravity_bottom, pressure_top):
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
            density.append(layer.material.evaluate(
                ['density'], pressures[layer.n_start:layer.n_end], temperatures[layer.n_start:layer.n_end]))
        return np.squeeze(np.hstack(density))

    def _evaluate_temperature(self, pressures):
        """
        Returns the temperatures of different layers for given pressures.
        Used by set_state()
        """
        temps = []
        temperature_top = None
        for layer in self.layers[::-1]:
            if temperature_top is None:
                temperature_top = layer.temperature_top
            temps.extend(layer._evaluate_temperature((pressures[layer.n_start:layer.n_end]), temperature_top)[::-1])
            temperature_top = temps[-1]
        return np.hstack(np.squeeze(temps))[::-1]

    def _compute_gravity(self, density, gravity_bottom):
        """
        Calculate the gravity of the planet, based on a density profile.  This integrates
        Poisson's equation in radius, under the assumption that the planet is laterally
        homogeneous.
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
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        Used to update the pressure profile in _evaluate_eos()
        """
        start_pressure = pressure_top
        press = []
        for layer in self.layers[::-1]:
            press.extend(layer._compute_pressure(
                                                 density[layer.n_start: layer.n_end], gravity[layer.n_start: layer.n_end], start_pressure)[::-1])
            start_pressure = press[-1]
        return np.array(press)[::-1]

    @property
    def mass(self):
        """
        calculates the mass of the entire planet [kg]
        """
        mass = 0.0
        for layer in self.layers:
            mass += layer.mass
        return mass

    @property
    def moment_of_inertia(self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment = 0.0
        for layer in self.layers:
            moment += layer.moment_of_inertia
        return moment

    @property
    def moment_of_inertia_factor(self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment_factor = self.moment_of_inertia / self.mass / \
            self.radius_planet / self.radius_planet
        return moment_factor

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
        Returns current pressure that was set with :func:`~burnman.material.Material.set_state`.


        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.P`.

        Returns
        -------
        pressure : array of floats
            Pressure in [Pa].
        """
        return self.pressures

    @property
    def temperature(self):
        """
        Returns current temperature that was set with :func:`~burnman.material.Material.set_state`.

        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.T`.

        Returns
        -------
        temperature : array of floats
            Temperature in [K].
        """
        return self.temperatures

    @material_property
    def internal_energy(self):
        """
        Returns the internal energy of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.energy`.

        Returns
        -------
        internal_energy : array of floats
            The internal energy in [J].
        """
        return self.evaluate(['internal_energy'])

    @material_property
    def molar_gibbs(self):
        """
        Returns the Gibbs free energy of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gibbs`.

        Returns
        -------
        molar_gibbs : array of floats
            Gibbs free energy in [J].
        """
        return self.evaluate(['molar_gibbs'])

    @material_property
    def molar_helmholtz(self):
        """
        Returns the Helmholtz free energy of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.helmholtz`.

        Returns
        -------
        molar_helmholtz : array of floats
            Helmholtz free energy in [J].
        """
        return self.evaluate(['molar_helmholtz'])

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.V`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.rho`.

        Returns
        -------
        density : array of floats
            The density of this material in [kg/m^3].
        """
        return self.evaluate(['density'])

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.S`.

        Returns
        -------
        entropy : array of floats
            Entropy in [J].
        """
        return self.evaluate(['molar_entropy'])

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.H`.

        Returns
        -------
        enthalpy : array of floats
            Enthalpy in [J].
        """
        return self.evaluate(['molar_enthalpy'])

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_T`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_S`.

        Returns
        -------
        adiabatic_bulk_modulus : array of floats
            Adiabatic bulk modulus in [Pa].
        """
        return self.evaluate(['adiabatic_bulk_modulus'])

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the planet (or inverse isothermal bulk modulus).

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_T`.

        Returns
        -------
        (K_T)^-1 : array of floats
            Compressibility in [1/Pa].
        """
        return self.evaluate(['istothermal_compressibility'])

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the planet (or inverse adiabatic bulk modulus).


        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_S`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_G`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_p`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_phi`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_s`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gr`.

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

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.alpha`.

        Returns
        -------
        alpha : array of floats
            Thermal expansivity in [1/K].
        """
        return self.evaluate(['thermal_expansivity'])

    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_v`.

        Returns
        -------
        heat_capacity_v : array of floats
            Heat capacity in [J/K/mol].
        """
        return self.evaluate(['heat_capacity_v'])

    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the planet.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_p`.

        Returns
        -------
        heat_capacity_p : array of floats
            Heat capacity in [J/K/mol].
        """
        return self.evaluate(['heat_capacity_p'])

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

