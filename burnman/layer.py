from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from burnman import averaging_schemes
from burnman import constants
from burnman import geotherm
import warnings


from .material import Material, material_property
from .mineral import Mineral
from .composite import Composite
from .seismic import Seismic1DModel

import matplotlib.pyplot as plt


class Layer(object):
    """
    A planetary layer class
    """

    def __init__(self, name=None, radius_planet=None,
                 min_depth=None, max_depth=None, n_slices=None, verbose=False):

        self.name = name
        self.min_depth = min_depth
        self.max_depth = max_depth
        self.thickness = self.max_depth - self.min_depth
        self.n_slices = n_slices
        self.depths = np.linspace(min_depth, max_depth, n_slices)
        self.radius_planet = radius_planet
        self.radii = self.radius_planet - self.depths
        self.outer_radius = max(self.radii)
        self.inner_radius = min(self.radii)
        self.verbose = verbose
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None
    
    def reset(self):
        """
        Resets all cached material properties.

        It is typically not required for the user to call this function.
        """
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None

    def set_composition(self, composition):
        """
        Set the composition of a Layer with a Material
        """
        assert(isinstance(composition, Material))
        self.composition = composition
        self.reset()

    def set_temperature_mode(self, temperature_mode='adiabat',
                             temperatures=None, temperature_top=None):
        """
        Sets temperature of the layer by user-defined values or as an (modified) adiabat.
        temperature_mode is 'user_defined','adiabatic', or 'modified_adiabat'

        Parameters
        ----------
        temperature_mode : string
        This can be set to 'user_defined','adiabatic', or 'modified_adiabat'
        temperatures : array of float
        The desired fixed temperatures in [K]. Should have same length as defined depths in layer.
        temperature_top : float
        Temperature at the top for an adiabat
        
        Note
        ---------
        'user_defined' = fixes the temperature with the profile input by the users
        'adiabat' = self-consistently computes the adiabat when setting the state of the layer
        'modified_adiabat' = adds the user input array to the adiabat, 
            e.g. allows to implement boundary layers
        """
        self.reset()
        assert(temperature_mode == 'user_defined' or temperature_mode ==
               'adiabat' or temperature_mode == 'modified_adiabat')

        self.temperature_mode = temperature_mode

        if temperature_mode == 'user_defined' or temperature_mode == 'modified_adiabat':
            assert(len(temperatures) == len(self.depths))
            self.usertemperatures = temperatures
        else:
            self.usertemperatures = np.zeros_like(self.depths)

        if temperature_mode == 'adiabat' or temperature_mode == 'modified_adiabat':
            self.temperature_top = temperature_top
        else:
            self.temperature_top = None

    def set_state( self, pressure_mode='selfconsistent', pressures=None, pressure_top=None,
            gravity_bottom=None, n_max_iterations=50):
        """
        Sets the pressure and temperature of the layer by user-defined values are in a self-consistent fashion.
        pressure_mode is 'user_defined' or 'selfconsistent'
        
        Parameters
        ----------
        pressure_mode : string
        This can be set to 'user_defined' or 'selfconsistent'
        pressures : array of floats
        Pressures (Pa) to set layer to ('user_defined'). This should be the same length as defined depths array for the layer
        pressure_top : float
        Pressure (Pa) at the top of the layer. 
        gravity_bottom : float
        gravity (m/s^2) at the bottom the layer
        n_max_iterations : int
        Maximum number of iterations to reach self-consistent pressures (default = 50)
        """
        self.reset()
        assert(pressure_mode == 'user_defined' or pressure_mode == 'selfconsistent')
        assert(self.temperature_mode is not None)

        self.gravity_bottom = gravity_bottom

        if pressure_mode == 'user_defined':
            assert(len(pressures) == len(self.depths))
            self._pressures = pressures
            warnings.warn(
                "By setting the pressures in Layer it is unlikely to be self-consistent")
            self._temperatures = self._evaluate_temperature(
                self._pressures, self.temperature_top)

        if pressure_mode == 'selfconsistent':
            self.pressure_top = pressure_top
            ref_press = np.zeros_like(pressures)
            new_press = self.pressure_top + \
                (self.depths - min(self.depths)) * \
                2.e5  # initial pressure curve guess
            temperatures = self._evaluate_temperature(
                new_press, self.temperature_top)
            # Make it self-consistent!!!
            i = 0

            while i < n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(
                    new_press, temperatures, gravity_bottom, pressure_top)
                temperatures = self._evaluate_temperature(
                    new_press, self.temperature_top)
                rel_err = abs(
                    (max(ref_press) - max(new_press)) / max(new_press))
                if self.verbose:
                    print(
                        "Iteration %i  maximum core pressure error between iterations: %e" %
                        (i, rel_err))

                if rel_err < 1e-5:
                    break

            self._pressures = new_press
            self._temperatures = temperatures

        self.sublayers = []
        for l in range(len(self.depths)):
            self.sublayers.append(self.composition.copy())
            self.sublayers[l].set_state(
                self._pressures[l], self._temperatures[l])

    def _evaluate_temperature(self, pressures=None, temperature_top=None):
        """
        Returns the temperatures of the layer for given pressures.
        Used by set_state()
        """
        if self.temperature_mode == 'adiabat' or self.temperature_mode == 'modified_adiabat':
            adiabat = geotherm.adiabatic(
                pressures, temperature_top, self.composition)
        else:
            adiabat = np.zeros_like(self.depths)
        return adiabat + self.usertemperatures

    def _evaluate_eos(self, pressures, temperatures, gravity_bottom, pressure_top):
        """
        Returns updated gravity and pressure 
        set_state() loops over this until consistency is achieved.
        """
        [density] = self.composition.evaluate(
            ['density'], pressures, temperatures)
        grav = self._compute_gravity(density, gravity_bottom)
        press = self._compute_pressure(density, grav, pressure_top)
        return grav, press

    # Functions needed to compute self-consistent depths-pressures
    def _compute_gravity(self, density, gravity_bottom):
        """
        Computes the gravity of a layer
        Used by _evaluate_eos()
        """
        radii = self.radii[::-1]
        density = np.squeeze(density)[::-1]
        # Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(radii, density)
        # Numerically integrate Poisson's equation

        def poisson(p, x): return 4.0 * np.pi * \
            constants.G * rhofunc(x) * x * x
        grav = np.ravel(
            odeint( poisson, gravity_bottom *
                radii[0] * radii[0], radii))
        
        if radii[0] == 0:
            grav[0] = 0
            grav[1:] = grav[1:] / radii[1:] / radii[1:]
        else:
            grav[:] = grav[:] / radii[:] / radii[:]
        return grav[::-1]

    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        Used by _evaluate_eos()
        """
        # convert radii to depths
        depth = self.depths
        # Make a spline fit of density as a function of depth
        rhofunc = UnivariateSpline(depth, density)
        # Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline(depth, gravity)

        # integrate the hydrostatic equation
        pressure = np.ravel(
            odeint((lambda p, x: gfunc(x) * rhofunc(x)), pressure_top, depth))

        return pressure

    @property
    def mass(self):
        """
        Calculates the mass of the layer [kg]
        """
        mass = 0.0
        radii = self.radii[::-1]
        density = self.density[::-1]
        rhofunc = UnivariateSpline(radii, density)
        mass = np.abs(quad(lambda r: 4 * np.pi * rhofunc(r) *
                           r * r, radii[0], radii[-1])[0])
        return mass

    @property
    def moment_of_inertia(self):
        """
        Returns the moment of inertia of the layer [kg m^2]
        """
        moment = 0.0
        radii = self.radii[::-1]
        density = self.density[::-1]
        rhofunc = UnivariateSpline(radii, density)
        moment = np.abs(quad(lambda r: 8.0 / 3.0 * np.pi * rhofunc(r)
                             * r * r * r * r, radii[0], radii[-1])[0])
        return moment

    @property
    def gravity(self):
        """
        Returns gravity of the layer [m s^(-2)]
        """
        return self._compute_gravity(self.density, self.gravity_bottom)

    @property
    def bullen(self):
        """
        Returns the Bullen parameter
        """
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        phi = self.bulk_sound_velocity * self.bulk_sound_velocity
        dkappadP = np.gradient(kappa, edge_order=2) / \
            np.gradient(self.pressure, edge_order=2)
        dphidz = np.gradient(phi,
                             edge_order=2) / np.gradient(self.depths,
                                                         edge_order=2) / self.gravity
        bullen = dkappadP - dphidz
        return bullen

    @property
    def brunt_vasala(self):
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        brunt_vasala = self.density * self.gravity * \
            self.gravity * (self.bullen - 1.) / kappa
        return brunt_vasala

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
        return self._pressures

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
        return self._temperatures

    @material_property
    def internal_energy(self):
        """
        Returns the internal energy of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.energy`.

        Returns
        -------
        internal_energy : array of floats
            The internal energy in [J].
        """
        return np.array(
            [self.sublayers[i].internal_energy for i in range(len(self.sublayers))])

    @material_property
    def molar_gibbs(self):
        """
        Returns the Gibbs free energy of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gibbs`.

        Returns
        -------
        molar_gibbs : array of floats
            Gibbs free energy in [J].
        """
        return np.array(
            [self.sublayers[i].molar_gibbs for i in range(len(self.sublayers))])

    @material_property
    def molar_helmholtz(self):
        """
        Returns the Helmholtz free energy of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.helmholtz`.

        Returns
        -------
        molar_helmholtz : array of floats
            Helmholtz free energy in [J].
        """
        return np.array(
            [self.sublayers[i].molar_helmholtz for i in range(len(self.sublayers))])

    @material_property
    def molar_mass(self):
        """
        Returns molar mass of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.

        Returns
        -------
        molar_mass : array of floats
            Molar mass in [kg/mol].
        """
        return np.array(
            [self.sublayers[i].molar_mass for i in range(len(self.sublayers))])

    @material_property
    def molar_volume(self):
        """
        Returns molar volume of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.V`.

        Returns
        -------
        molar_volume : array of floats
            Molar volume in [m^3/mol].
        """
        return np.array(
            [self.sublayers[i].molar_volume for i in range(len(self.sublayers))])

    @material_property
    def density(self):
        """
        Returns the density of this layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.rho`.

        Returns
        -------
        density : array of floats
            The density of this material in [kg/m^3].
        """
        return np.array(
            [self.sublayers[i].density for i in range(len(self.sublayers))])

    @material_property
    def molar_entropy(self):
        """
        Returns entropy of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.S`.

        Returns
        -------
        entropy : array of floats
            Entropy in [J].
        """
        return np.array(
            [self.sublayers[i].molar_entropy for i in range(len(self.sublayers))])

    @material_property
    def molar_enthalpy(self):
        """
        Returns enthalpy of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.H`.

        Returns
        -------
        enthalpy : array of floats
            Enthalpy in [J].
        """
        return np.array(
            [self.sublayers[i].molar_enthalpy for i in range(len(self.sublayers))])

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk modulus of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_T`.

        Returns
        -------
        isothermal_bulk_modulus : array of floats
            Bulk modulus in [Pa].
        """
        return np.array(
            [self.sublayers[i].isothermal_bulk_modulus for i in range(len(self.sublayers))])

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns the adiabatic bulk modulus of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_S`.

        Returns
        -------
        adiabatic_bulk_modulus : array of floats
            Adiabatic bulk modulus in [Pa].
        """
        return np.array(
            [self.sublayers[i].adiabatic_bulk_modulus for i in range(len(self.sublayers))])

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibility of the layer (or inverse isothermal bulk modulus).

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_T`.

        Returns
        -------
        (K_T)^-1 : array of floats
            Compressibility in [1/Pa].
        """
        return np.array(
            [self.sublayers[i].isothermal_compressibility for i in range(len(self.sublayers))])

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibility of the layer (or inverse adiabatic bulk modulus).


        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_S`.

        Returns
        -------
        adiabatic_compressibility : array of floats
            adiabatic compressibility in [1/Pa].
        """
        return np.array(
            [self.sublayers[i].adiabatic_compressibility for i in range(len(self.sublayers))])

    @material_property
    def shear_modulus(self):
        """
        Returns shear modulus of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_G`.

        Returns
        -------
        shear_modulus : array of floats
            Shear modulus in [Pa].
        """
        return np.array(
            [self.sublayers[i].shear_modulus for i in range(len(self.sublayers))])

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speed of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_p`.

        Returns
        -------
        p_wave_velocity : array of floats
            P wave speed in [m/s].
        """
        return np.array(
            [self.sublayers[i].p_wave_velocity for i in range(len(self.sublayers))])

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speed of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_phi`.

        Returns
        -------
        bulk sound velocity: array of floats
            Sound velocity in [m/s].
        """
        return np.array(
            [self.sublayers[i].bulk_sound_velocity for i in range(len(self.sublayers))])

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speed of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_s`.

        Returns
        -------
        shear_wave_velocity : array of floats
            Wave speed in [m/s].
        """
        return np.array(
            [self.sublayers[i].shear_wave_velocity for i in range(len(self.sublayers))])

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameter of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gr`.

        Returns
        -------
        gr : array of floats
            Grueneisen parameters [unitless].
        """
        return np.array(
            [self.sublayers[i].grueneisen_parameter for i in range(len(self.sublayers))])

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficient of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.alpha`.

        Returns
        -------
        alpha : array of floats
            Thermal expansivity in [1/K].
        """
        return np.array(
            [self.sublayers[i].thermal_expansivity for i in range(len(self.sublayers))])

    @material_property
    def heat_capacity_v(self):
        """
        Returns heat capacity at constant volume of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_v`.

        Returns
        -------
        heat_capacity_v : array of floats
            Heat capacity in [J/K/mol].
        """
        return np.array(
            [self.sublayers[i].heat_capacity_v for i in range(len(self.sublayers))])

    @material_property
    def heat_capacity_p(self):
        """
        Returns heat capacity at constant pressure of the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_p`.

        Returns
        -------
        heat_capacity_p : array of floats
            Heat capacity in [J/K/mol].
        """
        return np.array(
            [self.sublayers[i].heat_capacity_p for i in range(len(self.sublayers))])

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

