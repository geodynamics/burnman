from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline, interp1d
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
    The base class for a planetary layer. The user needs to set the following befor properties can be computed
    1. set_material(), set the material of the layer, e.g. a mineral, solid_solution, or composite
    2. set_temperature_mode(), either predefine, or set to an adiabatic profile
    3. set_pressure_mode(), to set the self-consistent pressure (with user-defined option the pressures can be overwritten).
    To set the self-consistent pressure the pressure at the top and the gravity at the bottom of the layer need to be set.
    4. make(), computes the self-consistent part of the layer and starts the settings to compute properties within the layer
    Note that not the entire planet this layer sits in is self-consistent, as the pressure at the top of the layer is a function of the density
    within the layer (through the gravity). Entire planets can be computed self-consistently with the planet class.
    Properties will be returned at the pre-defined radius array, although the evaluate() function
    can take a new defined depthlist and values are interpolated between these (sufficient sampling of the layer
    is needed for this to be accurate)
    """

    def __init__(self, name=None, radii=None, verbose=False):
        self.name = name
        assert np.all(np.diff(radii) > 0)
        self.radii = radii
        self.outer_radius = max(self.radii)
        self.inner_radius = min(self.radii)
        self.thickness = self.outer_radius - self.inner_radius
        self.n_slices = len(self.radii)
        self.verbose = verbose
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None
        self.material = None
        self.pressure_mode = 'self-consistent'
        self.temperature_mode = None
 
 
    def __str__(self):
        """
        Prints details of the layer
        """
        writing = 'The {0} is made of {1} with {2} temperatures and {3} pressures\n'.format(self.name,
                                                                                            self.material.name,
                                                                                            self.temperature_mode,
                                                                                            self.pressure_mode)
        return writing
            
            
    def reset(self):
        """
        Resets all cached material properties.
        It is typically not required for the user to call this function.
        """
        self._cached = {}
        self._pressures = None
        self._temperatures = None
        self.sublayers = None

    def set_material(self, material):
        """
        Set the material of a Layer with a Material
        """
        assert(isinstance(material, Material))
        self.material = material
        self.reset()

    def set_temperature_mode(self, temperature_mode='adiabatic',
                             temperatures=None, temperature_top=None):
        """
        Sets temperature of the layer by user-defined values or as an (perturbed) adiabatic.
        temperature_mode is 'user-defined','adiabatic', or 'perturbed-adiabatic'

        Parameters
        ----------
        temperature_mode : string
        This can be set to 'user-defined', 'adiabatic', or 'perturbed-adiabatic'
        temperatures : array of float
        The desired fixed temperatures in [K]. Should have same length as defined radii in layer.
        temperature_top : float
        Temperature at the top for an adiabatic
        
        Note
        ---------
        'user-defined' = fixes the temperature with the profile input by the users
        'adiabatic' = self-consistently computes the adiabat when setting the state of the layer
        'perturbed-adiabatic' = adds the user input array to the adiabat,
            e.g. allows to implement boundary layers
        """
        self.reset()
        assert(temperature_mode == 'user-defined'  or temperature_mode ==
               'adiabatic' or temperature_mode == 'perturbed-adiabatic')

        self.temperature_mode = temperature_mode

#if temperature_mode=='isothermal':
#            self.usertemperatures =
        if temperature_mode == 'user-defined' or temperature_mode == 'perturbed-adiabatic':
            assert(len(temperatures) == len(self.radii))
            self.usertemperatures = temperatures
        else:
            self.usertemperatures = np.zeros_like(self.radii)

        if temperature_mode == 'adiabatic' or temperature_mode == 'perturbed-adiabatic':
            self.temperature_top = temperature_top
        else:
            self.temperature_top = None

    def set_pressure_mode( self, pressure_mode='self-consistent', pressures=None,gravity_bottom=None,
                          pressure_top=None, n_max_iterations=50, max_delta = 1.e-5):
        """
        Sets the pressure mode of the layer, which can either be 'user-defined', or 'self-consistent'.
        
        Parameters
        ----------
        pressure_mode : string
        This can be set to 'user-defined' or 'self-consistent'
        pressures : array of floats
        Pressures [Pa] to set layer to ('user-defined'). This should be the same length as defined radii array for the layer
        pressure_top : float
        Pressure [Pa] at the top of the layer.
        gravity_bottom : float
        gravity [m/s^2] at the bottom the layer
        n_max_iterations : int
        Maximum number of iterations to reach self-consistent pressures (default = 50)
        max_delta : float
        Relative update to the highest pressure in the layer between iterations to stop iterations (default = 1.e-5)
        
        
        Note
        ---------
        'user-defined' = fixes the pressures with the profile input by the user in the 'pressures' argument. 
        'self-consistent' = user needs to give gravity_bottom [m/s^2] and pressure_top [Pa] to compute self-consistent pressures.
        """
        self.reset()
        assert(pressure_mode == 'user-defined' or pressure_mode == 'self-consistent')
        self.pressure_mode= pressure_mode
        
        assert(gravity_bottom is not None)
        self.gravity_bottom = gravity_bottom
        
        if pressure_mode == 'user-defined':
            assert(pressures is not None)
            assert(len(pressures) == len(self.radii))
            self.pressures = pressures
            warnings.warn("By setting the pressures in Layer it is unlikely to be self-consistent")
        elif pressure_mode == 'self-consistent':
            self.pressure_top = pressure_top
            self.n_max_iterations = n_max_iterations
            self.max_delta = max_delta
        else:
            raise NotImplementedError('pressure mode \"{0}\"not recognised'.format(pressure_mode))
            
            

    def make( self):
        """
        This routine needs to be called before evaluating any properties. If pressures and temperatures are self-consistent, they
        are computed here. Also initializes an array of materials to compute properties from.
        """
        self.reset()
        if not hasattr(self,'material'):
            raise AttributeError(' set_meterial() for layer before make()')
        if not hasattr(self,'temperature_mode'):
            raise AttributeError(' set_temperature_mode() for layer before make()')
        if not hasattr(self,'pressure_mode'):
            raise AttributeError(' set_pressure_mode() for layer before make() ')


        if self.pressure_mode == 'user-defined':
            self.temperatures = self._evaluate_temperature(
                self.pressures, self.temperature_top)
        elif self.pressure_mode == 'self-consistent':
            new_press = self.pressure_top + \
                (-self.radii + max(self.radii)) * \
                1.e3  # initial pressure curve guess
            temperatures = self._evaluate_temperature(
                new_press, self.temperature_top)
            # Make it self-consistent!!!
            i = 0

            while i < self.n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(
                    new_press, temperatures, self.gravity_bottom, self.pressure_top)
                temperatures = self._evaluate_temperature(
                    new_press, self.temperature_top)
                rel_err = abs(
                    (max(ref_press) - max(new_press)) / max(new_press))
                if self.verbose:
                    print('Iteration {0:0d} maximum relative pressure error: {1:.1f}'.format(i, rel_err))

                if rel_err < self.max_delta:
                    break

            self.pressures = new_press
            self.temperatures = temperatures
        else:
            raise NotImplementedError('pressure mode not recognised')
        
        self.sublayers = []
        for l in range(len(self.radii)):
            self.sublayers.append(self.material.copy())
            self.sublayers[l].set_state(
                self.pressures[l], self.temperatures[l])

    def evaluate(self, properties, radlist=None, radius_planet=None):
        """
        Function that is generally used to evaluate properties
        across the layer. If radlist is not defined, valuse are
        returned at the internal radlist.
        If asking for different radii than the internal radlist,
        pressure and temperature values are interpolated and the
        layer material evaluated at those pressures and 
        temperatures.
        
        Parameter
        ---------
        properties : list of strings
        List of properties to evaluate
        radlist : array of floats
        Radii to evaluate properties at. If left empty,
        internal radii list is used.
        planet_radius : float
        Planet outer radius. Used only to calculate depth.
        
        Returns
        -------
        1D or 2D array of requested properties 
        (1D if only one property was requested)
        """
        
        if radlist is None:
            values = np.empty([len(properties), len(self.radii)])
            for i, prop in enumerate(properties):
                if prop == 'depth':
                    values[i] = radius_planet - self.radii
                else:
                    try:
                        values[i] = getattr(self,prop)
                    except:
                        values[i]= np.array([getattr(self.sublayers[i],prop)
                                             for i in range(len(self.sublayers))])   
        else:
            func_p = interp1d(self.radii,self.pressures)
            pressures = func_p(radlist)
            func_t = interp1d(self.radii,self.temperatures)
            temperatures = func_t(radlist)
            values = np.empty([len(properties), len(radlist)])
            for i, prop in enumerate(properties):
                if prop == 'depth':
                    values[i] = radius_planet - radlist
                else:
                    try:
                        values[i] = self.material.evaluate([prop],pressures,temperatures)
                    except:
                        func_prop = interp1d(self.radii, getattr(self,prop))
                        values[i] = func_prop(radlist)
                    
        if values.shape[0] == 1:
            values = values[0]
        return values

    def _evaluate_temperature(self, pressures=None, temperature_top=None):
        """
        Returns the temperatures of the layer for given pressures.
        Used by make()
        """
        if self.temperature_mode == 'adiabatic' or self.temperature_mode == 'perturbed-adiabatic':
            adiabat = geotherm.adiabatic(pressures[::-1], temperature_top, self.material)[::-1]
        else:
            adiabat = np.zeros_like(self.radii)
        return adiabat + self.usertemperatures

    def _evaluate_eos(self, pressures, temperatures, gravity_bottom, pressure_top):
        """
        Returns updated gravity and pressure 
        make() loops over this until consistency is achieved.
        """
        [density] = self.material.evaluate(
            ['density'], pressures, temperatures)
        grav = self._compute_gravity(density, gravity_bottom)
        press = self._compute_pressure(density, grav, pressure_top)
        return grav, press

    # Functions needed to compute self-consistent radii-pressures
    def _compute_gravity(self, density, gravity_bottom):
        """
        Computes the gravity of a layer
        Used by _evaluate_eos()
        """
        # Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(self.radii, density)
        # Numerically integrate Poisson's equation

        def poisson(p, x): return 4.0 * np.pi * \
            constants.G * rhofunc(x) * x * x
        grav = np.ravel(
            odeint( poisson, gravity_bottom *
                self.radii[0] * self.radii[0], self.radii))
        
        if self.radii[0] == 0:
            grav[0] = 0
            grav[1:] = grav[1:] / self.radii[1:] / self.radii[1:]
        else:
            grav[:] = grav[:] / self.radii[:] / self.radii[:]
        return grav

    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        Used by _evaluate_eos()
        """
        # flip radius, density and gravity to increasing pressure
        depthfromtop = -self.radii[::-1] + max(self.radii)
        density = density[::-1]
        gravity = gravity[::-1]
        # Make a spline fit of density as a function of depth
        rhofunc = UnivariateSpline(depthfromtop, density)
        # Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline(depthfromtop, gravity)

        # integrate the hydrostatic equation
        pressure = np.ravel(
            odeint((lambda p, x: gfunc(x) * rhofunc(x)), pressure_top, depthfromtop))

        return pressure[::-1]

    @property
    def mass(self):
        """
        Calculates the mass of the layer [kg]
        """
        mass = 0.0
        radii = self.radii
        density = self.evaluate(['density'])
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
        radii = self.radii
        density = self.evaluate(['density'])
        rhofunc = UnivariateSpline(radii, density)
        moment = np.abs(quad(lambda r: 8.0 / 3.0 * np.pi * rhofunc(r)
                             * r * r * r * r, radii[0], radii[-1])[0])
        return moment

    @property
    def gravity(self):
        """
        Returns gravity profile of the layer [m s^(-2)]
        """
        return self._compute_gravity(self.density, self.gravity_bottom)

    @property
    def bullen(self):
        """
        Returns the Bullen parameter across the layer. 
        The Bullen parameter assess if compression as a function of pressure is 
        like homogeneous, adiabatic compression. 
        Bullen parameter =1  , homogeneous, adiabatic compression
        Bullen parameter > 1 , more compressed with pressure, e.g. across phase transitions
        Bullen parameter < 1, less compressed with pressure, e.g. across a boundary layer
        """
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        phi = self.bulk_sound_velocity * self.bulk_sound_velocity
        try:
            dkappadP = np.gradient(kappa, edge_order=2) / \
                       np.gradient(self.pressures, edge_order=2)
            dphidr = np.gradient(phi,edge_order=2) / np.gradient(self.radii,edge_order=2) / self.gravity
        except:
            dkappadP = np.gradient(kappa) / \
                       np.gradient(self.pressures)
            dphidr = np.gradient(phi) / np.gradient(self.radii) / self.gravity
        bullen = dkappadP + dphidr
        return bullen

    @property
    def brunt_vasala(self):
        """
        Returns the brunt-vasala (or buoyancy) frequency, N, across the layer.
        This frequency assess the stabilty of the layer:
        N < 0, fluid will convect
        N= 0, fluid is neutral
        N > 0, fluid is stabily stratified.
        """
        kappa = self.bulk_sound_velocity * self.bulk_sound_velocity * self.density
        brunt_vasala = self.density * self.gravity * \
            self.gravity * (self.bullen - 1.) / kappa
        return brunt_vasala

    @property
    def pressure(self):
        """
        Returns current pressures across the layer that was set with :func:`~burnman.material.Material.set_state`.


        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.P`.

        Returns
        -------
        pressure : array of floats
            Pressures in [Pa] at the predefined radii.
        """
        return self.pressures

    @property
    def temperature(self):
        """
        Returns current temperature  across the layer that was set with :func:`~burnman.material.Material.set_state`.

        Notes
        -----
        - Aliased with :func:`~burnman.material.Material.T`.

        Returns
        -------
        temperature : array of floats
            Temperatures in [K] at the predefined radii.
        """
        return self.temperatures


    @material_property
    def molar_internal_energy(self):
        """
        Returns the molar internal energies across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.energy`.

        Returns
        -------
        molar_internal_energy : array of floats
            The internal energies in [J/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_internal_energy for i in range(len(self.sublayers))])

    @material_property
    def molar_gibbs(self):
        """
        Returns the molar Gibbs free energies across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gibbs`.

        Returns
        -------
        molar_gibbs : array of floats
            Gibbs free energies in [J/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_gibbs for i in range(len(self.sublayers))])

    @material_property
    def molar_helmholtz(self):
        """
        Returns the molar Helmholtz free energies across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.helmholtz`.

        Returns
        -------
        molar_helmholtz : array of floats
            Helmholtz free energies in [J/mol] at the predefined radii.
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
        Returns molar volumes across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.V`.

        Returns
        -------
        molar_volume : array of floats
            Molar volumes in [m^3/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_volume for i in range(len(self.sublayers))])

    @material_property
    def density(self):
        """
        Returns the densities across this layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.rho`.

        Returns
        -------
        density : array of floats
            The densities of this material in [kg/m^3] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].density for i in range(len(self.sublayers))])

    @material_property
    def molar_entropy(self):
        """
        Returns molar entropies acroos the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.S`.

        Returns
        -------
        molar_entropy : array of floats
            Entropies in [J/K/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_entropy for i in range(len(self.sublayers))])

    @material_property
    def molar_enthalpy(self):
        """
        Returns molar enthalpies across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.H`.

        Returns
        -------
        molar_enthalpy : array of floats
            Enthalpies in [J/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_enthalpy for i in range(len(self.sublayers))])

    @material_property
    def isothermal_bulk_modulus(self):
        """
        Returns isothermal bulk moduli across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_T`.

        Returns
        -------
        isothermal_bulk_modulus : array of floats
            Bulk moduli in [Pa] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].isothermal_bulk_modulus for i in range(len(self.sublayers))])

    @material_property
    def adiabatic_bulk_modulus(self):
        """
        Returns the adiabatic bulk moduli across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.K_S`.

        Returns
        -------
        adiabatic_bulk_modulus : array of floats
            Adiabatic bulk modulus in [Pa] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].adiabatic_bulk_modulus for i in range(len(self.sublayers))])

    @material_property
    def isothermal_compressibility(self):
        """
        Returns isothermal compressibilities across the layer (or inverse isothermal bulk moduli).

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_T`.

        Returns
        -------
        (K_T)^-1 : array of floats
            Compressibilities in [1/Pa] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].isothermal_compressibility for i in range(len(self.sublayers))])

    @material_property
    def adiabatic_compressibility(self):
        """
        Returns adiabatic compressibilities across the layer (or inverse adiabatic bulk moduli).


        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_S`.

        Returns
        -------
        adiabatic_compressibility : array of floats
            adiabatic compressibilities in [1/Pa] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].adiabatic_compressibility for i in range(len(self.sublayers))])

    @material_property
    def shear_modulus(self):
        """
        Returns shear moduli across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.beta_G`.

        Returns
        -------
        shear_modulus : array of floats
            Shear modulie in [Pa]  at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].shear_modulus for i in range(len(self.sublayers))])

    @material_property
    def p_wave_velocity(self):
        """
        Returns P wave speeds across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_p`.

        Returns
        -------
        p_wave_velocity : array of floats
            P wave speeds in [m/s] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].p_wave_velocity for i in range(len(self.sublayers))])

    @material_property
    def bulk_sound_velocity(self):
        """
        Returns bulk sound speeds across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_phi`.

        Returns
        -------
        bulk sound velocity: array of floats
            Sound velocities in [m/s] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].bulk_sound_velocity for i in range(len(self.sublayers))])

    @material_property
    def shear_wave_velocity(self):
        """
        Returns shear wave speeds across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.v_s`.

        Returns
        -------
        shear_wave_velocity : array of floats
            Wave speeds in [m/s] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].shear_wave_velocity for i in range(len(self.sublayers))])

    @material_property
    def grueneisen_parameter(self):
        """
        Returns the grueneisen parameters across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.gr`.

        Returns
        -------
        gr : array of floats
            Grueneisen parameters [unitless] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].grueneisen_parameter for i in range(len(self.sublayers))])

    @material_property
    def thermal_expansivity(self):
        """
        Returns thermal expansion coefficients across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.alpha`.

        Returns
        -------
        alpha : array of floats
            Thermal expansivities in [1/K] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].thermal_expansivity for i in range(len(self.sublayers))])

    @material_property
    def molar_heat_capacity_v(self):
        """
        Returns molar heat capacity at constant volumes across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_v`.

        Returns
        -------
        molar_heat_capacity_v : array of floats
            Heat capacities in [J/K/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_heat_capacity_v for i in range(len(self.sublayers))])

    @material_property
    def molar_heat_capacity_p(self):
        """
        Returns molar_heat capacity at constant pressures across the layer.

        Notes
        -----
        - Needs to be implemented in derived classes.
        - Aliased with :func:`~burnman.material.Material.C_p`.

        Returns
        -------
        molar_heat_capacity_p : array of floats
            Heat capacities in [J/K/mol] at the predefined radii.
        """
        return np.array(
            [self.sublayers[i].molar_heat_capacity_p for i in range(len(self.sublayers))])

#
# Aliased properties
    @property
    def P(self):
        """Alias for :func:`~burnman.layer.Layer.pressure`"""
        return self.pressure
    
    @property
    def T(self):
        """Alias for :func:`~burnman.layer.Layer.temperature`"""
        return self.temperature
    
    @property
    def energy(self):
        """Alias for :func:`~burnman.layer.Layer.molar_internal_energy`"""
        return self.molar_internal_energy
    
    @property
    def helmholtz(self):
        """Alias for :func:`~burnman.layer.Layer.molar_helmholtz`"""
        return self.molar_helmholtz
    
    @property
    def gibbs(self):
        """Alias for :func:`~burnman.layer.Layer.molar_gibbs`"""
        return self.molar_gibbs
    
    @property
    def V(self):
        """Alias for :func:`~burnman.layer.Layer.molar_volume`"""
        return self.molar_volume
    
    @property
    def rho(self):
        """Alias for :func:`~burnman.layer.Layer.density`"""
        return self.density
    
    @property
    def S(self):
        """Alias for :func:`~burnman.layer.Layer.molar_entropy`"""
        return self.molar_entropy
    
    @property
    def H(self):
        """Alias for :func:`~burnman.layer.Layer.molar_enthalpy`"""
        return self.molar_enthalpy
    
    @property
    def K_T(self):
        """Alias for :func:`~burnman.layer.Layer.isothermal_bulk_modulus`"""
        return self.isothermal_bulk_modulus
    
    @property
    def K_S(self):
        """Alias for :func:`~burnman.layer.Layer.adiabatic_bulk_modulus`"""
        return self.adiabatic_bulk_modulus
    
    @property
    def beta_T(self):
        """Alias for :func:`~burnman.layer.Layer.isothermal_compressibility`"""
        return self.isothermal_compressibility
    
    @property
    def beta_S(self):
        """Alias for :func:`~burnman.layer.Layer.adiabatic_compressibility`"""
        return self.adiabatic_compressibility
    
    @property
    def G(self):
        """Alias for :func:`~burnman.layer.Layer.shear_modulus`"""
        return self.shear_modulus
    
    @property
    def v_p(self):
        """Alias for :func:`~burnman.layer.Layer.p_wave_velocity`"""
        return self.p_wave_velocity
    
    @property
    def v_phi(self):
        """Alias for :func:`~burnman.layer.Layer.bulk_sound_velocity`"""
        return self.bulk_sound_velocity
    
    @property
    def v_s(self):
        """Alias for :func:`~burnman.layer.Layer.shear_wave_velocity`"""
        return self.shear_wave_velocity
    
    @property
    def gr(self):
        """Alias for :func:`~burnman.layer.Layer.grueneisen_parameter`"""
        return self.grueneisen_parameter
    
    @property
    def alpha(self):
        """Alias for :func:`~burnman.layer.Layer.thermal_expansivity`"""
        return self.thermal_expansivity
    
    @property
    def C_v(self):
        """Alias for :func:`~burnman.material.Material.molar_heat_capacity_v`"""
        return self.molar_heat_capacity_v
    
    @property
    def C_p(self):
        """Alias for :func:`~burnman.material.Material.molar_heat_capacity_p`"""
        return self.molar_heat_capacity_p

