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


class Layer(Material):
    """
    A planetary layer class
    """

    def __init__(self, name=None, radius_planet=None,
                 min_depth=None, max_depth=None, n_slices=None, verbose=False):
        Material.__init__(self)

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

    def set_composition(self, composition):
        """
        Set the composition of a Layer with a Material
        """
        assert(isinstance(composition, Material))
        self.composition = composition
    
    def set_temperature_mode(self,  temperature_mode='adiabat',
                  temperatures=None, temperature_top=None):
        """
        Sets temperature of the layer by user-defined values or as an (modified) adiabat.
        temperature_mode is 'user_defined','adiabatic', or 'modified_adiabat'
        """
        assert(temperature_mode == 'user_defined' or temperature_mode ==
               'adiabat' or temperature_mode == 'modified_adiabat')
       
        self.temperature_mode= temperature_mode
        
        if temperature_mode == 'user_defined' or temperature_mode=='modified_adiabat':
            assert(len(temperatures) == len(self.depths))
            self.usertemperatures = temperatures
        else:
            self.usertemperatures = np.zeros_like(self.depths)

        if temperature_mode == 'adiabat' or temperature_mode=='modified_adiabat':
            self.temperature_top=temperature_top
        else:
            self.temperature_top=None


    def set_state(self, pressure_mode='selfconsistent',  pressures=None, pressure_top=None, gravity_bottom=None,n_max_iterations=50):
        """
        Sets the pressure and temperature of the layer by user-defined values are in a self-consistent fashion.
        temperature_mode is 'user_defined','adiabatic', or 'modified_adiabat'
        pressure_mode is 'user_defined' or 'selfconsistent'
        """
        assert(pressure_mode == 'user_defined' or pressure_mode == 'selfconsistent')
        assert(self.temperature_mode is not None)
        
        self.gravity_bottom = gravity_bottom

        if pressure_mode == 'user_defined':
            assert(len(pressures) == len(self.depths))
            self._pressures = pressures
            warnings.warn(
                "By setting the pressures in Layer it is unlikely to be self-consistent")
            self._temperatures = self._evaluate_temperature(self._pressures, self.temperature_top)

        if pressure_mode == 'selfconsistent':
            self.pressure_top = pressure_top
            ref_press = np.zeros_like(pressures)
            new_press = self.pressure_top + \
                (self.depths - min(self.depths)) * \
                2.e5  # initial pressure curve guess
            temperatures = self._evaluate_temperature(new_press, self.temperature_top)
            # Make it self-consistent!!!
            i = 0

            while i< n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(new_press, temperatures, gravity_bottom, pressure_top)
                temperatures = self._evaluate_temperature(new_press, self.temperature_top)
                rel_err = abs((max(ref_press) - max(new_press) )/ max(new_press))
                if self.verbose:
                    print("Iteration %i  maximum core pressure error between iterations: %e" % (i,rel_err))

                if rel_err < 1e-5:
                    break
            
            self._pressures = new_press
            self._temperatures = temperatures
            
        self.layer = []
        for l in range(len(self.depths)):
            self.layer.append(self.composition.copy())
            self.layer[l].set_state(
                self._pressures[l], self._temperatures[l])

    def _evaluate_temperature(self, pressures=None, temperature_top=None):
        if self.temperature_mode == 'adiabat' or self.temperature_mode == 'modified_adiabat' :
            adiabat = geotherm.adiabatic(
                pressures, temperature_top, self.composition)
        else:
            adiabat = np.zeros_like(self.depths)
        return adiabat + self.usertemperatures


    def _evaluate_eos(self, pressures, temperatures, gravity_bottom, pressure_top):
            [density] = self.composition.evaluate(['density'], pressures, temperatures)
            grav = self._compute_gravity(density , gravity_bottom)
            press =  self._compute_pressure( density, grav, pressure_top)
            return grav, press


    # Functions needed to compute self-consistent depths-pressures
    def _compute_gravity(self, density, gravity_bottom):
        """
        Computes the gravity of a layer
        """
        radii = self.radii[::-1]
        density = np.squeeze(density)[::-1]
        # Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(radii, density)
        # Numerically integrate Poisson's equation
        def poisson(p, x): return 4.0 * np.pi * \
            constants.G * rhofunc(x) * x * x
        grav = np.ravel(odeint(poisson, gravity_bottom * radii[0] * radii[0], radii))
        grav[:] = grav[:] / radii[:] / radii[:]
        if radii[0] ==0:
            grav[0]=0
        return grav[::-1]

    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
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
        radii=self.radii[::-1]
        density=self.density[::-1]
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
        radii=self.radii[::-1]
        density=self.density[::-1]
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
        dkappadP = np.gradient(kappa, edge_order=2) / np.gradient(self.pressure, edge_order=2)
        dphidz = np.gradient(phi, edge_order=2) / np.gradient(self.depths, edge_order=2) / self.gravity
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
        pressure : float
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
        temperature : float
            Temperature in [K].
        """
        return self._temperatures

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
        return np.array(
            [self.layer[i].internal_energy for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_gibbs for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_helmholtz for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_mass for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_volume for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].density for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_entropy for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].molar_enthalpy for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].isothermal_bulk_modulus for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].adiabatic_bulk_modulus for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].isothermal_compressibility for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].adiabatic_compressibility for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].shear_modulus for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].p_wave_velocity for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].bulk_sound_velocity for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].shear_wave_velocity for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].grueneisen_parameter for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].thermal_expansivity for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].heat_capacity_v for i in range(len(self.layer))])

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
        return np.array(
            [self.layer[i].heat_capacity_p for i in range(len(self.layer))])
