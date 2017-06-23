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
import warnings

from .material import Material, material_property
from .mineral import Mineral
from .composite import Composite
from .seismic import Seismic1DModel


class Layer(object):
    """
    A planetary layer class
    """
    def __init__(self, name=None, radius_planet=None, min_depth=None, max_depth=None, n_slices = None):
        
        self.name = name

        self.min_depth = min_depth
        self.max_depth = max_depth
        self.thickness = self.max_depth-self.min_depth
        self.n_slices = n_slices
        self.depths = np.linspace(min_depth, max_depth, n_slices)
        self.radius_planet = radius_planet

        # mass
        # moment_inertia
        
    def set_composition(self,composition):
        assert(isinstance(composition, Material) or isinstance(composition, Composite) or isinstance(composition, Mineral)  )
        self.composition = composition



    
    def set_state(self, pressure_mode = 'selfconsistent',temperature_mode='adiabat',pressures = None, temperatures = None, pressure_top= None, gravity_bottom=None, temperature_top=None):
        """
            temperature_mode is 'user_defined','adiabatic', or 'modified_adiabat'
            pressure_mode is 'user_defined' or 'selfconsistent'
        """
        if temperature_mode == 'user_defined':
            assert(len(temperatures)==len(self.depths))
            self.temperatures = temperatures
        if pressure_mode == 'user_defined':
            assert(len(pressures)==len(self.depths))
            self.pressures = pressures
            warnings.Warn("By setting the pressures in Layer it is unlikely to be self-consistent")
            if temperature_mode == 'adiabat':
                self.temperatures = burnman.geotherm.adiabatic(pressures, T0, rock)
            if  temperature_mode == 'modified_adiabat':
                self.temperatures = burnman.geotherm.adiabatic(pressures, T0, rock)+temperatures


        if pressure_mode == 'selfconsistent':
            self.pressure_top = pressure_top
            self.gravity_bottom = gravity_bottom
            
            ref_press=np.zeros_like(pressures)
            new_press=self.pressure_top + (self.depths-min(self.depths))*1.e5 # initial pressure curve
            i=0
            # Make it self-consistent!!!
            while np.sum(np.abs(new_press-ref_press))>1.e6*len(self.depths):
                i=i+1

                if temperature_mode == 'adiabat':
                    temperatures = burnman.geotherm.adiabatic(new_press, T0, rock)
                if  temperature_mode == 'modified_adiabat':
                    temperatures = burnman.geotherm.adiabatic(new_press, T0, rock)+temperatures
                mat_rho = self.compositon.evaluate(['density'], new_press, temperatures)
                grav=_compute_gravity(mat_rho[::-1]) #values with radius
                ref_press=new_press
                new_press=_compute_pressure(mat_rho[::-1],grav)#values with radius
                if i>50:
                    print('converged to ', str((new_press-ref_press/ref_press)*100.))
                    break
            self.pressures=new_press
            self.temperatures=temperatures
            self.gravity = grav
                          
    ##### Functions needed to compute self-consistent depths-pressures
    def _compute_gravity(self,density):
            self.radii=self.radius_planet - depths
            radii=self.radii[::-1]
            density =density[::-1]
            # Create a spline fit of density as a function of radius
            rhofunc = UnivariateSpline(radii, density)
            # Numerically integrate Poisson's equation
            poisson = lambda p, x: 4.0 * np.pi * G * rhofunc(x) * x * x
            grav = np.ravel(odeint(poisson, self_gravity_bottom*radii[0]*radii[0], radii))
            grav[:] = grav[:] / radii[:] / radii[:]
            return grav[::-1]
            
    def _compute_pressure(self,densities, gravity):
            # Calculate the pressure profile based on density and gravity.  This integrates
            # the equation for hydrostatic equilibrium  P = rho g z.
            
            # convert radii to depths
            depth = self.depths
            # Make a spline fit of density as a function of depth
            rhofunc = UnivariateSpline(depth, density)
            # Make a spline fit of gravity as a function of depth
            gfunc = UnivariateSpline(depth, gravity)
            
            # integrate the hydrostatic equation
            pressure = np.ravel(
                                odeint((lambda p, x: gfunc(x) * rhofunc(x)), self.pressure_top, depth))
            return pressure


    '''
    def _compute_mass( self):
        """
        calculates the mass of the layer [kg]
        """
        mass = 0.0
        radii = radius_planet-self.depths
        density = self.densities[layer.n_start: layer.n_end]
        rhofunc = UnivariateSpline(radii, density)
        layer.mass = quad(lambda r : 4*np.pi*rhofunc(r)*r*r,
                            radii[0], radii[-1])[0]
            mass += layer.mass
        return mass

    def _compute_moment_of_inertia( self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment = 0.0
        for layer in self.layers:
            radii = self.radial_slices[layer.n_start: layer.n_end]
            density = self.densities[layer.n_start: layer.n_end]
            rhofunc = UnivariateSpline(radii, density)
            moment += quad(lambda r : 8.0/3.0*np.pi*rhofunc(r)*r*r*r*r,
                           radii[0], radii[-1])[0]
        return moment

    '''

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
