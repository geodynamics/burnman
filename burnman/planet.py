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

class Planet(Material):
    """
    A planet class that find a self-consistent planet.
    """

    def __init__(self, name, layers, potential_temperature =0., n_max_iterations = 50, verbose = False):
        """
        Generate the planet based on the given layers (List of Layer)
        """
        Material.__init__(self)
        # sort layers
        self.layers = sorted(layers, key=lambda x: x.min_depth)
        #assert layers attach to one another
        if len(self.layers)>1:
            for l in range(1,len(self.layers)):
                assert(self.layers[l].outer_radius==self.layers[l-1].inner_radius)
        

        self.name = name
        
   
   
        self.depths=self.evaluate_planet(['depths'])
        self.min_depth = min(self.depths)
        self.max_depth = max(self.depths)
        self.n_slices = len(self.depths)
        self.radius_planet = max(self.depths)
        self.radii = self.radius_planet - self.depths
        
        for layer in self.layers:
            layer.n_start = np.where(self.depths == layer.min_depth)[0][-1]
            layer.n_end = np.where(self.depths == layer.max_depth)[0][0]+1
        self.potential_temperature = potential_temperature
        
        self.verbose = verbose


    def get_layer(self, name):
        for layer in self.layers:
            if layer.name == name:
                return layer
        raise LookupError()

    def get_layer_by_radius(self, radius):
        for layer in self.layers:
            if layer.outer_radius >= radius:
                return layer
        raise LookupError()



    def evaluate_planet(self,properties):
        all=None
        for prop in properties:
            oneprop=None
            for layer in self.layers:
                if oneprop is None:
                    oneprop = getattr(layer, prop)
                else:
                    oneprop = np.append(oneprop,getattr(layer, prop))
            if all is None:
                all=np.array(oneprop)
            else:
                np.append(all,oneprop, axis=0)
        return all




    def set_state(self,pressure_mode='selfconsistent',  pressures=None, pressure_top=0., gravity_bottom=0.,n_max_iterations=50):
        """
        Sets the pressure of the planet by user-defined values are in a self-consistent fashion.
        pressure_mode is 'user_defined' or 'selfconsistent'.
        The default for the planet is self-consistent, with zero pressure at the surface and zero pressure at the center.
        """
        assert(pressure_mode == 'user_defined' or pressure_mode == 'selfconsistent')
        for layer in self.layers:
            assert(layer.temperature_mode is not None)
        
        self.gravity_bottom = gravity_bottom

        if pressure_mode == 'user_defined':
            assert(len(pressures) == len(self.depths))
            self._pressures = pressures
            warnings.warn(
                "By setting the pressures in Planet it is unlikely to be self-consistent")
            self._temperatures = self._evaluate_temperature(self._pressures, self.potential_temperature)

        if pressure_mode == 'selfconsistent':
            self.pressure_top = pressure_top
            ref_press = np.zeros_like(pressures)
            new_press = self.pressure_top + \
                (self.depths - min(self.depths)) * \
                2.e5  # initial pressure curve guess
            temperatures = self._evaluate_temperature(new_press, self.potential_temperature)
            # Make it self-consistent!!!
            i = 0

            while i< n_max_iterations:
                i += 1
                ref_press = new_press
                new_grav, new_press = self._evaluate_eos(new_press, temperatures, gravity_bottom, pressure_top)
                temperatures = self._evaluate_temperature(new_press, self.potential_temperature)
                rel_err = abs((max(ref_press) - max(new_press) )/ max(new_press))
                if self.verbose:
                    print("Iteration %i  maximum core pressure error between iterations: %e" % (i,rel_err))

                if rel_err < 1e-5:
                    break
            
            self._pressures = new_press
            self._temperatures = temperatures
            self._gravity = new_grav
            
        
        for layer in self.layers:
            layer.layer = []
            layer._pressures = self._pressures[layer.n_start: layer.n_end]
            layer._temperatures = self._temperatures[layer.n_start: layer.n_end]
            for l in range(len(layer.depths)):
                layer.layer.append(layer.composition.copy())
                layer.layer[l].set_state(
                layer._pressures[l], self._temperatures[l])


    def _evaluate_eos(self, pressures,temperatures, gravity_bottom, pressure_top):
        density = self._evaluate_density(pressures,temperatures)
        grav = self._compute_gravity(density , gravity_bottom)
        press =  self._compute_pressure( density, grav, pressure_top)
        return grav, press
    
    
    def _evaluate_density(self, pressures,temperatures):
        density = []
        for layer in self.layers:
            density.append(layer.composition.evaluate(['density'], pressures[layer.n_start:layer.n_end], temperatures[layer.n_start:layer.n_end]))
        return np.squeeze(np.hstack(density))
        



    def _evaluate_temperature(self,pressures, temperature_top):
        temps= []
        for layer in self.layers:
            temps.append(layer._evaluate_temperature(pressures[layer.n_start:layer.n_end], temperature_top))
            temperature_top=temps[-1]
        return np.hstack(np.squeeze(temps))

    def _compute_gravity(self, density, gravity_bottom):
        """
        Calculate the gravity of the planet, based on a density profile.  This integrates
        Poisson's equation in radius, under the assumption that the planet is laterally
        homogeneous.
        """

        start_gravity = gravity_bottom
        grav = []
        for layer in self.layers[::-1]:
            grav.extend(layer._compute_gravity(density[layer.n_start: layer.n_end], start_gravity)[::-1])
            start_gravity = grav[-1]
        return np.array(grav)[::-1]


    def _compute_pressure(self, density, gravity, pressure_top):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        """

        start_pressure = pressure_top
        press=[]
        for layer in self.layers:
            press.extend(layer._compute_pressure(density[layer.n_start: layer.n_end], gravity[layer.n_start: layer.n_end], start_pressure))
        return np.array(press)
                         

    @property
    def mass( self):
        """
        calculates the mass of the entire planet [kg]
        """
        mass = 0.0
        for layer in self.layers:
            mass += layer.mass
        return mass

    @property
    def moment_of_inertia( self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment = 0.0
        for layer in self.layers:
            moment += layer.moment_of_inertia
        return moment

    @property
    def moment_of_inertia_factor( self):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        moment_factor = self.moment_of_inertia/self.mass/self.radius_planet/self.radius_planet
        return moment_factor

    @property
    def gravity(self):
        """
        Returns gravity of the layer [m s^(-2)]
        """
        return self._compute_gravity(self.density)

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
        return self.evaluate_planet(['density'])

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
