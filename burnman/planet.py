from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from burnman import averaging_schemes
from burnman import constants


class Planet(object):
    """
    A planet class that find a self-consistent planet.
    """

    class LayerBase(object):
        """
        Base class for a layer of a planet.
        """
        def __init__(self, name, rock, outer_radius, n_slices=100):
            self.name = name
            self.rock = rock
            self.inner_radius = None
            self.outer_radius = outer_radius
            self.thickness = None
            self.n_slices = n_slices
            self.n_start = None
            self.n_end = None
            self.mass = None
            assert(n_slices>=6)

    class Layer(LayerBase):
        """
        A layer with a constant temperature.
        """
        def __init__(self, name, rock, outer_radius, temperature = 300.0, n_slices = 100):
            Planet.LayerBase.__init__(self, name, rock, outer_radius, n_slices)
            self.constant_temperature = temperature

        def temperature(self, radius):
            return self.constant_temperature

    class LayerLinearTemperature(LayerBase):
        """
        A layer with a linear temperature profile.
        """
        def __init__(self, name, rock, outer_radius, temperature_inner, temperature_outer, n_slices = 100):
            Planet.LayerBase.__init__(self, name, rock, outer_radius, n_slices)
            self._temperature_inner = temperature_inner
            self._temperature_outer = temperature_outer

        def temperature(self, radius):
            return self._temperature_inner + (self._temperature_outer-self._temperature_inner)*(radius-self.inner_radius)/self.thickness

    def __init__(self, layers, n_max_iterations = 50, verbose = False):
        """
        Generate the planet based on the given layers (List of Layer)
        """

        # sort layers
        self.layers = sorted(layers, key=lambda x: x.outer_radius)

        # compute thickness, slices, indices, etc.
        self.radial_slices = np.array([])
        current_radius = 0.0
        n = 0
        for layer in self.layers:
            layer.inner_radius = current_radius
            slices = np.linspace(current_radius, layer.outer_radius, layer.n_slices)
            layer.thickness = layer.outer_radius - layer.inner_radius
            current_radius = layer.outer_radius
            slices[0] += 1 # make the boundaries one meter thick
            layer.n_start = n
            layer.n_end = n + layer.n_slices
            n += layer.n_slices
            self.radial_slices = np.append(self.radial_slices, slices)

        self.radial_slices[0] = 0.0  # overwrite inner-most radius to zero
        self.radius = current_radius

        # initial pressure guess: linear as a function of radius
        max_p = 350e9
        self.pressures = np.array([max_p * (1.0 - r / self.radius) for r in self.radial_slices])
        self.temperatures = [self.get_layer_by_radius(r).temperature(r) for r in self.radial_slices]
        self.densities = np.empty_like(self.pressures)
        self.gravity = np.empty_like(self.pressures)

        self.averaging_scheme = averaging_schemes.VoigtReussHill()

        last_center_pressure = 0.0
        for i in range(n_max_iterations):
            if verbose:
                print("on iteration %d" % (i+1))
            self._evaluate_eos()

            self._compute_gravity(self.densities, self.radial_slices)
            self._compute_pressure(self.densities, self.gravity, self.radial_slices)

            # compute relative error (pressure in the center of our planet)
            rel_err = abs(last_center_pressure - self.pressures[0]) / self.pressures[0]
            if verbose:
                print("  relative central core pressure error between iterations: %e" % rel_err)

            if rel_err < 1e-5:
                break
            last_center_pressure = self.pressures[0]

        self.mass = self._compute_mass()

        self.moment_of_inertia = self._compute_moment_of_inertia()
        self.moment_of_inertia_factor = self.moment_of_inertia / self.mass / self.radial_slices[-1] / self.radial_slices[-1]

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

    def _evaluate_eos(self):
        # evaluate each layer separately
        for layer in self.layers:
            mypressures = self.pressures[layer.n_start: layer.n_end]
            mytemperatures = self.temperatures[layer.n_start: layer.n_end]

            density = layer.rock.evaluate(['density'], mypressures, mytemperatures)

            self.densities[layer.n_start: layer.n_end] = density

    def _compute_gravity(self, density, radii):
        """
        Calculate the gravity of the planet, based on a density profile.  This integrates
        Poisson's equation in radius, under the assumption that the planet is laterally
        homogeneous.
        """

        start_gravity = 0.0
        for layer in self.layers:
            radii = self.radial_slices[layer.n_start: layer.n_end]
            density = self.densities[layer.n_start: layer.n_end]
            rhofunc = UnivariateSpline(radii, density)
            #Create a spline fit of density as a function of radius

            #Numerically integrate Poisson's equation
            poisson = lambda p, x: 4.0 * np.pi * constants.G * rhofunc(x) * x * x
            grav = np.ravel(odeint( poisson, start_gravity, radii))
            start_gravity = grav[-1]
            self.gravity[layer.n_start: layer.n_end] = grav

        # we need to skip scaling gravity[0]
        self.gravity[1:] = self.gravity[1:]/self.radial_slices[1:]/self.radial_slices[1:]

    def _compute_pressure(self, density, gravity, radii):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        """

        start_pressure = 0.0
        for layer in self.layers[::-1]:
            radii = self.radial_slices[layer.n_start: layer.n_end]
            density = self.densities[layer.n_start: layer.n_end]
            gravity = self.gravity[layer.n_start: layer.n_end]

            # convert radii to depths
            depths = radii[-1]-radii

            # Make a spline fit of density as a function of depth
            rhofunc = UnivariateSpline(depths[::-1], density[::-1])
            # Make a spline fit of gravity as a function of depth
            gfunc = UnivariateSpline(depths[::-1], gravity[::-1])

            # integrate the hydrostatic equation
            pressure = np.ravel(odeint((lambda p, x : gfunc(x)* rhofunc(x)), start_pressure, depths[::-1]))
            start_pressure = pressure[-1]

            self.pressures[layer.n_start: layer.n_end] = pressure[::-1]

    def _compute_mass( self):
        """
        calculates the mass of the entire planet [kg]
        """
        mass = 0.0
        for layer in self.layers:
            radii = self.radial_slices[layer.n_start: layer.n_end]
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
