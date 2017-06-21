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

from .material import Material
from .mineral import Mineral
from .composite import Composite
from .seismic import Seismic1DModel


class Layer(object):
    """
    A planetary layer class
    """
    def __init__(self, name=None,  min_depth=None, max_depth=None, n_slices = None):
        
        self.name = name

        self.min_depth = min_depth
        self.max_depth = max_depth
        self.thickness = self.max_depth-self.min_depth
        self.n_slices = n_slices
        self.depths = np.arange(min_depth, max_depth, n_slices)

    
        self.pressures = None
        self.temperatures = None
        self.composition = None
        #planet radius
        # mass
        # moment_inertia
        
    def set_composition(self,composition):
        assert(isinstance(composition, Material) or isinstance(composition, Composite) or isinstance(composition, Mineral)  )
        self.composition = composition
    
    def set_pressures(self, seismicmodel_or_array):
        if isinstance(seismicmodel_or_array, Seismic1DModel):
            seismicmodel = seismicmodel_or_array
            self.pressures = seismicmodel.get_pressures(self.depths)
        else:
            self.pressures = seismicmodel_or_array
        warnings.Warn("By setting the pressures in your layer by hand, you are making the")

    def set_temperature(self, temperatures)

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
