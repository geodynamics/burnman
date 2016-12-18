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
    def __init__(self, compositions,radii = None, temperatures = None, n_layers = None,n_iterations = None):
        """
        Generate the density, gravity, pressure, vphi, and vs profiles for the planet.
        """
        if n_layers == None:
            n_layers = 1000
        if n_iterations == None:
            n_iterations = 7
        if radii == None:
            Planet_radius = sum([k[1] for k in compositions])
            self.radii = np.linspace(0., Planet_radius, n_layers)
        if temperatures == None:
            self.temperatures = [300. for i in self.radii]

        self.pressures = np.linspace(350.0e9, 0.0, len(self.radii)) # initial guess at pressure profile
        self.compositions = compositions

        averaging_scheme = averaging_schemes.VoigtReussHill()

        for i in range(n_iterations):
            print("on iteration #", str(i+1)+"/"+str(n_iterations))
            self.densities, self.bulk_sound_speed, self.shear_velocity = self._evaluate_eos(self.compositions,self.pressures, self.temperatures, self.radii,averaging_scheme)
            self.gravity = self._compute_gravity(self.densities, self.radii)
            self.pressures = self._compute_pressure(self.densities, self.gravity, self.radii)

        self.mass = self._compute_mass(self.densities, self.radii)

        self.moment_of_inertia = self._compute_moment_of_inertia(self.densities, self.radii)
        self.moment_of_inertia_factor = self.moment_of_inertia / self.mass / self.radii[-1] / self.radii[-1]

    def _evaluate_eos(self, compositions, pressures, temperatures, radii,averaging_scheme):
        """
        Evaluates the equation of state for each radius slice of the model.
        Returns density, bulk sound speed, and shear speed.
        """
        rho = np.empty_like(radii)
        bulk_sound_speed = np.empty_like(radii)
        shear_velocity = np.empty_like(radii)
        for i in range(len(radii)):
            dummy_depth = 0.
            for j in range(len(compositions)):
                dummy_depth+=compositions[j][1]

                if radii[i] <= dummy_depth:
                    rock = compositions[j][0]
                    break

            density, vs, vphi = rock.evaluate(['density', 'v_s', 'v_phi'], np.array([pressures[i]]),np.array([temperatures[i]]))

            rho[i] = density
            bulk_sound_speed[i] = vphi
            shear_velocity[i] = vs


        return rho, bulk_sound_speed, shear_velocity


    def _compute_gravity(self, density, radii):
        """
        Calculate the gravity of the planet, based on a density profile.  This integrates
        Poisson's equation in radius, under the assumption that the planet is laterally
        homogeneous.
        """
        #Create a spline fit of density as a function of radius
        rhofunc = UnivariateSpline(radii, density )

        #Numerically integrate Poisson's equation
        poisson = lambda p, x : 4.0 * np.pi * constants.G * rhofunc(x) * x * x
        grav = np.ravel(odeint( poisson, 0.0, radii ))
        grav[1:] = grav[1:]/radii[1:]/radii[1:]
        grav[0] = 0.0 #Set it to zero a the center, since radius = 0 there we cannot divide by r^2
        return grav

    def _compute_pressure(self, density, gravity, radii):
        """
        Calculate the pressure profile based on density and gravity.  This integrates
        the equation for hydrostatic equilibrium  P = rho g z.
        """
        #convert radii to depths
        depth = radii[-1]-radii

        #Make a spline fit of density as a function of depth
        rhofunc = UnivariateSpline( depth[::-1], density[::-1] )
        #Make a spline fit of gravity as a function of depth
        gfunc = UnivariateSpline( depth[::-1], gravity[::-1] )

        #integrate the hydrostatic equation
        pressure = np.ravel(odeint( (lambda p, x : gfunc(x)* rhofunc(x)), 0.0,depth[::-1]))
        return pressure[::-1]

    def _compute_mass( self, density, radii):
        """
        calculates the mass of the entire planet [kg]
        """
        rhofunc = UnivariateSpline(radii, density )
        mass = quad( lambda r : 4*np.pi*rhofunc(r)*r*r,
                                 radii[0], radii[-1] )[0]
        return mass


    def _compute_moment_of_inertia( self, density, radii):
        """
        #Returns the moment of inertia of the planet [kg m^2]
        """
        rhofunc = UnivariateSpline(radii, density )
        moment = quad( lambda r : 8.0/3.0*np.pi*rhofunc(r)*r*r*r*r,
                                 radii[0], radii[-1] )[0]
        return moment