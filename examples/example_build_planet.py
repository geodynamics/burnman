# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


'''
example_build_planet
--------------------

For Earth we have well-constrained one-dimensional density models.  This allows us to
calculate pressure as a funcion of depth.  Furthermore, petrologic data and assumptions
regarding the convective state of the planet allow us to estimate the temperature.

For planets other than Earth we have much less information, and in particular we
know almost nothing about the pressure and temperature in the interior.  Instead, we tend
to have measurements of things like mass, radius, and moment-of-inertia.  We would like
to be able to make a model of the planet's interior that is consistent with those
measurements.

However, there is a difficulty with this.  In order to know the density of the planetary
material, we need to know the pressure and temperature.  In order to know the pressure,
we need to know the gravity profile.  And in order to the the gravity profile, we need
to know the density.  This is a nonlinear problem which requires us to iterate to find
a self-consistent solution.

Here we show an example that does this, using the planet Mercury as motivation.


*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :func:`burnman.material.Material.evaluate`
'''
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

import burnman
import burnman.minerals as minerals


if __name__ == "__main__":

    # gravitational constant
    G = 6.67e-11

    # A basic set of EoS parameters for solid iron
    class iron(burnman.Mineral):

        def __init__(self):
            # Parameters for gamma - Fe are from Tsujino et al. 2013
            # G_0 G_Prime0 from Mao et al. 2001 (fig. 3)
            self.params = {
                'equation_of_state': 'slb3',
                'T_0': 1273.,
                'V_0': 7.381e-06,
                'K_0': 111.5e9,
                'Kprime_0': 5.2,
                'G_0': 83.2e9,  # Shear modulus and derivative from Gleason and Mao, 2013
                'Gprime_0': 2.04,
                'molar_mass': 55.845 / 1000.,
                'n': 1,
                'Debye_0': 340.,
                'grueneisen_0': 2.28,
                'q_0': 0.21,
                'eta_s_0': 2.0  # Wholly invented value
            }
            burnman.Mineral.__init__(self)

    # Create a class for Mercury that will do the heavy lifting.
    class Mercury(object):

        def __init__(self, n_slices):
            # The constructor takes the number of depth slices which will
            # be used for the calculation.  More slices will generate more
            # accurate profiles, but it will take longer.

            self.cmb = 2020.e3  # Guess for the radius of the core-mantle-boundary
            self.outer_radius = 2440.e3  # Outer radius of the planet

            self.radii = np.linspace(
                0.e3, self.outer_radius, n_slices)  # Radius list
            self.pressures = np.linspace(
                35.0e9, 0.0, n_slices)  # initial guess at pressure profile
            self.temperatures = np.ones_like(
                self.pressures) * 1000.  # Assume an isothermal interior (not a great assumption, but we do this for simplicity's sake).

            # The planet will be represented by a two layer model, mantle and core.
            # The top layer will be a composite of 80% perovskite and 20%
            # periclase.
            amount_olivine = 0.8
            self.mantle = burnman.Composite([minerals.SLB_2011.forsterite(),
                                             minerals.SLB_2011.enstatite()],
                                            [amount_olivine, 1.0 - amount_olivine])
             # The core will be represented by solid iron.
            self.core = iron()

        def generate_profiles(self, n_iterations):
            # Generate the density, gravity, pressure, vphi, and vs profiles for the planet.
            # The n_iterations parameter sets how many times to iterate over density, pressure,
            # and gravity.  Empirically, five-ish iterations is adequate.  After the iterations,
            # this also calculates mass and moment of inertia of the planet.

            for i in range(n_iterations):
                self.densities, self.bulk_sound_speed, self.shear_velocity = self._evaluate_eos(
                    self.pressures, self.temperatures, self.radii)
                self.gravity = self._compute_gravity(
                    self.densities, self.radii)
                self.pressures = self._compute_pressure(
                    self.densities, self.gravity, self.radii)

            self.mass = self._compute_mass(self.densities, self.radii)
            self.moment_of_inertia = self._compute_moment_of_inertia(
                self.densities, self.radii)
            self.moment_of_inertia_factor = self.moment_of_inertia / \
                self.mass / self.outer_radius / self.outer_radius

        def _evaluate_eos(self, pressures, temperatures, radii):
            # Evaluates the equation of state for each radius slice of the model.
            # Returns density, bulk sound speed, and shear speed.

            rho = np.empty_like(radii)
            bulk_sound_speed = np.empty_like(radii)
            shear_velocity = np.empty_like(radii)

            for i in range(len(radii)):
                density = vs = vphi = 0.

                if radii[i] > self.cmb:
                    density, vs, vphi = self.mantle.evaluate(
                        ['density', 'v_s', 'v_phi'], [pressures[i]], [temperatures[i]])
                else:
                    density, vs, vphi = self.core.evaluate(
                        ['density', 'v_s', 'v_phi'], [pressures[i]], [temperatures[i]])

                rho[i] = density
                bulk_sound_speed[i] = vphi
                shear_velocity[i] = vs

            return rho, bulk_sound_speed, shear_velocity

        def _compute_gravity(self, density, radii):
            # Calculate the gravity of the planet, based on a density profile.  This integrates
            # Poisson's equation in radius, under the assumption that the planet is laterally
            # homogeneous.

            # Create a spline fit of density as a function of radius
            rhofunc = UnivariateSpline(radii, density)

            # Numerically integrate Poisson's equation
            poisson = lambda p, x: 4.0 * np.pi * G * rhofunc(x) * x * x
            grav = np.ravel(odeint(poisson, 0.0, radii))
            grav[1:] = grav[1:] / radii[1:] / radii[1:]
            grav[
                0] = 0.0  # Set it to zero a the center, since radius = 0 there we cannot divide by r^2
            return grav

        def _compute_pressure(self, density, gravity, radii):
            # Calculate the pressure profile based on density and gravity.  This integrates
            # the equation for hydrostatic equilibrium  P = rho g z.

            # convert radii to depths
            depth = radii[-1] - radii

            # Make a spline fit of density as a function of depth
            rhofunc = UnivariateSpline(depth[::-1], density[::-1])
            # Make a spline fit of gravity as a function of depth
            gfunc = UnivariateSpline(depth[::-1], gravity[::-1])

            # integrate the hydrostatic equation
            pressure = np.ravel(
                odeint((lambda p, x: gfunc(x) * rhofunc(x)), 0.0, depth[::-1]))
            return pressure[::-1]

        def _compute_mass(self, density, radii):
            # Returns a list of moments of inertia of the planet [kg m^2]
            rhofunc = UnivariateSpline(radii, density)
            mass = quad(lambda r: 4 * np.pi * rhofunc(r) * r * r,
                        radii[0], radii[-1])[0]
            return mass

        def _compute_moment_of_inertia(self, density, radii):
            # Returns the moment of inertia of the planet [kg m^2]

            rhofunc = UnivariateSpline(radii, density)
            moment = quad(
                lambda r: 8.0 / 3.0 * np.pi * rhofunc(r) * r * r * r * r,
                radii[0], radii[-1])[0]
            return moment

    # Here we actually do the interation.  We make an instance
    # of our Mercury planet, then call generate_profiles.
    # Emprically, 300 slices and 5 iterations seem to do
    # a good job of converging on the correct profiles.
    n_slices = 300
    n_iterations = 5
    merc = Mercury(n_slices)
    merc.generate_profiles(n_iterations)

    # These are the actual observables
    # from the model, that is to say,
    # the total mass of the planet and
    # the moment of inertia factor,
    # or C/MR^2
    observed_mass = 3.02e23
    observed_moment = 0.346  # From Margot. et al, 2012

    print(("Total mass of the planet: %.2e, or %.0f%% of the observed mass" %
          (merc.mass, merc.mass / observed_mass * 100.)))
    print(("Moment of inertia factor of the planet: %.3g, or %0.f%% of the observed factor" %
          (merc.moment_of_inertia_factor, merc.moment_of_inertia_factor / observed_moment * 100.)))

    # As we can see by running this, the calculated mass of the planet is much too large.
    # One could do a better job of fitting this by using a more complicated interior model,
    # with a liquid outer core, light alloying elements in the core, and a more realistic
    # temperature profile.  That, however, is outside of the scope of this
    # example.

    import matplotlib.gridspec as gridspec

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
    plt.rc('font', family='sans-serif')

    # Come up with axes for the final plot
    figure = plt.figure(figsize=(12, 10))
    ax1 = plt.subplot2grid((5, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((5, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((5, 3), (4, 0), colspan=3, rowspan=1)

    # Plot density, vphi, and vs for the planet.
    ax1.plot(merc.radii / 1.e3, merc.densities /
             1.e3, label=r'$\rho$', linewidth=2.)
    ax1.plot(merc.radii / 1.e3, merc.bulk_sound_speed /
             1.e3, label=r'$V_\phi$', linewidth=2.)
    ax1.plot(merc.radii / 1.e3, merc.shear_velocity /
             1.e3, label=r'$V_S$', linewidth=2.)

    # Also plot a black line for the CMB
    ylimits = [3., 10.]
    ax1.plot([merc.cmb / 1.e3, merc.cmb / 1.e3], ylimits, 'k', linewidth=6.)

    ax1.legend()
    ax1.set_ylabel("Velocities (km/s) and Density (kg/m$^3$)")

    # Make a subplot showing the calculated pressure profile
    ax2.plot(merc.radii / 1.e3, merc.pressures / 1.e9, 'k', linewidth=2.)
    ax2.set_ylabel("Pressure (GPa)")

    # Make a subplot showing the calculated gravity profile
    ax3.plot(merc.radii / 1.e3, merc.gravity, 'k', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlabel("Radius (km)")

    plt.show()
