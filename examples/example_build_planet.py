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

import burnman
from burnman import planet
from burnman import mineral_helpers as helpers


def core_mass(density, radii):
    # returns a mass for a layer or sum of layers smaller than the whole planet

    rhofunc = UnivariateSpline(radii, density)
    mass = quad(lambda r: 4. * np.pi * rhofunc(r) * r * r,
                radii[0], radii[-1])[0]
    return mass

if __name__ == "__main__":

    # gravitational constant
    Radius_of_earth = 6371.e3
    # A basic set of EoS parameters for solid iron

    Core = [burnman.minerals.other.Liquid_Fe_Anderson(),3485e3]
    LM = burnman.minerals.SLB_2011.mg_bridgmanite()
    UM = burnman.minerals.SLB_2011.forsterite()
    mantle_rock = helpers.HelperLowHighPressureRockTransition(25.0e9, LM, UM)
    Mantle = [mantle_rock, 2886e3]
    Planet_radius = (Core[1]+Mantle[1])
    #temperatures = [300. for i in radii]

    compositions = [Core,Mantle]

    #Plan = planet.Planet(compositions, radii, temperatures, n_layers ,n_iterations)
    Plan = planet.Planet(compositions)

    #Plan.generate_profiles(radii, n_iterations)
    import matplotlib.gridspec as gridspec

    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
    plt.rc('font', family='sanserif')

    #Come up with axes for the final plot
    figure = plt.figure( figsize = (12,10) )
    ax1 = plt.subplot2grid( (5,3) , (0,0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid( (5,3) , (3,0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid( (5,3) , (4,0), colspan=3, rowspan=1)
         #Plot density, vphi, and vs for the planet.

    #Also plot a black line for the icb, cmb and upper mantle
    ylimits = [2., (Plan.densities[0]/1e3)+1.]
    ax1.plot(Plan.radii/1.e3,Plan.densities/1.e3,'k', linewidth = 2.,color='r')
    #ax1.plot( [trans_radii/1.e3, trans_radii/1.e3], ylimits, 'k', linewidth = 2.,color='g')

    ax1.set_ylabel("Density (kg/m$^3$)")

    #Make a subplot showing the calculated pressure profile
    ax2.plot( Plan.radii/1.e3, Plan.pressures/1.e9, 'k', linewidth=2.)
    ax2.set_ylabel("Pressure (GPa)")

    #Make a subplot showing the calculated gravity profile
    ax3.plot( Plan.radii/1.e3, Plan.gravity, 'k', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlabel("Radius (km)")

    plt.show()
