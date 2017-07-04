# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function

import os
import sys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import matplotlib.pyplot as plt
import numpy as np

import burnman
import burnman.mineral_helpers as helpers

if __name__ == "__main__":

    '''
    example_build_planet
    --------------------

    For Earth we have well-constrained one-dimensional density models.  This allows us to
    calculate pressure as a function of depth.  Furthermore, petrologic data and assumptions
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

    This example allows the user to define layers of planets of known outer radius and self-
    consistently solve for the density, pressure and gravity profiles. The calculation will
    iterate until the difference between central pressure calculations are less than 1e-5.
    The planet class in BurnMan (../burnman/planet.py) allows users to call multiple
    properties of the model planet after calculations, such as the mass of an individual layer,
    the total mass of the planet and the moment if inertia. See planets.py for information
    on each of the parameters which can be called.

    *Uses:*

    * :doc:`mineral_database`
    * :class:`burnman.planet.Planet`
    * :class: `burnman.layer.Layer`
    
    *Demonstrates:*
    
    * setting up a planet
    * computing its self-consistent state
    * computing various parameters for the planet
    * seismic comparison
    '''

    # FIRST: we must define the composition of the planet as individual layers.
    # A layer is defined by 4 parameters: Name, min_depth, max_depth,and number of slices within the layer.
    # Separately the composition and the temperature_mode need to set.
    radius_planet = 6371.e3
    # inner_core
    inner_core = burnman.Layer("inner core", radius_planet=radius_planet,
        max_depth=6371.e3, min_depth=5151.e3, n_slices=10)
    inner_core.set_composition(burnman.minerals.other.Fe_Dewaele())
    
    # The minerals that make up our core do not currently implement the thermal equation of state, so we will define the temperature with the Brown & Shankland geotherm.
    inner_core.set_temperature_mode( 'user_defined',
        burnman.geotherm.brown_shankland(inner_core.depths))

    # outer_core
    outer_core = burnman.Layer( "outer core", radius_planet=radius_planet,
        max_depth=5151.e3, min_depth=2891.e3, n_slices=10)
    outer_core.set_composition(burnman.minerals.other.Liquid_Fe_Anderson())
    # The minerals that make up our core do not currently implement the thermal equation of state, so we will define the temperature with the Brown & Shankland geotherm.
    outer_core.set_temperature_mode('user_defined',
        burnman.geotherm.brown_shankland(outer_core.depths))

    # Next the Mantle.
    lower_mantle = burnman.Layer( "lower_mantle", radius_planet=radius_planet,
        max_depth=2891.e3, min_depth=660.e3, n_slices=10)
    lower_mantle.set_composition(burnman.minerals.SLB_2011.mg_bridgmanite())
    lower_mantle.set_temperature_mode('adiabat')
    upper_mantle = burnman.Layer( "upper_mantle", radius_planet=radius_planet,
        max_depth=660e3, min_depth=0., n_slices=10)
    upper_mantle.set_composition(burnman.minerals.SLB_2011.forsterite())
    upper_mantle.set_temperature_mode('adiabat')



    # Now we calculate the planet.
    Plan = burnman.Planet('earth_like',
                          [inner_core, outer_core, lower_mantle, upper_mantle],
                          potential_temperature=1200, verbose=True)
    # Here we compute its state. Go BurnMan Go!
    # (If we were to change composition of one of the layers, we would have to
    # recompute the state)
    Plan.set_state()

    # Now we output the mass of the planet and moment of inertia
    print()
    print("mass/Earth= %.3f, moment of inertia= %.3f" %
          (Plan.mass / 5.97e24, Plan.moment_of_inertia_factor))

    # And here's the mass of the individual layers:
    for layer in Plan.layers:
        print("%s mass fraction of planet %.3f" %
              (layer.name, layer.mass / Plan.mass))
    print()

    # Let's get PREM to compare everything to as we are trying
    # to imitate Earth
    prem = burnman.seismic.PREM()
    premradii = 6371.e3 - prem.internal_depth_list()
    [premdensity, prempressure, premgravity] = prem.evaluate(
        ['density', 'pressure', 'gravity'])

    # Now let's plot everything up

    # Do you like pretty, but slow plotting figures? Try Burnman's patented pretty_plot package.
    # otherwise this line can be commented out
    #burnman.tools.pretty_plot()

    # Come up with axes for the final plot

    figure = plt.figure(figsize=(11, 7))
    figure.suptitle(
        'Your planet is %.3f Earth Masses with Average Density of %.1f kg/m$^3$' %
        ((Plan.mass /5.97e24),(Plan.mass /(4. /3. * np.pi * Plan.radius_planet *
              Plan.radius_planet * Plan.radius_planet))), fontsize=20)

    ax1 = plt.subplot2grid((6, 3), (0, 0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid((6, 3), (3, 0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid((6, 3), (4, 0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid((6, 3), (5, 0), colspan=3, rowspan=1)

    ax1.plot(Plan.radii / 1.e3, Plan.density / 1.e3, 'k', linewidth=2.,
        label='your planet')
    ax1.plot( premradii / 1.e3, premdensity / 1.e3, '--k', linewidth=1.,
        label='PREM')
    ax1.set_ylim(0., (max(Plan.density) / 1.e3) + 1.)
    ax1.set_xlim(0., max(Plan.radii) / 1.e3)
    ax1.set_xticklabels([])
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")
    ax1.legend()

    # Make a subplot showing the calculated pressure profile
    ax2.plot(Plan.radii / 1.e3, Plan.pressure / 1.e9, 'b', linewidth=2.)
    ax2.plot(premradii / 1.e3, prempressure / 1.e9, '--b', linewidth=1.)
    ax2.set_ylim(0., (max(Plan.pressure) / 1e9) + 10.)
    ax2.set_xlim(0., max(Plan.radii) / 1.e3)
    ax2.set_xticklabels([])
    ax2.set_ylabel("Pressure (GPa)")

    # Make a subplot showing the calculated gravity profile
    ax3.plot(Plan.radii / 1.e3, Plan.gravity, 'g', linewidth=2.)
    ax3.plot(premradii / 1.e3, premgravity, '--g', linewidth=1.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0., max(Plan.radii) / 1.e3)
    ax3.set_ylim(0., max(Plan.gravity) + 0.5)
    ax3.set_xticklabels([])

    # Make a subplot showing the calculated temperature profile
    ax4.plot(Plan.radii / 1.e3, Plan.temperature, 'r', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0., max(Plan.radii) / 1.e3)
    ax4.set_ylim(0., max(Plan.temperature) + 100)

    plt.show()
