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

import matplotlib.pyplot as plt
import numpy as np

import burnman
import burnman.planet as planet
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
    * :class:`burnman.planet`
    '''

    #FIRST: we must define the composition of the planet as individual layers.
    #A layer is defined by 5 parameters: Name, mineral or rock composition, outer radius,
    #temperature profile and number of slices within the layer.

    #The minerals that make up our core do not currently implement the thermal equation of state, so we will
    #define the temperature as None. For an isothermal profile this can be change to a singular value (e.g. 300)

    inner_core = planet.Planet.Layer("inner core", burnman.minerals.other.Fe_Dewaele(), 1220e3, n_slices = 100)
    outer_core = planet.Planet.Layer("outer core", burnman.minerals.other.Liquid_Fe_Anderson(), 3485e3, n_slices = 500)

    #Next the Mantle.
    LM = burnman.minerals.SLB_2011.mg_bridgmanite()
    UM = burnman.minerals.SLB_2011.forsterite()

    #We don't want to have an ad hoc transition from the upper mantle to lower mantle, so we'll stitch them
    #together with a helper to switch between them at 25 GPa. So we define a class of "mantle":

    class mantle_rock(helpers.HelperLowHighPressureRockTransition):
        def __init__(self):
            helpers.HelperLowHighPressureRockTransition.__init__(self,25.e9,UM,LM)

    #Isothermal mantles are boring, here we define a mantle that has a linear temperature profile, all we need
    #
    mantle = planet.Planet.LayerLinearTemperature("mantle", mantle_rock(), 6371.e3, 3000.0, 300.0, 200)

    #Now we calculate the planet. Go BurnMan Go!
    Plan = planet.Planet([inner_core, outer_core, mantle], verbose=True)

    #Now we output the mass of the planet and moment of inertia
    print()
    print("mass/Earth= %.3f, moment of inertia= %.3f" % (Plan.mass/5.97e24, Plan.moment_of_inertia_factor))

    #And here's the mass of the individual layers:
    for layer in Plan.layers:
        print("%s mass fraction of planet %.3f" %(layer.name, layer.mass/Plan.mass))
    print()

    #Now let's plot everything up

    # Do you like pretty, but slow plotting figures? Try Burnman's patented pretty_plot package.
    # otherwise this line can be commented out
    #burnman.tools.pretty_plot()

    #Come up with axes for the final plot

    figure = plt.figure(figsize = (12,10))
    figure.suptitle('Your planet is %.3f Earth Masses with Average Density of %.1f kg/m$^3$' %((Plan.mass/5.97e24), \
                    (Plan.mass/(4.*np.pi*Plan.radial_slices[-1]*Plan.radial_slices[-1]*Plan.radial_slices[-1]))),\
                     fontsize=20)

    ax1 = plt.subplot2grid( (6,3) , (0,0), colspan=3, rowspan=3)
    ax2 = plt.subplot2grid( (6,3) , (3,0), colspan=3, rowspan=1)
    ax3 = plt.subplot2grid( (6,3) , (4,0), colspan=3, rowspan=1)
    ax4 = plt.subplot2grid( (6,3) , (5,0), colspan=3, rowspan=1)

    ax1.plot(Plan.radial_slices/1.e3,Plan.densities/1.e3,'k', linewidth = 2.)
    ax1.set_ylim(0.,(max(Plan.densities)/1.e3)+1.)
    ax1.set_xlim(0.,max(Plan.radial_slices)/1.e3)
    ax1.set_ylabel("Density ( $\cdot 10^3$ kg/m$^3$)")

    #Make a subplot showing the calculated pressure profile
    ax2.plot( Plan.radial_slices/1.e3, Plan.pressures/1.e9, 'b', linewidth=2.)
    ax2.set_ylim(0.,(max(Plan.pressures)/1e9)+10.)
    ax2.set_xlim(0.,max(Plan.radial_slices)/1.e3)
    ax2.set_ylabel("Pressure (GPa)")

    #Make a subplot showing the calculated gravity profile
    ax3.plot( Plan.radial_slices/1.e3, Plan.gravity, 'r', linewidth=2.)
    ax3.set_ylabel("Gravity (m/s$^2)$")
    ax3.set_xlim(0.,max(Plan.radial_slices)/1.e3)
    ax3.set_ylim(0.,max(Plan.gravity)+0.5)

    #Make a subplot showing the calculated temperature profile
    ax4.plot( Plan.radial_slices/1.e3, Plan.temperatures, 'g', linewidth=2.)
    ax4.set_ylabel("Temperature ($K$)")
    ax4.set_xlabel("Radius (km)")
    ax4.set_xlim(0.,max(Plan.radial_slices)/1.e3)
    ax4.set_ylim(0.,max(Plan.temperatures)+100)

    plt.show()
