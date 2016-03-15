# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.
# Ian Rose ian.r.rose@gmail.com

"""
CIDER 2014 BurnMan Tutorial --- step 1
--------------------------------------

In this first part of the tutorial we will acquaint ourselves with a basic
script for calculating the elastic properties of a mantle mineralogical
model.

In general, there are three portions of this script:

1) Define a set of pressures and temperatures at which we want to
calculate elastic properties

2) Setup a composite of minerals (or "rock") and calculate its
elastic properties at those pressures and temperatures.

3) Plot those elastic properties, and compare them to a seismic
model, in this case PREM


The script is basically already written, and should run as is by typing:

    python step_1.py

on the command line.  However, the mineral model for the rock is not
very realistic, and you will want to change it to one that is more
in accordance with what we think the bulk composition of Earth's lower mantle is.

"""
from __future__ import absolute_import


# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
import os
import sys
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../../'))

# Here we import the relevant modules from BurnMan.  The burnman
# module imports several of the most important functionalities of
# the library, including the ability to make composites, and compute
# thermoelastic properties of them.  The minerals module includes
# the mineral physical parameters for the predefined minerals in
# BurnMan
import burnman
from burnman import minerals


if __name__ == '__main__':
    """
    Part (1) --- Defining the pressures and temperatures.
    We get the pressures from the PREM density model for a list of
    depths in Earth, and a model for the temperature profile with
    depth from Brown and Shankland (1981).  In addition, we get the
    density and elastic properties from PREM, which we will use later
    for comparison purposes
    """

    # Here we create and load the PREM seismic velocity model, which will be
    # used for comparison with the seismic velocities of the "rock" composite
    seismic_model = burnman.seismic.PREM()

    # We create an array of 20 depths at which we want to evaluate PREM, and then
    # query the seismic model for the pressure, density, P wave speed, S wave
    # speed, and bulk sound velocity at those depths
    n_depths = 20
    min_depth = 850.e3
    max_depth = 2800.e3
    depths = np.linspace(min_depth, max_depth, n_depths)
    pressure, seis_rho, seis_vp, seis_vs = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s'], depths)

    # Now we get an array of temperatures at which will be used for computing
    # the seismic properties of the rock.  Here we use the Brown+Shankland (1981)
    # geotherm for mapping pressure to temperature
    temperature = burnman.geotherm.brown_shankland(pressure)

    """
    Part (2) -- Defining the rock

    The object burnman.Composite expects two lists, one with the molar
    fractions of the different minerals making up the rock, and one with
    the minerals themselves.

    Here we setup a rock made from stishovite (SiO_2) and wustite (FeO), which is
    probably not the most realistic mineralogical model for Earth's lower
    mantle. This is the place where you will want to make changes.

    For our simplified mantle model, consider a three element model, with
    Mg, Si, and O.  This will make a mantle with magnesium perovskite (MgSiO_3)
    and periclase (MgO). We can replace minerals.SLB_2011.stishovite() and
    minerals.SLB_2011.wuestite() with minerals.SLB_2011.mg_perovskite()
    and minerals.SLB_2011.periclase() and play with the relative fraction
    of the two phases.

    """

    # ---------------------------------------------------------#
    # ------------- MAKE MODIFICATIONS HERE -------------------#
    # ---------------------------------------------------------#

    phase_1_fraction = 0.5
    phase_2_fraction = 1.0 - phase_1_fraction
    rock = burnman.Composite(
        [minerals.SLB_2011.stishovite(), minerals.SLB_2011.wuestite()], [phase_1_fraction, phase_2_fraction])

    # ---------------------------------------------------------#
    # ---------------------------------------------------------#
    # ---------------------------------------------------------#

    # At this point we want to tell the rock which equation of state to use for
    # its thermoelastic calculations. In general, we recommend the 'slb3'
    # equation of state as the most self-consistent model.  The parameters from
    # the SLB_2011 mineral library are fit using this model.
    rock.set_method('slb3')

    # Here is the step which does the heavy lifting.  burnman.velocities_from_rock
    # sets the state of the rock at each of the pressures and temperatures defined,
    # then calculates the elastic moduli and density of each individual phase.  After that,
    # it performs elastic averaging on the phases to get a single bulk and shear
    # modulus for the rock.  This averaging scheme defaults to Voigt-Reuss-Hilli,
    # but see example_averaging.py for other options.  Finally, it calculates the seismic
    # wave speeds for the whole rock.  It returns a tuple of density, p-wave velocity
    # s-wave velocity, bulk sound speed, bulk modulus, and shear modulus.
    density, vp, vs = rock.evaluate(
        ['density', 'v_p', 'v_s'], pressure, temperature)

    """
    Part (3) --- Plotting and comparison

    This should not need to be changed, though of course you may if you like.  It plots the
    resulting densities, shear wave speeds and bulk sound speeds against pressure, as
    well as plotting the relevant values from PREM for comparison.

    You can try modifying the phase_1_fraction above to try to get a closer fit, or
    a more complicated mineralogical model.

    You will probably have a hard time getting a close match with our simple model.
    In particular, the densities will be too low and the wave speeds will be too fast.
    This is likely due to the fact that we have not included iron, which will make the
    rocks denser, which will in turn slow down the wave speeds.

    """

    # All the work is done except the plotting!  Here we want to plot the seismic wave
    # speeds and the density against PREM using the matplotlib plotting tools.  We make
    # a 2x2 array of plots.  The fourth subplot plots the geotherm used for
    # this calculation.

    # First, we plot the s-wave speed verses the PREM s-wave speed
    plt.subplot(2, 2, 1)
    plt.plot(pressure / 1.e9, vs / 1.e3, color='b', linestyle='-',
             marker='o', markerfacecolor='b', markersize=4, label='computation')
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4, label='PREM')
    plt.title("S wave speed (km/s)")
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.legend(loc='lower right')
    plt.ylim(5, 8.0)

    # Next, we plot the p-wave speed verses the PREM p-wave speed
    plt.subplot(2, 2, 2)
    plt.plot(pressure / 1.e9, vp / 1.e3, color='b', linestyle='-',
             marker='o', markerfacecolor='b', markersize=4)
    plt.plot(pressure / 1.e9, seis_vp / 1.e3, color='k',
             linestyle='-', marker='o', markerfacecolor='k', markersize=4)
    plt.title("P wave speed (km/s)")
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.ylim(10, 16)

    # Next, we plot the density verses the PREM density
    plt.subplot(2, 2, 3)
    plt.plot(pressure / 1.e9, density / 1.e3, color='b',
             linestyle='-', marker='o', markerfacecolor='b', markersize=4)
    plt.plot(pressure / 1.e9, seis_rho / 1.e3, color='k',
             linestyle='-', marker='o', markerfacecolor='k', markersize=4)
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.xlabel("Pressure (GPa)")
    plt.title("density (kg/m^3)")

    # Finally, we plot the goetherm used
    plt.subplot(2, 2, 4)
    plt.plot(pressure / 1e9, temperature, color='r', linestyle='-',
             marker='o', markerfacecolor='r', markersize=4)
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.xlabel("Pressure (GPa)")
    plt.title("temperature (K)")

    # At long last, we show the results!  We are done!
    plt.show()
