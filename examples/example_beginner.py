# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_beginner
----------------

This example script is intended for absolute beginners to BurnMan.
We cover importing BurnMan modules, creating a composite material,
and calculating its seismic properties at lower mantle pressures and
temperatures.  Afterwards, we plot it against a 1D seismic model
for visual comparison.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.material.Material.evaluate`


*Demonstrates:*

* creating basic composites
* calculating thermoelastic properties
* seismic comparison
"""
from __future__ import absolute_import

# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


# Here we import the relevant modules from BurnMan.  The burnman
# module imports several of the most important functionalities of
# the library, including the ability to make composites, and compute
# thermoelastic properties of them.  The minerals module includes
# the mineral physical parameters for the predefined minerals in
# BurnMan
import burnman
from burnman import minerals


if __name__ == "__main__":

    # This is the first actual work done in this example.  We define
    # composite object and name it "rock".  A composite is made by
    # giving burnman.composite a list of minerals and their molar fractions.
    # Here "rock" has two constituent minerals: it is 80% Mg perovskite
    # and 20% periclase.  More minerals may be added by simply extending
    # the list given to burnman.composite
    # For the preset minerals from the SLB_2011, the equation of state
    # formulation from Stixrude and Lithgow-Bertolloni (2005) will be used.
    rock = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                              minerals.SLB_2011.periclase()],
                             [0.8, 0.2])

    # Here we create and load the PREM seismic velocity model, which will be
    # used for comparison with the seismic velocities of the "rock" composite
    seismic_model = burnman.seismic.PREM()

    # We create an array of 20 depths at which we want to evaluate PREM, and then
    # query the seismic model for the pressure, density, P wave speed, S wave
    # speed, and bulk sound velocity at those depths
    depths = np.linspace(750e3, 2800e3, 20)
    pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    # Now we get an array of temperatures at which will be used for computing
    # the seismic properties of the rock.  Here we use the Brown+Shankland (1981)
    # geotherm for mapping pressure to temperature
    temperature = burnman.geotherm.brown_shankland(pressure)

    # Here is the step which does the heavy lifting.  burnman.evaluate
    # sets the state of the rock at each of the pressures and temperatures defined,
    # then calculates the elastic moduli and density of each individual phase.  After that,
    # it calcalates all the variables asked for. For the composite a default averaging scheme
    # of Voigt-Reuss-Hill is used (this can be changes with set_averaging
    # scheme)
    density, vp, vs, vphi = rock.evaluate(
        ['density', 'v_p', 'v_s', 'v_phi'], pressure, temperature)

    # All the work is done except the plotting!  Here we want to plot the seismic wave
    # speeds and the density against PREM using the matplotlib plotting tools.  We make
    # a 2x2 array of plots.  The fourth subplot plots the geotherm used for
    # this calculation.
    # First, we plot the s-wave speed verses the PREM s-wave speed
    plt.subplot(2, 2, 1)
    plt.plot(pressure / 1.e9, vs / 1.e3, color='b', linestyle='-',
             marker='o', markerfacecolor='b', markersize=4, label='computation')
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4, label='reference')
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
