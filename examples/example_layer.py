# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_layer
----------------

This example script is building on example_beginner.py.
In this case we use the layer Class. Instead of computing properties at
pressures defined by the seismic PREM model (as is done in many examples), 
we compute properties at pressures self-consistent to the layer.
Layer can be used to evaluate geophysical parameter, such as
the Bullen parameter or the Brunt_vasala frequency. Through the 'modified_adiabat'
temperature setting it allows for inclusions of thermal boundary layers. 
Layers can also be used to build an entire planet (see example_build_planet.py)

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.layer.Layer`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`


*Demonstrates:*

* creating basic layer
* calculating thermoelastic properties with selfconsistent pressures
* seismic comparison
"""
from __future__ import absolute_import

# Here we import standard python modules that are required for
# usage of BurnMan.
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))


# Here we import the relevant modules from BurnMan.
import burnman
from burnman import minerals


if __name__ == "__main__":

    # This is the first actual work done in this example.  We define
    # composite object and name it "rock".
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
    pressure, gravity, seis_rho, seis_vp, seis_vs, seis_bullen = seismic_model.evaluate(
        ['pressure', 'gravity', 'density', 'v_p', 'v_s', 'bullen'], depths)

    # Here we define the lower mantle as a Layer(). The layer needs various
    # parameters to set a depth array and radius array.
    lower_mantle = burnman.Layer( name='Lower Mantle', radius_planet=6371.e3,
        min_depth=750.e3, max_depth=2800.e3,  n_slices=20)
    # Here we set the composition of the layer as the above defined 'rock'.
    lower_mantle.set_composition(rock)

    # Now we set the temperature mode of the layer. , and a
    # self-consistent pressure. The pressure at the top of the layer and
    # gravity at the bottom of the layer are given by the PREM.
    
    # Here we use an adiabatic temperature and set the temperature at the top of the layer
    lower_mantle.set_temperature_mode(temperature_mode='adiabat', temperature_top=1900.)
    #Alternatively, we choose a user-defined temperature, given by the Brown & Shankland geotherm
    #lower_mantle.set_temperature_mode(temperature_mode ='user_defined',
                                          temperatures =burnman.geotherm.brown_shankland(depths))
                           
    lower_mantle.set_state(pressure_mode='selfconsistent',
                               pressure_top=pressure[0], gravity_bottom=gravity[-1])

       
    # All the work is done, now we can plot various properties!
    # First, we plot the s-wave speed verses the PREM s-wave speed
    plt.figure(figsize = (10,6))
    plt.subplot(1, 2, 1)
    plt.plot(lower_mantle.pressure / 1.e9, lower_mantle.v_s / 1.e3, color='b', linestyle='-',
             marker='o', markerfacecolor='b', markersize=4, label='Vs computed')
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, color='b', linestyle='--',
             marker='o', markerfacecolor='w', markersize=4, label='Vs ref')
    plt.ylabel("Wave speeds (km/s)  and density (g/cm^3)")
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.legend(loc='lower right')

    # Next, we plot the p-wave speed verses the PREM p-wave speed
    plt.plot(lower_mantle.pressure / 1.e9, lower_mantle.v_p / 1.e3, color='r', linestyle='-',
             marker='o', markerfacecolor='r', markersize=4, label='Vp computed')
    plt.plot(pressure / 1.e9, seis_vp / 1.e3, color='r',
             linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='Vp ref')


    # Next, we plot the density verses the PREM density
    plt.plot(lower_mantle.pressure / 1.e9, lower_mantle.density / 1.e3, color='g',
             linestyle='-', marker='o', markerfacecolor='g', markersize=4, label='rho computed')
    plt.plot(pressure / 1.e9, seis_rho / 1.e3, color='g',
             linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='rho ref')
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.xlabel("Pressure (GPa)")
    plt.xlim(min(lower_mantle.pressure) / 1.e9, max(lower_mantle.pressure) / 1.e9)
    
    plt.legend()
    
    plt.subplot(2, 2, 2)
    plt.plot(lower_mantle.pressure / 1e9, lower_mantle.bullen, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4, label='bullen computed')
    plt.plot(pressure / 1.e9, seis_bullen , color='k',
                 linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='bullen ref')
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.legend()
    plt.ylabel("Bullen parameter")

    # Finally, we plot the used geotherm
    plt.subplot(2, 2, 4)
    plt.plot(lower_mantle.pressure / 1e9, lower_mantle.temperature, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4)
    plt.xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("temperature (K)")

    # At long last, we show the results!  We are done!
    plt.show()
