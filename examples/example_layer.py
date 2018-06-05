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
* calculating thermoelastic properties with self-consistent pressures
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
    depths = np.linspace(2890e3, 670e3, 20)
    pressure, gravity, seis_rho, seis_vp, seis_vs, seis_bullen = seismic_model.evaluate(
        ['pressure', 'gravity', 'density', 'v_p', 'v_s', 'bullen'], depths)

    # Here we define the lower mantle as a Layer(). The layer needs various
    # parameters to set a depth array and radius array.
    lower_mantle = burnman.Layer( name='Lower Mantle', radii=6371.e3-depths)
    # Here we set the composition of the layer as the above defined 'rock'.
    lower_mantle.set_material(rock)

    # Now we set the temperature mode of the layer.
    # Here we use an adiabatic temperature and set the temperature at the top of the layer
    lower_mantle.set_temperature_mode(temperature_mode='adiabatic', temperature_top=1900.)
    #Alternatively, we choose a user-defined temperature, given by the Brown & Shankland geotherm
    #lower_mantle.set_temperature_mode(temperature_mode ='user_defined',
    #                                      temperatures =burnman.geotherm.brown_shankland(depths))

    # And we set a self-consistent pressure. The pressure at the top of the layer and
    # gravity at the bottom of the layer are given by the PREM.
    lower_mantle.set_pressure_mode(pressure_mode='self-consistent',
                                   pressure_top=pressure[-1], gravity_bottom=gravity[0])
    # Alternatively, we set a user-defined pressure given by PREM
    #lower_mantle.set_pressure_mode(pressure_mode='user-defined',
    #                               pressures = pressure, gravity_bottom=gravity[0])

    lower_mantle.make()
    
    # All the work is done, now we can plot various properties!
    fig = plt.figure(figsize = (10,6))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]
    
    # First, we plot the p-wave speed verses the PREM p-wave speed
    ax[0].plot(lower_mantle.pressures / 1.e9, lower_mantle.v_p / 1.e3, color='r', linestyle='-',
               marker='o', markerfacecolor='r', markersize=4, label='V$_P$ (computed)')
    ax[0].plot(pressure / 1.e9, seis_vp / 1.e3, color='r',
               linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='V$_P$ (reference)')

    # Next, we plot the s-wave speed verses the PREM s-wave speed
    ax[0].plot(lower_mantle.pressures / 1.e9, lower_mantle.v_s / 1.e3, color='b', linestyle='-',
             marker='o', markerfacecolor='b', markersize=4, label='V$_S$ (computed)')
    ax[0].plot(pressure / 1.e9, seis_vs / 1.e3, color='b', linestyle='--',
             marker='o', markerfacecolor='w', markersize=4, label='V$_S$ (reference)')

    ax[0].set_ylabel("Wave speeds (km/s)")
    
    # Next, we plot the density versus the PREM density
    ax[1].plot(lower_mantle.pressures / 1.e9, lower_mantle.density / 1.e3, color='g',
             linestyle='-', marker='o', markerfacecolor='g', markersize=4, label='computed')
    ax[1].plot(pressure / 1.e9, seis_rho / 1.e3, color='g',
             linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='reference')
    ax[1].set_ylabel("Density (g/cm$^3$)")

    # And the Bullen parameter
    ax[2].plot(lower_mantle.pressures / 1e9, lower_mantle.bullen, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4, label='computed')
    ax[2].plot(pressure / 1.e9, seis_bullen , color='k',
                 linestyle='--', marker='o', markerfacecolor='w', markersize=4, label='reference')
    ax[2].set_ylabel("Bullen parameter")

    # Finally, we plot the used geotherm
    ax[3].plot(lower_mantle.pressures / 1e9, lower_mantle.temperatures, color='k', linestyle='-',
             marker='o', markerfacecolor='k', markersize=4, label='used geotherm')
    ax[3].set_ylabel("Temperature (K)")

    for i in range(4):
        ax[i].set_xlabel("Pressure (GPa)")
        ax[i].set_xlim(min(pressure) / 1.e9, max(pressure) / 1.e9)
        ax[i].legend(loc='best')
        
    # At long last, we show the results!  We are done!
    plt.show()
