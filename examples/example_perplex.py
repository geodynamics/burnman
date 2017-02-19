# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_perplex
---------------

This minimal example demonstrates how burnman can be used to read and interrogate
a PerpleX tab file (as produced by burnman/misc/create_burnman_readable_perplex_table.py 
It also demonstrates how we can smooth a given property on a given P-T grid.

*Uses:*

* :doc:`PerplexMaterial`


*Demonstrates:*

* Use of PerplexMaterial
* Smoothing gridded properties


"""
import sys
import os
import numpy as np

sys.path.insert(1, os.path.abspath('..'))

import burnman
import matplotlib.pyplot as plt


rock=burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')
P = 1.e9
T = 1650.
rock.set_state(P, T)
print('P: {0:.1f} GPa, T: {1:.1f} K, density: {2:.1f} kg/m^3'.format(P/1.e9, T, rock.rho))
                                                                    
pressures = np.linspace(1.e9, 100.e9, 101)
temperatures = [1650.] * len(pressures)
densities = rock.evaluate(['rho'], pressures, temperatures)[0]
plt.plot(pressures/1.e9, densities)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Density (kg/m^3)')
plt.show()


pressures = np.linspace(10.e9, 30.e9, 101)
temperatures = np.linspace(1600., 1700., 3)

pressure_stdev = 1.e9
temperature_stdev = 0.
pp, TT, entropies, smoothed_entropies = burnman.tools.smooth_gridded_property(rock, 'S', pressures, temperatures,
                                                                              pressure_stdev, temperature_stdev)


plt.plot(pp[0]/1.e9, entropies[0], label='entropies')
plt.plot(pp[0]/1.e9, smoothed_entropies[0], label='smoothed entropies')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Entropy (J/K/mol)')
plt.legend(loc='upper right')
plt.show()

