# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_perplex
---------------

This minimal example demonstrates how burnman can be used to read and interrogate
a PerpleX tab file (as produced by burnman/misc/create_burnman_readable_perplex_table.py 

*Uses:*

* :doc:`PerplexMaterial`



*Demonstrates:*

* Use of PerplexMaterial


"""
import sys
import os
import numpy as np

sys.path.insert(1, os.path.abspath('..'))

import burnman
import matplotlib.pyplot as plt


rock=burnman.PerplexMaterial('../burnman/data/input_perplex/in23_1.tab')
rock.set_state(1.e9, 300.)
pressures = np.linspace(1.e9, 100.e9, 101)
temperatures = [1650.] * len(pressures)
densities = rock.evaluate(['rho'], pressures, temperatures)[0]
plt.plot(pressures/1.e9, densities)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Density (kg/m^3)')
plt.show()

