from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))
import numpy as np
import matplotlib.pyplot as plt
import burnman

import scipy.integrate
import time


def old_thermal(T, debye_T, n):
    if T == 0:
        return 0
    return 3. * n * constants.R * T * debye_fn(debye_T / T)


def old_heat(T, debye_T, n):
    if T == 0:
        return 0
    deb = scipy.integrate.quad(
        lambda x: pow(x, 4.) * np.exp(x) / pow((np.exp(x) - 1.), 2.), 0.0, debye_T / T)
    return 9. * n * burnman.constants.gas_constant * deb[0] / pow(debye_T / T, 3.)

temperatures = np.linspace(100, 5000, 10000)
Debye_T = 1000.
old = np.empty_like(temperatures)
start = time.clock()
for i in range(len(temperatures)):
    old[i] = old_heat(temperatures[i], Debye_T, 1.0)
time_old = time.clock() - start

new = np.empty_like(temperatures)
start = time.clock()
for i in range(len(temperatures)):
    new[i] = burnman.eos.debye.heat_capacity_v(temperatures[i], Debye_T, 1.0)
time_new = time.clock() - start

print("error %e" % np.linalg.norm((old - new) / new))
print("time old %g, time new %g" % (time_old, time_new))


temperatures = np.linspace(0, 5000, 200)
vibrational_energy = np.empty_like(temperatures)
heat_capacity = np.empty_like(temperatures)
Debye_T = 1000.
for i in range(len(temperatures)):
    vibrational_energy[i] = burnman.eos.debye.thermal_energy(
        temperatures[i], Debye_T, 1.0)
    heat_capacity[i] = burnman.eos.debye.heat_capacity_v(
        temperatures[i], Debye_T, 1.0)

plt.subplot(121)
plt.plot(temperatures, vibrational_energy)
plt.subplot(122)
plt.plot(temperatures, heat_capacity)
plt.show()
