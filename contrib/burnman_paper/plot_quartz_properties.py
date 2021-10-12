import matplotlib.pyplot as plt
import numpy as np


import burnman_path  # adds the local burnman directory to the path
from burnman import minerals

assert burnman_path  # silence pyflakes warning

# The calls to burnman
qtz = minerals.SLB_2011.quartz()
temperatures = np.linspace(300., 1300., 1001)
pressures = 1.e5*np.ones_like(temperatures)
isobaric_heat_capacities, bulk_sound_velocities = qtz.evaluate(['molar_heat_capacity_p', 'v_phi'], pressures, temperatures)

# The following lines do the plotting
fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

ax[0].plot(temperatures, isobaric_heat_capacities)
ax[1].plot(temperatures, bulk_sound_velocities/1000.)

for i in range(2):
  ax[i].set_xlabel('Temperature (K)')
ax[0].set_ylabel('Isobaric heat capacity (J/K/mol)')
ax[1].set_ylabel('Bulk sound velocity (km/s)')
ax[0].set_ylim(50., 150.)
ax[1].set_ylim(1.5, 5.)
fig.tight_layout()
fig.savefig('figures/quartz_properties.pdf')
plt.show()
