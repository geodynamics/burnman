import matplotlib.pyplot as plt
import numpy as np

import burnman_path  # adds the local burnman directory to the path
from burnman import minerals

assert burnman_path  # silence pyflakes warning


# Bridgmanite properties
bdg = minerals.SLB_2011.mg_fe_bridgmanite()
bdg.set_composition([0.8, 0.1, 0.1])
bdg.set_state(30.e9, 2000.)

pressures = np.linspace(30.e9, 130.e9, 101)
temperatures = 2000.*np.ones_like(pressures)
densities, p_wave_velocities = bdg.evaluate(['rho', 'v_p'], pressures, temperatures)

# The following lines do the plotting
fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

ax[0].plot(pressures/1.e9, densities)
ax[1].plot(pressures/1.e9, p_wave_velocities/1000.)

for i in range(2):
    ax[i].set_xlabel('Pressure (GPa)')
ax[0].set_ylabel('Densities (kg/m^3)')
ax[1].set_ylabel('P wave velocity (km/s)')
fig.tight_layout()

fig.savefig('figures/bdg_80_10_10_properties.pdf')
plt.show()
