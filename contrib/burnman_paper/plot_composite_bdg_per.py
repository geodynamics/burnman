import matplotlib.pyplot as plt
import numpy as np

import burnman_path  # adds the local burnman directory to the path
from burnman import minerals
from burnman import Composite
assert burnman_path  # silence pyflakes warning

bdg = minerals.SLB_2011.mg_fe_bridgmanite()
fper = minerals.SLB_2011.ferropericlase()

assemblage = Composite(phases=[bdg, fper],
                       fractions=[0.5, 0.5],
                       fraction_type='molar',
                       name='rock')

bdg = minerals.SLB_2011.mg_fe_bridgmanite()
fper = minerals.SLB_2011.ferropericlase()
assemblage = Composite(phases=[bdg, fper],
                       fractions=[0.5, 0.5],
                       fraction_type='molar',
                       name='mantle rock')

bdg.set_composition([0.9, 0.1, 0.0]) # MgSiO3, FeSiO3, AlAlO3
fper.set_composition([0.8, 0.2]) # MgO, FeO
assemblage.set_state(30.e9, 2000.)

print(f'Assemblage entropy: {assemblage.S:.1f} J/K/mol')


pressures = np.linspace(30.e9, 130.e9, 101)
temperatures = np.ones(101)*2000.

assemblage.set_averaging_scheme('VoigtReussHill')
densities, Vp_VRH, Vs_VRH = assemblage.evaluate(['rho', 'v_p', 'v_s'],
                                                pressures, temperatures)

assemblage.set_averaging_scheme('Reuss')
Vp_R, Vs_R = assemblage.evaluate(['v_p', 'v_s'], pressures, temperatures)

assemblage.set_averaging_scheme('Voigt')
Vp_V, Vs_V = assemblage.evaluate(['v_p', 'v_s'], pressures, temperatures)


fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]
ax[0].plot(pressures/1.e9, densities)
ax[1].fill_between(pressures/1.e9, Vp_R/1.e3, Vp_V/1.e3, alpha=0.2)
ax[1].fill_between(pressures/1.e9, Vs_R/1.e3, Vs_V/1.e3, alpha=0.2)
ax[1].plot(pressures/1.e9, Vs_R/1.e3, color='orange', linewidth=0.5)
ax[1].plot(pressures/1.e9, Vs_V/1.e3, color='orange', linewidth=0.5,
           label='$V_S$')
ax[1].plot(pressures/1.e9, Vp_R/1.e3, color='blue', linewidth=0.5)
ax[1].plot(pressures/1.e9, Vp_V/1.e3, color='blue', linewidth=0.5,
           label='$V_P$')

for i in range(2):
    ax[i].set_xlabel('Pressure (GPa)')

ax[0].set_ylabel('Density (kg/m$^3$)')
ax[1].set_ylabel('Velocities (km/s)')

ax[1].legend()
fig.tight_layout()
fig.savefig('figures/composite_bdg_per.pdf')
plt.show()
