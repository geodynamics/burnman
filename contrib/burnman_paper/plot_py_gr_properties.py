import matplotlib.pyplot as plt
import numpy as np


import burnman_path  # adds the local burnman directory to the path
from burnman import minerals
from burnman.tools.chemistry import formula_to_string
assert burnman_path  # silence pyflakes warning

gt = minerals.JH_2015.garnet()
print(gt.endmember_names)

gt.set_composition([0.4, 0., 0.6, 0., 0.])
gt.set_state(5.e9, 1000.)

np.set_printoptions(precision=0)
print(f'Formula: {formula_to_string(gt.formula)}')
print(f'Molar volume: {gt.V*1.e6:.2f} cm^3/mol')
print(f'Partial entropies: {gt.partial_entropies}')



comp = np.linspace(1e-5, 1.-1e-5, 101)
temperatures = [600., 1000., 1400.]
colors = ['r', 'g', 'b']
gt_excess_gibbs = np.empty((3, 101))
gt_activities = np.empty((3, 101, 5))

pressure = 1.e9
for i, c in enumerate(comp):
    molar_fractions = [1.0 - c, 0., c, 0., 0.]
    gt.set_composition(molar_fractions)
    for j, temperature in enumerate(temperatures):
        gt.set_state(pressure, temperature)
        gt_excess_gibbs[j][i] = gt.excess_gibbs
        gt_activities[j][i] = gt.activities

fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

for i, temperature in enumerate(temperatures):
    ax[0].plot(comp, gt_excess_gibbs[i]/1000., f'{colors[i]}-',
               linewidth=1., label=f'{int(temperature)} K')

    ax[1].plot(comp, gt_activities[i, :, 0], f'{colors[i]}-',
               linewidth=1., label=f'pyrope ({int(temperature)} K)')
    ax[1].plot(comp, gt_activities[i, :, 2], f'{colors[i]}:',
               linewidth=1., label=f'grossular ({int(temperature)} K)')


ax[0].set_ylabel("Excess Gibbs energy of solution (kJ/mol)")
ax[1].set_ylabel("Endmember activities")

for i in range(2):
    ax[i].set_xlabel("Molar grossular fraction")
    ax[i].legend(loc='lower left')

fig.tight_layout()
fig.savefig('figures/mg_ca_gt_properties.pdf')
plt.show()
