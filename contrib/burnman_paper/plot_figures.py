import burnman
from burnman import minerals
import numpy as np
import matplotlib.pyplot as plt

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


import numpy as np
import matplotlib.pyplot as plt

comp = np.linspace(1e-5, 1.-1e-5, 101)

bdg_excess_gibbs_400 = np.empty_like(comp)
bdg_excess_gibbs_800 = np.empty_like(comp)
bdg_excess_gibbs_1200 = np.empty_like(comp)

bdg_activities_400 = np.empty((101, 3))
bdg_activities_800 = np.empty((101, 3))
bdg_activities_1200 = np.empty((101, 3))

pressure = 30.e9
for i, c in enumerate(comp):
    molar_fractions = [1.0 - c, c, 0.]
    bdg.set_composition(molar_fractions)
    bdg.set_state(pressure, 400.)
    bdg_excess_gibbs_400[i] = bdg.excess_gibbs
    bdg_activities_400[i] = bdg.activities
    bdg.set_state(pressure, 800.)
    bdg_excess_gibbs_800[i] = bdg.excess_gibbs
    bdg_activities_800[i] = bdg.activities
    bdg.set_state(pressure, 1200.)
    bdg_excess_gibbs_1200[i] = bdg.excess_gibbs
    bdg_activities_1200[i] = bdg.activities

fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]
ax[0].plot(comp, bdg_excess_gibbs_400/1000., 'r-', linewidth=1., label='400 K')
ax[0].plot(comp, bdg_excess_gibbs_800/1000., 'g-', linewidth=1., label='800 K')
ax[0].plot(comp, bdg_excess_gibbs_1200/1000., 'b-', linewidth=1., label='1200 K')

ax[1].plot(comp, bdg_activities_1200[:,0], 'b-', linewidth=1., label='MgSiO$_3$')
ax[1].plot(comp, bdg_activities_1200[:,1], 'b:', linewidth=1., label='FeSiO$_3$')


ax[0].set_ylabel("Excess Gibbs energy of solution (kJ/mol)")
ax[1].set_ylabel("Endmember activities")

for i in range(2):
  ax[i].set_xlabel("Molar FeSiO$_3$ fraction")
  ax[i].legend(loc='lower left')

fig.tight_layout()
fig.savefig('figures/fe_mg_bdg_properties.pdf')
plt.show()



from burnman import Composite

bdg = minerals.SLB_2011.mg_fe_bridgmanite()
fper = minerals.SLB_2011.ferropericlase()

assemblage = Composite(phases=[bdg, fper],
                       fractions=[0.5, 0.5],
                       fraction_type='molar',
                       name='rock')

bdg.set_composition([0.9, 0.1, 0.0])
fper.set_composition([0.8, 0.2])

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
ax[1].plot(pressures/1.e9, Vs_V/1.e3, color='orange', linewidth=0.5, label='$V_S$')
ax[1].plot(pressures/1.e9, Vp_R/1.e3, color='blue', linewidth=0.5)
ax[1].plot(pressures/1.e9, Vp_V/1.e3, color='blue', linewidth=0.5, label='$V_P$')

for i in range(2):
  ax[i].set_xlabel('Pressure (GPa)')

ax[0].set_ylabel('Density (kg/m$3$)')
ax[1].set_ylabel('Velocities (km/s)')

ax[1].legend()
fig.tight_layout()
fig.savefig('figures/composite_bdg_per.pdf')
plt.show()
