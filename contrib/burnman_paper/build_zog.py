from __future__ import absolute_import
from __future__ import print_function

from burnman import minerals
from burnman import Layer, Planet, Mineral, PerplexMaterial

import matplotlib.pyplot as plt
import numpy as np

from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline

import burnman_path  # adds the local burnman directory to the path
import burnman

import warnings

assert burnman_path  # silence pyflakes warning

import urllib.request

from burnman import Layer, Planet, Mineral, PerplexMaterial
from burnman import minerals

icb_radius = 1220.e3
inner_core = Layer('inner core', radii=np.linspace(0., icb_radius, 21))
Fe_Dewaele = minerals.other.Fe_Dewaele()
params = Fe_Dewaele.params
params['name'] = 'modified solid iron'
params['formula'] = {'Fe': 1.0}
inner_core_material = Mineral(params=params,
                              property_modifiers=[['linear',
                                                   {'delta_E': 0.,
                                                    'delta_S': 0.,
                                                    'delta_V': 3.95e-7}]])

inner_core.set_material(inner_core_material)
inner_core.set_temperature_mode('user-defined',
                                np.nan*np.ones_like(inner_core.radii))


cmb_radius = 3480.e3
outer_core = Layer('outer core', radii=np.linspace(icb_radius, cmb_radius, 21))
EPOC = minerals.ICL_2018.EPOC_vinet()
params = EPOC.params
params['name'] = 'modified EPOC'
outer_core_material = Mineral(params=params,
                              property_modifiers=[['linear',
                                                   {'delta_E': 0.,
                                                    'delta_S': 0.,
                                                    'delta_V': 4.e-8}]])
outer_core.set_material(outer_core_material)
outer_core.set_temperature_mode('user-defined',
                                np.nan*np.ones_like(outer_core.radii))

urllib.request.urlretrieve("https://raw.githubusercontent.com/"
                           "bobmyhill/bobmyhill.github.io/master/"
                           "files/pyrolite_perplex_table.dat",
                           "pyrolite_perplex_table.dat")


lab_radius = 6171.e3
pyrolite = PerplexMaterial('./pyrolite_perplex_table.dat', name='pyrolite')
convecting_mantle = Layer('convecting mantle',
                          radii=np.linspace(cmb_radius, lab_radius, 101))
convecting_mantle.set_material(pyrolite)
convecting_mantle.set_temperature_mode('adiabatic')


moho_radius = 6341.e3
lab_temperature = 1600.
moho_temperature = 620.
surface_temperature = 300.

dunite = minerals.SLB_2011.mg_fe_olivine(molar_fractions=[0.915, 0.085])
lithospheric_mantle = Layer('lithospheric mantle',
                            radii=np.linspace(lab_radius, moho_radius, 31))
lithospheric_mantle.set_material(dunite)
lithospheric_mantle.set_temperature_mode('user-defined',
                                         np.linspace(lab_temperature,
                                                     moho_temperature, 31))


planet_radius = 6371.e3
andesine = minerals.SLB_2011.plagioclase(molar_fractions=[0.4, 0.6])
crust = Layer('crust', radii=np.linspace(moho_radius, planet_radius, 11))
crust.set_material(andesine)
crust.set_temperature_mode('user-defined',
                           np.linspace(moho_temperature,
                                       surface_temperature, 11))

planet_zog = Planet('Planet Zog',
                    [inner_core, outer_core,
                     convecting_mantle, lithospheric_mantle,
                     crust], verbose=True)
planet_zog.make()


earth_mass = 5.972e24
earth_moment_of_inertia_factor = 0.3307

print(f'mass = {planet_zog.mass:.3e} (Earth = {earth_mass:.3e})')
print(f'moment of inertia factor= {planet_zog.moment_of_inertia_factor:.4f} '
      f'(Earth = {earth_moment_of_inertia_factor:.4f})')

print('Layer mass fractions:')
for layer in planet_zog.layers:
    print(f'{layer.name}: {layer.mass / planet_zog.mass:.3f}')


# Let's get PREM to compare everything to as we are trying
# to imitate Earth
prem = burnman.seismic.PREM()
premradii = 6371.e3 - prem.internal_depth_list()

with warnings.catch_warnings(record=True) as w:
    eval = prem.evaluate(['density', 'pressure', 'gravity', 'v_s', 'v_p'])
    premdensity, prempressure, premgravity, premvs, premvp = eval
    print(w[-1].message)



fig = plt.figure(figsize=(8, 5))
ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

ax[0].plot(planet_zog.radii / 1.e3, planet_zog.density / 1.e3,
           label=planet_zog.name)
ax[0].plot(premradii / 1.e3, premdensity / 1.e3, linestyle=':', label='PREM')
ax[0].set_ylabel('Density ($10^3$ kg/m$^3$)')
ax[0].legend()

# Make a subplot showing the calculated pressure profile
ax[1].plot(planet_zog.radii / 1.e3, planet_zog.pressure / 1.e9)
ax[1].plot(premradii / 1.e3, prempressure / 1.e9, linestyle=':')
ax[1].set_ylabel('Pressure (GPa)')

# Make a subplot showing the calculated gravity profile
ax[2].plot(planet_zog.radii / 1.e3, planet_zog.gravity)
ax[2].plot(premradii / 1.e3, premgravity, linestyle=':')
ax[2].set_ylabel('Gravity (m/s$^2)$')
ax[2].set_xlabel('Radius (km)')

# Make a subplot showing the calculated temperature profile
ax[3].plot(planet_zog.radii / 1.e3, planet_zog.temperature)
ax[3].set_ylabel('Temperature ($K$)')
ax[3].set_xlabel('Radius (km)')
ax[3].set_ylim(0.,)

for i in range(2):
    ax[i].set_xticklabels([])
for i in range(4):
    ax[i].set_xlim(0., max(planet_zog.radii) / 1.e3)

fig.tight_layout()
plt.show()
