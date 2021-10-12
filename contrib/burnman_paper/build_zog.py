from __future__ import absolute_import
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np


import burnman_path  # adds the local burnman directory to the path
import burnman
from burnman import minerals
from burnman import Layer, Planet, Mineral, PerplexMaterial
from burnman import BoundaryLayerPerturbation
import warnings
import urllib.request

assert burnman_path  # silence pyflakes warning

icb_radius = 1220.e3
inner_core = Layer('inner core', radii=np.linspace(0., icb_radius, 21))
Fe_Dewaele = minerals.other.Fe_Dewaele()
params = Fe_Dewaele.params

hcp_iron = minerals.SE_2015.hcp_iron()
params = hcp_iron.params

params['name'] = 'modified solid iron'
params['formula'] = {'Fe': 1.0}
inner_core_material = Mineral(params=params,
                              property_modifiers=[['linear',
                                                   {'delta_E': 0.,
                                                    'delta_S': 0.,
                                                    'delta_V': 2.45e-7}]])

inner_core.set_material(inner_core_material)
#inner_core.set_temperature_mode('user-defined',
#                                np.nan*np.ones_like(inner_core.radii))

inner_core.set_temperature_mode('adiabatic')

cmb_radius = 3480.e3
outer_core = Layer('outer core', radii=np.linspace(icb_radius, cmb_radius, 21))
EPOC = minerals.ICL_2018.EPOC_vinet()
params = EPOC.params
params['name'] = 'modified EPOC'


liq_iron = minerals.SE_2015.liquid_iron()
params = liq_iron.params

params['name'] = 'modified liquid iron'
params['formula'] = {'Fe': 1.0}

outer_core_material = Mineral(params=params,
                              property_modifiers=[['linear',
                                                   {'delta_E': 0.,
                                                    'delta_S': 0.,
                                                    'delta_V': 2.28e-7}]])
outer_core.set_material(outer_core_material)
#outer_core.set_temperature_mode('user-defined',
#                                np.nan*np.ones_like(outer_core.radii))

outer_core.set_temperature_mode('adiabatic')

urllib.request.urlretrieve("https://raw.githubusercontent.com/"
                           "bobmyhill/bobmyhill.github.io/master/"
                           "files/pyrolite_perplex_table.dat",
                           "pyrolite_perplex_table.dat")


lab_radius = 6001.e3
moho_radius = 6341.e3
lab_temperature = 1645. - 119.5  # to be consistent with TBL, defined below
moho_temperature = 620.
surface_temperature = 300.

convecting_mantle_radii = np.linspace(cmb_radius, lab_radius, 101)
pyrolite = PerplexMaterial('./pyrolite_perplex_table.dat', name='pyrolite')
convecting_mantle = Layer('convecting mantle',
                          radii=convecting_mantle_radii)
convecting_mantle.set_material(pyrolite)

# Here we add a thermal boundary layer perturbation, assuming that the
# lower mantle has a Rayleigh number of 1.e7, and that there
# is an 840 K jump across the basal thermal boundary layer and an
# 84 K jump at the top of the lower mantle.
tbl_perturbation = BoundaryLayerPerturbation(radius_bottom=cmb_radius,
                                             radius_top=lab_radius,
                                             rayleigh_number=1.e7,
                                             temperature_change=0.,  # tweaked
                                             boundary_layer_ratio=0.)
dTdr_top = ((moho_temperature - lab_temperature)
            / (moho_radius - lab_radius))
tbl_perturbation.set_model_thermal_gradients(-19.e-3, dTdr_top)

lower_mantle_tbl = tbl_perturbation.temperature(convecting_mantle_radii)

print(lower_mantle_tbl[0], lower_mantle_tbl[-1])
convecting_mantle.set_temperature_mode('perturbed-adiabatic',
                                       temperatures=lower_mantle_tbl)

dunite = minerals.SLB_2011.mg_fe_olivine(molar_fractions=[0.92, 0.08])
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


bounds = np.array([[layer.radii[0]/1.e3, layer.radii[-1]/1.e3]
                   for layer in planet_zog.layers])
maxy = [15, 400, 12, 5000]
for bound in bounds:
    for i in range(4):
        ax[i].fill_betweenx([0., maxy[i]],
                            [bound[0], bound[0]],
                            [bound[1], bound[1]], alpha=0.2)

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
ax[3].set_ylabel('Temperature (K)')
ax[3].set_xlabel('Radius (km)')
ax[3].set_ylim(0.,)

for i in range(2):
    ax[i].set_xticklabels([])
for i in range(4):
    ax[i].set_xlim(0., max(planet_zog.radii) / 1.e3)
    ax[i].set_ylim(0., maxy[i])

fig.tight_layout()
fig.savefig('figures/zog.pdf')
plt.show()
