from __future__ import absolute_import
from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

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
                                                    'delta_V': 1.95e-7}]])

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
                                                    'delta_V': 1.78e-7}]])
outer_core.set_material(outer_core_material)
#outer_core.set_temperature_mode('user-defined',
#                                np.nan*np.ones_like(outer_core.radii))

outer_core.set_temperature_mode('adiabatic')

try:
    pyrolite = PerplexMaterial('./pyrolite_perplex_table.dat', name='pyrolite')
except FileNotFoundError:
    urllib.request.urlretrieve("https://raw.githubusercontent.com/"
                               "bobmyhill/bobmyhill.github.io/master/"
                               "files/pyrolite_perplex_table.dat",
                               "pyrolite_perplex_table.dat")
    pyrolite = PerplexMaterial('./pyrolite_perplex_table.dat', name='pyrolite')

lm_radius = 5701.e3
lab_radius = 6171.e3 # 200 km thick lithosphere

lower_mantle_radii = np.linspace(cmb_radius, lm_radius, 101)
lower_mantle = Layer('lower mantle',
                     radii=lower_mantle_radii)
lower_mantle.set_material(pyrolite)

# Here we add a thermal boundary layer perturbation, assuming that the
# lower mantle has a Rayleigh number of 1.e7, and that the basal thermal
# boundary layer has a thermal gradient equal to dTdr_bottom
dTdr_bottom = -33.e-3 # this figure roughly matches the Anzellini geotherm.
dTdr_bottom = -19.e-3
dTdr_top = -2.e-3 # a small TBL at the UM/LM boundary
dTdr_sa = -300. / (cmb_radius - lab_radius) # following Anderson (1982), with a somewhat higher gradient
dT_superadiabatic = -dTdr_sa * (lower_mantle_radii - lower_mantle_radii[-1])

tbl_perturbation = BoundaryLayerPerturbation(radius_bottom=cmb_radius,
                                             radius_top=lm_radius,
                                             rayleigh_number=1.e7,
                                             temperature_change=0.,  # tweaked
                                             boundary_layer_ratio=0.)

tbl_perturbation.set_model_thermal_gradients(dTdr_bottom-dTdr_sa, dTdr_top-dTdr_sa)

lower_mantle_tbl = (tbl_perturbation.temperature(lower_mantle_radii)
                    + dT_superadiabatic)

lower_mantle.set_temperature_mode('perturbed-adiabatic',
                                  temperatures=lower_mantle_tbl)

moho_radius = 6341.e3
lab_temperature = 1645. - 49.  # to be consistent with TBL, defined below
moho_temperature = 620.
surface_temperature = 300.

upper_convecting_mantle_radii = np.linspace(lm_radius, lab_radius, 31)
upper_convecting_mantle = Layer('upper convecting mantle',
                                radii=upper_convecting_mantle_radii)
upper_convecting_mantle.set_material(pyrolite)

# Here we add a thermal boundary layer perturbation, assuming that the
# lower mantle has a Rayleigh number of 1.e7, and that the basal thermal
# boundary layer has a thermal gradient equal to dTdr_bottom
dTdr_bottom = dTdr_top # gradient must be the same as the layer above
dTdr_top = ((moho_temperature - lab_temperature)
            / (moho_radius - lab_radius))

dT_superadiabatic = -dTdr_sa * (upper_convecting_mantle_radii - upper_convecting_mantle_radii[-1])

tbl_perturbation = BoundaryLayerPerturbation(radius_bottom=lm_radius,
                                             radius_top=lab_radius,
                                             rayleigh_number=1.e7,
                                             temperature_change=0.,  # tweaked
                                             boundary_layer_ratio=0.)

tbl_perturbation.set_model_thermal_gradients(dTdr_bottom-dTdr_sa, dTdr_top-dTdr_sa)

upper_convecting_mantle_tbl = (tbl_perturbation.temperature(upper_convecting_mantle_radii)
                               + dT_superadiabatic)

print(upper_convecting_mantle_tbl[0], upper_convecting_mantle_tbl[-1])
upper_convecting_mantle.set_temperature_mode('perturbed-adiabatic',
                                             temperatures=upper_convecting_mantle_tbl)

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
                     lower_mantle, upper_convecting_mantle,
                     lithospheric_mantle, crust],
                    verbose=True)
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
maxy = [15, 400, 12, 7000]
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

labels = ['Stacey (1977)',
          'Brown and Shankland (1981)',
          'Anderson (1982)',
          'Alfe et al. (2007)',
          'Anzellini et al. (2013)']

labels = ['S1977',
          'BS1981',
          'A1982',
          'A2007',
          'A2013']


ax[3].plot(planet_zog.radii / 1.e3,
burnman.geotherm.stacey_continental(planet_zog.depth),
linestyle=':', label=labels[0])
mask = planet_zog.depth > 269999.
ax[3].plot(planet_zog.radii[mask] / 1.e3,
           burnman.geotherm.brown_shankland(planet_zog.depth[mask]),
           linestyle=':', label=labels[1])
ax[3].plot(planet_zog.radii / 1.e3,
           burnman.geotherm.anderson(planet_zog.depth),
           linestyle=':', label=labels[2])

ax[3].scatter([planet_zog.layers[0].radii[-1] / 1.e3,
               planet_zog.layers[1].radii[-1] / 1.e3],
              [5400., 4000.],
              linestyle=':', label=labels[3])

d = np.loadtxt('Anzellini_2013_geotherm.dat')
Anz_interp = interp1d(d[:,0]*1.e9, d[:,1])
mask = planet_zog.pressure < 330.e9
temperatures = Anz_interp(planet_zog.pressure[mask])
ax[3].plot(planet_zog.radii[mask] / 1.e3, temperatures,
           linestyle=':', label=labels[4])

ax[3].legend()

for i in range(2):
    ax[i].set_xticklabels([])
for i in range(4):
    ax[i].set_xlim(0., max(planet_zog.radii) / 1.e3)
    ax[i].set_ylim(0., maxy[i])

fig.tight_layout()
fig.savefig('figures/zog.pdf')
plt.show()
