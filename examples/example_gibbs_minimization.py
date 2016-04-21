# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

'''
example_gibbs_minimization
--------------------

Quite often in the Earth Sciences, we are interested in finding the
minimum gibbs free energy of a composite material.
Here we demonstrate how burnman can be used to find the compositions
and phase fractions of a composite which minimize the free energy.


*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :func:`burnman.equilibriumassemblage.gibbs_minimizer`
'''
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman.minerals import HP_2011_ds62, SLB_2011
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from burnman.equilibriumassemblage import gibbs_minimizer


# Solve 1: Here's the classic aluminosilicate phase diagram
andalusite = HP_2011_ds62.andalusite()
sillimanite = HP_2011_ds62.sill()
kyanite = HP_2011_ds62.ky()

composition = andalusite.params['formula']

assemblage = burnman.Composite([andalusite, sillimanite, kyanite])
constraints= [['X', [1., 0., 0.], [1., 1., 1.], 0.],
              ['X', [0., 1., 0.], [1., 1., 1.], 0.]]
sol = gibbs_minimizer(composition, assemblage, constraints)
P_inv = sol[0]
T_inv = sol[1]

# Kyanite-sillimanite appears at positive P
lo_pressures = np.linspace(1.e5, P_inv, 21)
hi_pressures = np.linspace(P_inv, 2.*P_inv, 21)

and_ky_temperatures = np.empty_like(lo_pressures)
and_sill_temperatures = np.empty_like(lo_pressures)
sill_ky_temperatures = np.empty_like(hi_pressures)
 
for i, P in enumerate(lo_pressures):
    constraints= [['P', P], ['X', [1., 0.], [1., 1.], 0.]]
    
    assemblage = burnman.Composite([andalusite, kyanite])
    and_ky_temperatures[i] = gibbs_minimizer(composition, assemblage, constraints)[1]
    
    assemblage = burnman.Composite([andalusite, sillimanite])
    and_sill_temperatures[i] = gibbs_minimizer(composition, assemblage, constraints)[1]

for i, P in enumerate(hi_pressures):
    constraints= [['P', P], ['X', [1., 0.], [1., 1.], 0.]]
    
    assemblage = burnman.Composite([sillimanite, kyanite])
    sill_ky_temperatures[i] = gibbs_minimizer(composition, assemblage, constraints)[1]

plt.plot(and_ky_temperatures, lo_pressures/1.e9, label='and-ky')
plt.plot(and_sill_temperatures, lo_pressures/1.e9, label='and-sill')
plt.plot(sill_ky_temperatures, hi_pressures/1.e9, label='sill-ky')


plt.title('Aluminosilicate phase diagram')
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (GPa)')
plt.legend(loc='upper left')
plt.show()




# Solve 2: Let's plot the olivine-wadsleyite phase loop at 1400 K
ol = SLB_2011.mg_fe_olivine()
wad = SLB_2011.mg_fe_wadsleyite()
rw = SLB_2011.mg_fe_ringwoodite()

def delta_invariant_temperature(arg, T):
    x_ol = arg[0]
    assemblage = burnman.Composite([ol, wad, rw])
    composition = { 'Mg': 2.*(1. - x_ol), 'Fe': 2.*x_ol, 'O': 4., 'Si': 1.}
    constraints= [['X', [0., 0., 1., 0., 0., 0.], [1., 0., 1., 0., 1., 0.], 0.],
                  ['X', [0., 0., 0., 0., 1., 0.], [1., 0., 1., 0., 1., 0.], 0.]]
    inv = gibbs_minimizer(composition, assemblage, constraints)
    return inv[1] - T

T = 1400.
opt.fsolve(delta_invariant_temperature, [0.25], T)

p_inv = ol.pressure
x_ol_inv = ol.molar_fractions[1]
x_wad_inv = wad.molar_fractions[1]
x_rw_inv = rw.molar_fractions[1]

plt.plot([x_ol_inv, x_rw_inv], [p_inv/1.e9, p_inv/1.e9])

assemblages = [burnman.Composite([ol, wad]),
               burnman.Composite([ol, rw]),
               burnman.Composite([wad, rw])]

c_bounds = [[0.000001, x_ol_inv, 21], [x_ol_inv, 0.999999, 21], [0.0000001, x_wad_inv, 21]]

for i, assemblage in enumerate(assemblages):
    x1 = np.linspace(*c_bounds[i])
    x2 = np.empty_like(x1)
    pressures = np.empty_like(x1)
    
    for i, x in enumerate(x1):
        composition = { 'Mg': 2.*(1. - x), 'Fe': 2.*x, 'O': 4., 'Si': 1.}
        constraints= [['T', T], ['X', [1., 0., 0., 0.], [1., 0., 1., 0.], 0.999999]]
        sol = gibbs_minimizer(composition, assemblage, constraints)
        pressures[i] = sol[0]
        x2[i] = sol[5]
    plt.plot(x1, pressures/1.e9, label=assemblage.phases[1].name+'-in')
    plt.plot(x2, pressures/1.e9, label=assemblage.phases[0].name+'-in')
    
plt.title('Mg2SiO4-Fe2SiO4 phase diagram at '+str(T)+' K')
plt.xlabel('p (Fe2SiO4)')
plt.ylabel('Pressure (GPa)')
plt.legend(loc='upper right')
plt.show()



