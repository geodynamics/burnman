# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
cubic_fitting
-------------
"""

from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from tools import print_table_for_mineral_constants

import burnman_path  # adds the local burnman directory to the path
import burnman

from burnman import AnisotropicMineral

assert burnman_path  # silence pyflakes warning

per = burnman.minerals.SLB_2011.periclase()
per.set_state(1.e5, 300.)
a = np.cbrt(per.params['V_0'])
cell_parameters = np.array([a, a, a, 90, 90, 90])

beta_RT = per.beta_T

per_data = np.loadtxt('data/isentropic_stiffness_tensor_periclase.dat')


def make_cubic_mineral_from_parameters(x):
    f_order = 3
    Pth_order = 1
    constants = np.zeros((6, 6, f_order+1, Pth_order+1))

    S11_0 = x[0]
    dS11df = x[1]
    d2S11df2 = x[2]
    dS11dPth = x[3] * 1.e-11
    d2S11dfdPth = x[4] * 1.e-11

    S44_0 = x[5]
    dS44df = x[6]
    d2S44df2 = x[7]
    dS44dPth = x[8] * 1.e-11
    d2S44dfdPth = x[9] * 1.e-11

    S12_0 = (1. - 3.*S11_0)/6.
    dS12df = -dS11df/2.
    d2S12df2 = -d2S11df2/2.
    dS12dPth = -dS11dPth/2.
    d2S12dfdPth = -d2S11dfdPth/2.

    constants[:3, :3, 1, 0] = S12_0
    constants[:3, :3, 2, 0] = dS12df
    constants[:3, :3, 3, 0] = d2S12df2
    constants[:3, :3, 0, 1] = dS12dPth
    constants[:3, :3, 1, 1] = d2S12dfdPth
    for i in range(3):
        constants[i, i, 1, 0] = S11_0
        constants[i, i, 2, 0] = dS11df
        constants[i, i, 3, 0] = d2S11df2

        constants[i, i, 0, 1] = dS11dPth
        constants[i, i, 1, 1] = d2S11dfdPth

    for i in range(3, 6):
        constants[i, i, 1, 0] = S44_0
        constants[i, i, 2, 0] = dS44df
        constants[i, i, 3, 0] = d2S44df2

        constants[i, i, 0, 1] = dS44dPth
        constants[i, i, 1, 1] = d2S44dfdPth

    return AnisotropicMineral(per, cell_parameters, constants)


parameters = []

run_fitting = False
if run_fitting:
    def cubic_misfit(x):
        m = make_cubic_mineral_from_parameters(x)

        chisqr = 0.
        for d in per_data:
            T, PGPa, Perr, rho, rhoerr = d[:5]
            C11S, C11Serr, C12S, C12Serr, C44S, C44Serr = d[5:]

            P = PGPa * 1.e9

            m.set_state(P, T)
            chisqr += np.power((m.isentropic_stiffness_tensor[0, 0]/1.e9
                                - C11S)/C11Serr, 2.)
            chisqr += np.power((m.isentropic_stiffness_tensor[0, 1]/1.e9
                                - C12S)/C12Serr, 2.)
            chisqr += np.power((m.isentropic_stiffness_tensor[3, 3]/1.e9
                                - C44S)/C44Serr, 2.)

        print(x)
        print(chisqr)

        return chisqr

    sol = minimize(cubic_misfit, [1./3., 0., 0., 0., 0., 1., 0., 0., 0., 0.],
                   method='COBYLA')

    parameters = sol.x

if not run_fitting:
    parameters = [0.64434719,  0.97982023,  2.28703418,  0.04069744,
                  0.83313498, 1.02999379, -3.39390829, -2.02738898,
                  0.06480835,  0.52939447]

m = make_cubic_mineral_from_parameters(parameters)

m.set_state(1.e5, 300.)

fig = plt.figure(figsize=(8, 8))
ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

pressures = np.linspace(1.e6, 30.e9, 101)
G_iso = np.empty_like(pressures)
G_aniso = np.empty_like(pressures)
C11 = np.empty_like(pressures)
C12 = np.empty_like(pressures)
C44 = np.empty_like(pressures)
G_slb = np.empty_like(pressures)
G_V = np.empty_like(pressures)
G_R = np.empty_like(pressures)

f = np.empty_like(pressures)
dXdf = np.empty_like(pressures)

temperatures = [300., 500., 700., 900.]
for T in temperatures:
    for i, P in enumerate(pressures):

        per.set_state(P, T)
        m.set_state(P, T)
        C11[i] = m.isentropic_stiffness_tensor[0,0]
        C12[i] = m.isentropic_stiffness_tensor[0,1]
        C44[i] = m.isentropic_stiffness_tensor[3,3]
        G_V[i] = m.isentropic_shear_modulus_voigt
        G_R[i] = m.isentropic_shear_modulus_reuss

        G_slb[i] = per.shear_modulus

    ax[0].plot(pressures/1.e9, C11/1.e9, label=f'{T} K')
    ax[1].plot(pressures/1.e9, C12/1.e9, label=f'{T} K')
    ax[2].plot(pressures/1.e9, C44/1.e9, label=f'{T} K')

    ln = ax[3].plot(pressures/1.e9, G_R/1.e9)
    ax[3].plot(pressures/1.e9, G_V/1.e9, color=ln[0].get_color(), label=f'{T} K')

    ax[3].fill_between(pressures/1.e9, G_R/1.e9, G_V/1.e9,
                       alpha = 0.3, color=ln[0].get_color())

    ax[3].plot(pressures/1.e9, G_slb/1.e9, label=f'{T} K (SLB2011)',
               linestyle='--', color=ln[0].get_color())

    # T, PGPa, Perr, rho, rhoerr, C11S, C11Serr, C12S, C12Serr, C44S, C44Serr
    T_data = np.array([[d[1], d[5], d[7], d[9]]
                       for d in per_data if np.abs(d[0] - T) < 1])

    ax[0].scatter(T_data[:, 0], T_data[:, 1])
    ax[1].scatter(T_data[:, 0], T_data[:, 2])
    ax[2].scatter(T_data[:, 0], T_data[:, 3])

for i in range(4):
    ax[i].set_xlabel('Pressure (GPa)')

ax[0].set_ylabel('$C_{N 11}$ (GPa)')
ax[1].set_ylabel('$C_{N 12}$ (GPa)')
ax[2].set_ylabel('$C_{N 44}$ (GPa)')
ax[3].set_ylabel('$G$ (GPa)')

for i in range(4):
    ax[i].legend()

fig.set_tight_layout(True)
fig.savefig('periclase_properties.pdf')
plt.show()

fig = plt.figure(figsize=(4, 4))
ax = [fig.add_subplot(1, 1, 1)]

temperatures = np.linspace(10., 2000., 101)
P = 1.e5
for i, T in enumerate(temperatures):
    m.set_state(P, T)
    per.set_state(P, T)
    C11[i] = m.isentropic_stiffness_tensor[0,0]
    C12[i] = m.isentropic_stiffness_tensor[0,1]
    C44[i] = m.isentropic_stiffness_tensor[3,3]
    G_V[i] = m.isentropic_shear_modulus_voigt
    G_R[i] = m.isentropic_shear_modulus_reuss

    G_slb[i] = per.shear_modulus
ax[0].plot(temperatures, C11/1.e9, label='$C_{N 11}$')
ax[0].plot(temperatures, C12/1.e9, label='$C_{N 12}$')
ax[0].plot(temperatures, C44/1.e9, label='$C_{44}$')
ln = ax[0].plot(temperatures, G_R/1.e9, color=ln[0].get_color(), label=f'$G$')
ax[0].plot(temperatures, G_V/1.e9, color=ln[0].get_color())
ax[0].fill_between(temperatures, G_R/1.e9, G_V/1.e9,
                   alpha = 0.3, color=ln[0].get_color(), zorder=2)
ax[0].plot(temperatures, G_slb/1.e9, label='$G$ (SLB2011)',
           color=ln[0].get_color(), linestyle='--', zorder=1)

# T, PGPa, Perr, rho, rhoerr, C11S, C11Serr, C12S, C12Serr, C44S, C44Serr
LP_data = np.array([[d[0], d[5], d[7], d[9]] for d in per_data if d[1] < 0.1])

ax[0].scatter(LP_data[:, 0], LP_data[:, 1])
ax[0].scatter(LP_data[:, 0], LP_data[:, 2])
ax[0].scatter(LP_data[:, 0], LP_data[:, 3])

ax[0].set_xlabel('Temperature (K)')
ax[0].set_ylabel('Elastic modulus (GPa)')

print_table_for_mineral_constants(m, [(1, 1), (1, 2), (4, 4)])

ax[0].legend()

fig.set_tight_layout(True)
fig.savefig('periclase_properties_1bar.pdf')
plt.show()
