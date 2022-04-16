# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
cubic_fitting
-------------

This script creates an AnisotropicMineral object corresponding to
periclase (a cubic mineral). If run_fitting is set to True, the script uses
experimental data to find the optimal anisotropic parameters.
If set to False, it uses pre-optimized parameters.
The data is used only to optimize the anisotropic parameters;
the isotropic equation of state is taken from
Stixrude and Lithgow-Bertelloni (2011).

The script ends by making three plots; one with elastic moduli
at high pressure, one with the corresponding shear moduli,
and one with the elastic moduli at 1 bar.
"""

from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from tools import print_table_for_mineral_constants

import burnman

from burnman import AnisotropicMineral


run_fitting = False

per = burnman.minerals.SLB_2011.periclase()
V0 = per.params['V_0']
K0 = per.params['K_0']
Kprime0 = per.params['Kprime_0']

per.set_state(1.e5, 300.)
beta_RT = per.beta_T


def psi_func(f, Pth, params):
    dPsidf = params['a'] + params['b_1']*params['c_1']*np.exp(params['c_1']*f) + params['b_2']*params['c_2']*np.exp(params['c_2']*f)
    Psi = 0. + params['a']*f + params['b_1']*np.exp(params['c_1']*f) + params['b_2']*np.exp(params['c_2']*f)
    dPsidPth = 0.*params['b_1']
    return (Psi, dPsidf, dPsidPth)


def make_cubic_mineral_from_parameters(x):

    anisotropic_parameters = {'a': np.zeros((6, 6)),
                              'b_1': np.zeros((6, 6)),
                              'c_1': np.ones((6, 6)),
                              'b_2': np.zeros((6, 6)),
                              'c_2': np.ones((6, 6))}

    per.params['V_0'] = x[0]*V0
    per.params['K_0'] = x[1]*K0
    per.params['Kprime_0'] = x[2]*Kprime0

    a = np.cbrt(per.params['V_0'])
    cell_parameters = np.array([a, a, a, 90, 90, 90])

    i = 3
    for p in ['a', 'b_1', 'c_1', 'b_2', 'c_2']:
        if p == 'a':
            anisotropic_parameters[p][:3,:3] = (1. - 3.*x[i])/6.
        elif p[0] == 'b' or p == 'd':
            anisotropic_parameters[p][:3,:3] = -3.*x[i]/6.
        else:
            anisotropic_parameters[p][:3,:3] = x[i]

        anisotropic_parameters[p][0, 0] = x[i]
        anisotropic_parameters[p][1, 1] = x[i]
        anisotropic_parameters[p][2, 2] = x[i]
        i = i + 1
        anisotropic_parameters[p][3, 3] = x[i]
        anisotropic_parameters[p][4, 4] = x[i]
        anisotropic_parameters[p][5, 5] = x[i]
        i = i + 1
    assert len(x) == i

    return AnisotropicMineral(per, cell_parameters, anisotropic_parameters, psi_func, orthotropic=True)


per_data = np.loadtxt('data/Karki_2000_periclase_CSijs.dat')
parameters = []

if run_fitting:
    def cubic_misfit(x, imin):
        m = make_cubic_mineral_from_parameters(x)

        chisqr = 0.
        try:
            for d in per_data[::1]:
                T, P, C11S, C12S, C44S, KT, V, betaTminusbetaS = d

                C11Serr = 0.01*C11S
                C12Serr = 0.01*C12S
                C44Serr = 0.01*C44S
                Verr = 0.01*V

                if T < 500.:
                    m.set_state(P, T)

                    chisqr += np.power((m.V - V)/Verr, 2.)
                    chisqr += np.power((m.isentropic_stiffness_tensor[0, 0]
                                        - C11S)/C11Serr, 2.)
                    chisqr += np.power((m.isentropic_stiffness_tensor[0, 1]
                                        - C12S)/C12Serr, 2.)
                    chisqr += np.power((m.isentropic_stiffness_tensor[3, 3]
                                        - C44S)/C44Serr, 2.)

            if np.isnan(chisqr):
                print(d, "Noooo, there was a nan")
                chisqr = 1.e7

        except Exception as e:
            print(e)
            print('There was an exception')
            chisqr = 1.e7
        imin[0][0] += 1
        if chisqr < imin[0][1]:
            imin[0][1] = chisqr
            print(imin[0])
            print(repr(x))

        return chisqr

    i = 0
    min = 1.e10
    guesses = np.array([1., 1., 1.,
                        1./3., 1.,
                        0., 0.,
                        1., 1.,
                        0., 0.,
                        10., 10.])

    guesses = np.array([ 1.00762096,  0.99900194,  1.08547729,  0.47944268,  0.04445275,
       -0.77294538, -0.17422384,  5.81473511,  1.13934342,  0.79772947,
       -0.53002167,  5.81567299, -2.37105322])

    #sol = minimize(cubic_misfit, guesses, method='COBYLA', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 10000})
    #sol = minimize(cubic_misfit, guesses, method='Powell', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 10000}) # good for rapid convergence
    sol = minimize(cubic_misfit, guesses, method='Nelder-Mead', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 100000}) # fine tuning
    print(sol)
    parameters = sol.x

if not run_fitting:

    parameters = np.array([ 1.00762096,  0.99900194,  1.08547729,  0.47944268,  0.04445275,
       -0.77294538, -0.17422384,  5.81473511,  1.13934342,  0.79772947,
       -0.53002167,  5.81567299, -2.37105322])


m = make_cubic_mineral_from_parameters(parameters)

assert burnman.tools.eos.check_anisotropic_eos_consistency(m)

print(repr(m.anisotropic_params['a']))
print(repr(m.anisotropic_params['b_1']))
print(repr(m.anisotropic_params['c_1']))
print(repr(m.anisotropic_params['b_2']))
print(repr(m.anisotropic_params['c_2']))
exit()

m.set_state(1.e5, 300.)

fig = plt.figure(figsize=(4, 12))
fig2 = plt.figure(figsize=(4, 4))
ax = [fig.add_subplot(3, 1, i) for i in range(1, 4)]
ax.append(fig2.add_subplot(1, 1, 1))

pressures = np.linspace(1.e6, 150.e9, 101)
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

temperatures = [300., 1000., 2000.] #, 3000.]
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

    # T, P, C11S, C12S, C44S, KT, V, betaTminusbetaS
    T_data = np.array([[d[1], d[2], d[3], d[4]]
                       for d in per_data if np.abs(d[0] - T) < 1])/1.e9

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
fig.savefig('periclase_stiffness_tensor.pdf')
fig2.set_tight_layout(True)
fig2.savefig('periclase_shear_modulus.pdf')
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
