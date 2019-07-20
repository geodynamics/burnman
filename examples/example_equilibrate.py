# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

'''
example_equilibrate
--------------------

This example demonstrates how burnman may be used to calculate the
equilibrium phase proportions and compositions for an assemblage
of a fixed bulk composition.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :func:`burnman.equilibrate.equilibrate`
'''
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')

import burnman
from burnman import equilibrate
from burnman.minerals import HP_2011_ds62, SLB_2011, JH_2015
from burnman.solutionbases import feasible_solution_in_component_space

# Examples:
#           ordering: calculates the state of order of the
#                     Jennings and Holland orthopyroxene
#                     in the simple en-fs binary at 1 bar.
#   aluminosilicates: plots the andalusite-sillimanite-kyanite phase diagram
#          gt_solvus: demonstrates the shape of the pyrope-grossular solvus
#       lower_mantle: calculates temperatures and assemblage properties
#                     along an isentrope in the lower mantle
#       upper_mantle: calculates a 2D grid of the ol-opx-gt field
# olivine_polymorphs: example produces a P-T pseudosection
#                     for a fo90 composition

fper_ol = False  # WORKS
upper_mantle = False  # WORKS (though not at highest T, 1500K)
ordering = False  # WORKS (better with decreasing temperatures)
aluminosilicates = False  # WORKS
gt_solvus = True  # WORKS
lower_mantle = False  # WORKS
olivine_polymorphs = False  # partially WORKS (last loop doesn't work)

gt = SLB_2011.garnet()
ol = SLB_2011.mg_fe_olivine()
wad = SLB_2011.mg_fe_wadsleyite()
rw = SLB_2011.mg_fe_ringwoodite()
bdg = SLB_2011.mg_fe_bridgmanite()
ppv = SLB_2011.post_perovskite()
fper = SLB_2011.ferropericlase()
opx = SLB_2011.orthopyroxene()
stv = SLB_2011.stishovite()
coe = SLB_2011.coesite()
cpv = SLB_2011.ca_perovskite()

ol.guess = np.array([0.93, 0.07])
wad.guess = np.array([0.91, 0.09])
rw.guess = np.array([0.93, 0.07])
opx.guess = np.array([0.8, 0.1, 0.05, 0.05])
gt.guess = np.array([0.8, 0.1, 0.05, 0.03, 0.02])
bdg.guess = np.array([0.86, 0.1, 0.04])
ppv.guess = np.array([0.86, 0.1, 0.04])
fper.guess = np.array([0.9, 0.1])

# ol fper

temperatures = np.linspace(1000., 1500., 11)
pressures = np.linspace(0.e9, 20.e9, 5)

if fper_ol:
    assemblage = burnman.Composite([ol, fper], [0.7, 0.3])
    for ph in assemblage.phases:
        if isinstance(ph, burnman.SolidSolution):
            ph.set_composition(ph.guess)

    assemblage.set_state(pressures[0], temperatures[0])
    equality_constraints = [('P', pressures), ('T', temperatures)]

    composition = {'Mg': 1., 'Fe': 0.5, 'Si': 0.5, 'O': 2.5}

    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False,
                           store_assemblage=True)

    for sols in sol:
        PGPa = sols[0].assemblage.pressure / 1.e9
        T = [s.assemblage.temperature for s in sols]
        x_fa = [s.assemblage.phases[0].molar_fractions[1] for s in sols]
        x_wus = [s.assemblage.phases[1].molar_fractions[1] for s in sols]
        p = plt.plot(T, x_fa,
                     label='p(fa), {0} GPa'.format(PGPa))
        plt.plot(T, x_wus, c=p[0].get_color(), linestyle='--',
                 label='p(wus), {0} GPa'.format(PGPa))
    plt.xlabel('Temperature (K)')
    plt.ylabel('proportion (mol fraction)')
    plt.legend()
    plt.show()


if upper_mantle:
    temperatures = np.linspace(800., 1500., 8)
    pressures = np.linspace(1.e9, 14.e9, 11)

    composition = {'Na': 0.02, 'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9,
                   'Ca': 0.2, 'Al': 0.4, 'O': 6.81}

    assemblage = burnman.Composite([ol, opx, gt], [0.7, 0.2, 0.1])
    for ph in assemblage.phases:
        if isinstance(ph, burnman.SolidSolution):
            ph.set_composition(ph.guess)

    assemblage.set_state(pressures[0], temperatures[0])
    equality_constraints = [('P', pressures), ('T', temperatures)]

    print('WARNING')
    composition = assemblage.formula

    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)

if ordering:
    orthopyroxene = feasible_solution_in_component_space(JH_2015.orthopyroxene(),
                                                         ['MgO', 'FeO', 'SiO2'],
                                                         solution_name='opx')

    orthopyroxene.guess = np.array([1./3., 1./3., 1./3.])
    temperatures = np.linspace(2000., 300., 41)

    Mg_numbers = np.linspace(10., 50., 5)
    assemblage = burnman.Composite([orthopyroxene])
    equality_constraints = [('P', 1.e5), ('T', temperatures)]

    for Mg_number in Mg_numbers:
        composition = {'Mg': Mg_number/100.*2.,
                       'Fe': (1.-Mg_number/100.)*2.,
                       'Si': 2., 'O': 6.}
        sols, prm = equilibrate(composition, assemblage,
                                equality_constraints, store_iterates=False)
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        p_fms = np.array([sol.x[-1] for sol in sols if sol.success])
        plt.plot(Ts, p_fms, label='Mg# = {0}'.format(Mg_number))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Proportion of ordered orthopyroxene")
    plt.legend(loc='best')
    plt.show()

if aluminosilicates:
    sillimanite = HP_2011_ds62.sill()
    andalusite = HP_2011_ds62.andalusite()
    kyanite = HP_2011_ds62.ky()

    composition = sillimanite.formula
    assemblage = burnman.Composite([sillimanite, andalusite, kyanite])
    equality_constraints = [('phase_proportion', (kyanite, np.array([0.0]))),
                            ('phase_proportion',
                             (sillimanite, np.array([0.0])))]

    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    P_inv, T_inv = sol.x[0:2]
    print('invariant point found at {0} GPa, {1} K'.format(P_inv/1e9, T_inv))

    low_pressures = np.linspace(1.e5, P_inv, 21)
    high_pressures = np.linspace(P_inv, 1.e9, 21)

    for pair, pressures in [([andalusite, kyanite], low_pressures),
                            ([sillimanite, andalusite], low_pressures),
                            ([sillimanite, kyanite], high_pressures)]:

        assemblage = burnman.Composite(pair)
        equality_constraints = [('P', pressures),
                                ('phase_proportion',
                                 (pair[0], np.array([0.0])))]
        sols, prm = equilibrate(composition, assemblage, equality_constraints,
                                store_iterates=False)
        Ps = np.array([sol.x[0] for sol in sols if sol.success])
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        plt.plot(Ts, Ps/1.e9)

    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')

    plt.show()

if gt_solvus:
    composition = {'Mg': 1.5, 'Ca': 1.5, 'Al': 2., 'Si': 3.0, 'O': 12.0}

    gt1 = feasible_solution_in_component_space(SLB_2011.garnet(),
                                               ['MgO', 'CaO', 'Al2Si3O9'],
                                               solution_name='gt1')

    gt2 = feasible_solution_in_component_space(SLB_2011.garnet(),
                                               ['MgO', 'CaO', 'Al2Si3O9'],
                                               solution_name='gt2')

    assemblage = burnman.Composite([gt1, gt2], [0.5, 0.5])
    gt1.set_composition([0.01, 0.99])
    gt2.set_composition([0.99, 0.01])

    pressure = 1.e5
    temperatures = np.linspace(100., 601.4, 41)
    equality_constraints = [('P', pressure), ('T', temperatures)]

    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False,
                            store_assemblage=True)

    x1s = np.array([sol.assemblage.phases[0].molar_fractions[1]
                    for sol in sols if sol.code == 0])
    x2s = np.array([sol.assemblage.phases[1].molar_fractions[1]
                    for sol in sols if sol.code == 0])
    Ts = np.array([sol.assemblage.temperature
                   for sol in sols if sol.code == 0])

    plt.plot(x1s, Ts)
    plt.plot(x2s, Ts)
    plt.xlabel('Molar proportion of pyrope')
    plt.ylabel('Temperature (K)')
    plt.show()


if lower_mantle:
    P0 = 25.e9
    T0 = 1600.

    composition = {'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9, 'Ca': 0.2, 'Al': 0.4,
                   'O': 6.8}
    assemblage = burnman.Composite([bdg, fper, cpv])
    equality_constraints = [('P', P0), ('T', T0)]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)

    S = np.array([assemblage.molar_entropy*assemblage.n_moles])

    assemblage = burnman.Composite([bdg, fper, ppv, cpv])
    assemblage.set_state(sol.x[0], sol.x[1])
    equality_constraints = [('S', S),
                            ('phase_proportion', (bdg, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    P_bdg_in = assemblage.pressure

    assemblage.set_state(sol.x[0], sol.x[1])
    equality_constraints = [('S', S),
                            ('phase_proportion', (ppv, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    P_ppv_in = assemblage.pressure
    T_ppv_in = assemblage.temperature

    pressures = np.linspace(P_ppv_in, P_bdg_in, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols1, prm1 = equilibrate(composition, assemblage, equality_constraints,
                              store_iterates=False)
    p1 = np.array([sol.x for sol in sols1 if sol.success]).T

    assemblage = burnman.Composite([bdg, fper, cpv])
    pressures = np.linspace(25.e9, P_ppv_in, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols2, prm2 = equilibrate(composition, assemblage,
                              equality_constraints,
                              store_iterates=False)

    p2 = np.array([sol.x for sol in sols2 if sol.success]).T

    assemblage = burnman.Composite([ppv, fper, cpv])
    pressures = np.linspace(P_bdg_in, 140.e9, 21)
    equality_constraints = [('P', pressures), ('S', S)]
    sols3, prm3 = equilibrate(composition, assemblage,
                              equality_constraints,
                              store_iterates=False)

    p3 = np.array([sol.x for sol in sols3 if sol.success]).T

    fig = plt.figure(figsize=(16, 8))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]
    phases = [bdg, fper, cpv, ppv]
    colors = ['red', 'green', 'blue', 'orange', 'purple']
    for p, prm in [(p1, prm1), (p2, prm2), (p3, prm3)]:

        pressures, temperatures = p[[0, 1], :]
        ax[0].plot(pressures/1.e9, temperatures, color='black')
        ax[0].set_xlabel('Pressure (GPa)')
        ax[0].set_ylabel('Temperature (K)')
        for i, phase in enumerate(phases):
            try:
                idx = prm.parameter_names.index('x({0})'.format(phase.name))
                x_phase = p[idx, :]
                ax[i+1].plot(pressures/1.e9, x_phase, color=colors[i])
                ax[i+1].set_ylim(0, 2)
                ax[i+1].set_xlim(0, 140)
                ax[i+1].set_xlabel('Pressure (GPa)')
                ax[i+1].set_ylabel('x({0})'.format(phase.name))
            except ValueError:
                pass

        try:
            pv_idx = prm.parameter_names.index('x({0})'.format(phases[0].name))
            per_idx = prm.parameter_names.index('x({0})'.format(phases[1].name))
            KD = (p[pv_idx+1, :] * (1. - p[per_idx+1, :])
                  / ((1. - p[pv_idx+1, :] - p[pv_idx+2, :])*p[per_idx+1, :]))
            ax[5].plot(pressures/1.e9, KD, color='red', label='pv KD')
        except:
            pass

        try:
            ppv_idx = prm.parameter_names.index('x({0})'.format(phases[3].name))
            per_idx = prm.parameter_names.index('x({0})'.format(phases[1].name))
            KD = (p[ppv_idx+1, :] * (1. - p[per_idx+1, :])
                  / ((1. - p[ppv_idx+1, :] - p[ppv_idx+2, :])*p[per_idx+1, :]))
            ax[5].plot(pressures/1.e9, KD, color='blue', label='ppv KD')
        except:
            pass

        ax[5].set_ylim(0., 1.)
        ax[5].set_xlabel('Pressure (GPa)')
        ax[5].set_ylabel('[FeSiO3/MgSiO3]/[FeO/MgO]')

    plt.show()


# This example produces a P-T pseudosection for a fo90 composition
if olivine_polymorphs:

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    T0 = 573.15
    T1 = 1773.15

    composition = {'Fe': 0.2, 'Mg': 1.8, 'Si': 1.0, 'O': 4.0}
    assemblage = burnman.Composite([ol, wad, rw])
    equality_constraints = [('phase_proportion', (ol, 0.0)),
                            ('phase_proportion', (rw, 0.0))]

    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    Pinv1, Tinv1 = sol.x[0:2]

    if sol.code != 0:
        raise Exception("Couldn't find ol-wad-rw invariant")

    assemblage = burnman.Composite([ol, wad, rw])
    equality_constraints = [('phase_proportion', (wad, 0.0)),
                            ('phase_proportion', (rw, 0.0))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           store_iterates=False)
    Pinv2, Tinv2 = sol.x[0:2]

    temperatures = np.linspace(T0, Tinv1, 8)
    assemblage = burnman.Composite([ol, wad, rw])

    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (ol, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(T0, Tinv2, 8)
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (wad, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(Tinv2, Tinv1, 8)
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(T0, Tinv2, 8)
    assemblage = burnman.Composite([ol, rw])
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(Tinv2, T1, 8)
    assemblage = burnman.Composite([ol, wad])
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (wad, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(Tinv1, T1, 8)
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (wad, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    temperatures = np.linspace(Tinv1, T1, 8)
    assemblage = burnman.Composite([wad, rw])
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    assemblage.set_fractions([0., 1.])
    temperatures = np.linspace(T0, T1, 8)
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    assemblage = burnman.Composite([rw, bdg, fper])
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 0.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    """
    equality_constraints = [('T', temperatures),
                            ('phase_proportion', (rw, 1.))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints,
                            store_iterates=False)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')
    """
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Pressure (GPa)')
    plt.show()
