# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

'''
example_equilibrate
--------------------

This example demonstrates how BurnMan may be used to calculate the
equilibrium state for an assemblage of a fixed bulk composition given
two constraints. Each constraint has the form
[<constraint type>, <constraint>], where
<constraint type> is one of the strings: P, T, S, V, X, PT_ellipse,
phase_fraction, or phase_composition. The <constraint> object should
either be a float or an array of floats for P, T, S, V
(representing the desired pressure, temperature,
entropy or volume of the material). If the constraint type is X
(a generic constraint on the solution vector) then the constraint c is
represented by the following equality:
np.dot(c[0], x) - c[1]. If the constraint type is
PT_ellipse, the equality is given by
norm(([P, T] - c[0])/c[1]) - 1.
The constraint_type phase_fraction assumes a tuple of the phase object
(which must be one of the phases in the burnman.Composite) and a float
or vector corresponding to the phase fractions. Finally, a
phase_composition constraint has the format (site_names, n, d, v),
where site names dictates the sites involved in the equality constraint.
The equality constraint is given by n*x/d*x = v,
where x are the site occupancies and n and d are fixed vectors of
site coefficients. So, one could for example choose a constraint
([Mg_A, Fe_A], [1., 0.], [1., 1.], [0.5]) which would
correspond to equal amounts Mg and Fe on the A site.

This script provides a number of examples, which can be turned on and off
with a series of boolean variables. In order of complexity:

* run_aluminosilicates: Creates the classic aluminosilicate diagram involving
  univariate reactions between andalusite, sillimanite and kyanite.
* run_ordering: Calculates the state of order of Jennings and Holland (2015)
  orthopyroxene in the simple en-fs binary at 1 bar.
* run_gt_solvus: Demonstrates the shape of the pyrope-grossular solvus.
* run_fper_ol: Calculates the equilibrium Mg-Fe partitioning between
  ferropericlase and olivine.
* run_fixed_ol_composition: Calculates the composition of wadsleyite in
  equilibrium with olivine of a fixed composition at a fixed pressure.
* run_upper_mantle: Calculates the equilibrium compositions and
  phase proportions for an ol-opx-gt composite in an NCFMAS bulk composition.
* run_lower_mantle: Calculates temperatures and assemblage properties
  along an isentrope in the lower mantle. Includes calculations of the
  post-perovskite-in and bridgmanite-out lines.
* run_olivine_polymorphs: Produces a P-T pseudosection for a fo90 composition.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Composite`
* :func:`burnman.equilibrate.equilibrate`
'''
from __future__ import absolute_import
from __future__ import print_function

from copy import copy
import numpy as np
import matplotlib.pyplot as plt

import burnman_path
import burnman
from burnman import equilibrate
from burnman.minerals import HP_2011_ds62, SLB_2011, JH_2015
from burnman.tools.polytope import simplify_composite_with_composition

assert burnman_path  # silence pyflakes warning

run_aluminosilicates = True
run_ordering = True
run_gt_solvus = True
run_fper_ol = True
run_fixed_ol_composition = True
run_upper_mantle = True
run_lower_mantle = True
run_olivine_polymorphs = True

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

if __name__ == "__main__" and run_aluminosilicates:
    """
    Creates the classic aluminosilicate diagram involving
    univariate reactions between andalusite, sillimanite and kyanite.
    """
    sillimanite = HP_2011_ds62.sill()
    andalusite = HP_2011_ds62.andalusite()
    kyanite = HP_2011_ds62.ky()

    # First, find the pressure and temperature of the invariant point
    composition = sillimanite.formula
    assemblage = burnman.Composite([sillimanite, andalusite, kyanite])
    equality_constraints = [('phase_fraction', (kyanite, np.array([0.0]))),
                            ('phase_fraction', (sillimanite, np.array([0.0])))]

    sol, prm = equilibrate(composition, assemblage, equality_constraints)
    P_inv, T_inv = sol.x[0:2]
    print(f'invariant point found at {P_inv/1e9:.2f} GPa, {T_inv:.2f} K')

    # Now we can find the univariant lines which all converge on the
    # invariant point. In this case, we assume we know which side of each line
    # is stable (because the aluminosilicate diagram is so ubiquitous in
    # metamorphic textbooks), but if we didn't know,
    # we could also calculate which field is stable around the invariant point
    # by checking to see which had the minimum Gibbs energy.
    low_pressures = np.linspace(1.e5, P_inv, 21)
    high_pressures = np.linspace(P_inv, 1.e9, 21)

    for pair, pressures in [([andalusite, kyanite], low_pressures),
                            ([sillimanite, andalusite], low_pressures),
                            ([sillimanite, kyanite], high_pressures)]:

        assemblage = burnman.Composite(pair)
        equality_constraints = [('P', pressures),
                                ('phase_fraction',
                                 (pair[0], np.array([0.0])))]
        sols, prm = equilibrate(composition, assemblage,
                                equality_constraints)
        Ps = np.array([sol.x[0] for sol in sols if sol.success])
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        plt.plot(Ts, Ps/1.e9)

    # Now we plot each stable field. Here, we do explicitly check
    # which polymorph is stable at each point.
    assemblage = burnman.Composite([sillimanite, andalusite, kyanite])
    for (P, T) in [[0.2e9, 800.],
                   [0.6e9, 650],
                   [0.5e9, 1000]]:
        assemblage.set_state(P, T)

        stable_phase = kyanite
        min_gibbs = kyanite.gibbs

        for phase in [sillimanite, andalusite]:
            if phase.gibbs < min_gibbs:
                min_gibbs = phase.gibbs
                stable_phase = phase

        plt.text(T, P/1.e9, stable_phase.name)

    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')

    plt.show()

if __name__ == "__main__" and run_ordering:
    """
    Calculates the state of order of Jennings and Holland (2015)
    orthopyroxene in the simple en-fs binary at 1 bar.
    """
    assemblage = burnman.Composite([JH_2015.orthopyroxene()])

    composition = {'Mg': 50/100.*2.,
                   'Fe': (1.-50/100.)*2.,
                   'Si': 2., 'O': 6.}

    # The JH2015 orthopyroxene has some components which do not
    # exist in our simple bulk composition, so here we call a
    # helper function to simplify any phases that cannot be represented
    # by the bulk composition. This can either remove endmembers, or in some
    # cases redefine the independent set of endmembers to create the simplest
    # solutions possible.
    assemblage = simplify_composite_with_composition(
        assemblage, composition)

    assemblage.phases[0].set_composition([1./3., 1./3., 1./3.])

    temperatures = np.linspace(2000., 300., 41)

    Mg_numbers = np.linspace(10., 50., 5)
    equality_constraints = [('P', 1.e5), ('T', temperatures)]

    # Loop over Mg numbers
    for Mg_number in Mg_numbers:
        composition = {'Mg': Mg_number/100.*2.,
                       'Fe': (1.-Mg_number/100.)*2.,
                       'Si': 2., 'O': 6.}
        sols, prm = equilibrate(composition, assemblage,
                                equality_constraints)
        # We know that we've made a simpler orthopyroxene solution.
        # In this case, the endmembers are the same, and the names
        # have just been appended with "in child solution".
        # We get the index of the ordered endmember here.
        idx = prm.parameter_names.index('p(ordered ferroenstatite '
                                        'in child solution)')

        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        p_fms = np.array([sol.x[idx] for sol in sols if sol.success])
        plt.plot(Ts, p_fms, label=f'Mg# = {Mg_number}')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Proportion of ordered orthopyroxene")
    plt.legend(loc='best')
    plt.show()

if __name__ == "__main__" and run_gt_solvus:
    """
    Demonstrates the shape of the pyrope-grossular solvus.
    """
    composition = {'Mg': 1.5, 'Ca': 1.5, 'Al': 2., 'Si': 3.0, 'O': 12.0}

    # The assemblage in this case is made of two distinct copies of the
    # SLB_2011 garnet solution.
    # In the bulk composition of interest, there are only
    # two stable endmembers of the solution. Here, we simplify each phase
    # so that it only contains those endmembers
    a = burnman.Composite([SLB_2011.garnet()])
    a = simplify_composite_with_composition(a, composition)
    py_gr = a.phases[0]

    assemblage = burnman.Composite([py_gr, copy(py_gr)],
                                   [0.5, 0.5])

    # We must set the compositions as distinct so that the equilibrium
    # function approaches the minimum with two distinct compositions,
    # rather than the metastable equilibrium where both phases have
    # the same composition.
    assemblage.phases[0].set_composition([0.05, 0.95])
    assemblage.phases[1].set_composition([0.95, 0.05])

    # Run the equilibration
    pressure = 1.e5
    temperatures = np.linspace(300., 601.4, 41)
    equality_constraints = [('P', pressure), ('T', temperatures)]

    sols, prm = equilibrate(composition, assemblage, equality_constraints)

    # Interrogate the stable assemblages for phase compositions.
    x1s = np.array([sol.assemblage.phases[0].molar_fractions[1]
                    for sol in sols if sol.code == 0])
    x2s = np.array([sol.assemblage.phases[1].molar_fractions[1]
                    for sol in sols if sol.code == 0])
    Ts = np.array([sol.assemblage.temperature
                   for sol in sols if sol.code == 0])

    plt.plot(x1s, Ts)
    plt.plot(x2s, Ts)

    plt.text(0.5, 400., 'miscibility gap', horizontalalignment='center')
    plt.xlabel('Molar proportion of pyrope')
    plt.ylabel('Temperature (K)')
    plt.show()

if __name__ == "__main__" and run_fper_ol:
    """
    Calculates the equilibrium Mg-Fe partitioning between
    ferropericlase and olivine.
    """
    temperatures = np.linspace(1000., 1500., 11)
    pressures = np.linspace(0.e9, 10.e9, 3)

    assemblage = burnman.Composite([ol, fper], [0.7, 0.3])

    ol.set_composition([0.93, 0.07])
    fper.set_composition([0.9, 0.1])

    assemblage.set_state(pressures[0], temperatures[0])
    equality_constraints = [('P', pressures), ('T', temperatures)]

    composition = {'Mg': 1., 'Fe': 0.5, 'Si': 0.5, 'O': 2.5}

    sol, prm = equilibrate(composition, assemblage, equality_constraints)

    for sols in sol:
        PGPa = sols[0].assemblage.pressure / 1.e9
        T = [s.assemblage.temperature for s in sols]
        x_fa = [s.assemblage.phases[0].molar_fractions[1] for s in sols]
        x_wus = [s.assemblage.phases[1].molar_fractions[1] for s in sols]
        p = plt.plot(T, x_fa, label=f'p(fa), {PGPa} GPa')
        plt.plot(T, x_wus, c=p[0].get_color(), linestyle='--',
                 label=f'p(wus), {PGPa} GPa')
    plt.xlabel('Temperature (K)')
    plt.ylabel('proportion (mol fraction)')
    plt.legend()
    plt.show()

if __name__ == "__main__" and run_fixed_ol_composition:
    """
    Calculates the composition of wadsleyite in equilibrium with olivine
    of a fixed composition at a fixed pressure.
    """
    assemblage = burnman.Composite([ol, wad], [0.7, 0.3])
    ol.set_composition([0.5, 0.5])
    wad.set_composition([0.6, 0.4])

    assemblage.set_state(10.e9, 1200.)
    equality_constraints = [('P', 10.e9),
                            ('phase_composition',
                             (ol, (['Mg_A', 'Fe_A'],
                                   [0., 1.], [1., 1.],
                                   0.45)))]
    composition = {'Mg': 1., 'Fe': 1., 'Si': 1., 'O': 4.}

    sol, prm = equilibrate(composition, assemblage, equality_constraints)

    print(assemblage)

if __name__ == "__main__" and run_upper_mantle:
    """
    Calculates the equilibrium compositions and phase proportions for
    an ol-opx-gt composite in an NCFMAS bulk composition.
    """
    pressure = 10.e9
    temperature = 1500.

    composition = {'Na': 0.02, 'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9,
                   'Ca': 0.2, 'Al': 0.4, 'O': 6.81}

    assemblage = burnman.Composite([ol, opx, gt], [0.7, 0.1, 0.2])

    ol.set_composition([0.93, 0.07])
    opx.set_composition([0.8, 0.1, 0.05, 0.05])
    gt.set_composition([0.8, 0.1, 0.05, 0.03, 0.02])

    assemblage.set_state(pressure, temperature)
    equality_constraints = [('P', 10.e9), ('T', 1500.)]

    sol, prm = equilibrate(composition, assemblage, equality_constraints,
                           max_iterations=20)
    print(assemblage)

if __name__ == "__main__" and run_lower_mantle:
    """
    Calculates temperatures and assemblage properties along an isentrope
    in the lower mantle. Includes calculations of the
    post-perovskite-in and bridgmanite-out lines.
    """
    P0 = 25.e9
    T0 = 1600.

    bdg.name = 'bdg'
    ppv.name = 'ppv'
    fper.name = 'fper'
    cpv.name = 'cpv'

    bdg.set_composition([0.86, 0.1, 0.04])
    ppv.set_composition([0.86, 0.1, 0.04])
    fper.set_composition([0.9, 0.1])

    composition = {'Fe': 0.2, 'Mg': 2.0, 'Si': 1.9, 'Ca': 0.2, 'Al': 0.4,
                   'O': 6.8}
    assemblage = burnman.Composite([bdg, fper, cpv])

    # We start by calculating the assemblage entropy at a reference
    # pressure and temperature.
    equality_constraints = [('P', P0), ('T', T0)]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)

    S = np.array([assemblage.molar_entropy*assemblage.n_moles])

    # Now we calculate the assemblages and states at the post-perovskite-in and
    # bridgmanite-out boundaries at the desired entropy.
    assemblage = burnman.Composite([bdg, fper, ppv, cpv])
    assemblage.set_state(sol.x[0], sol.x[1])
    equality_constraints = [('S', S),
                            ('phase_fraction', (ppv, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)
    P_ppv_in = assemblage.pressure
    T_ppv_in = assemblage.temperature

    assemblage.set_state(sol.x[0], sol.x[1])
    equality_constraints = [('S', S),
                            ('phase_fraction', (bdg, np.array([0.])))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)
    P_bdg_in = assemblage.pressure
    T_bdg_in = assemblage.temperature

    # Now we calculate the properties of each assemblage
    # along its stable pressure segments.
    prp_list = []
    prm_list = []
    for (pressures, phases) in [[np.linspace(25.e9, P_ppv_in, 21), [bdg, fper, cpv]],
                                [np.linspace(P_ppv_in, P_bdg_in, 21), [bdg, fper, ppv, cpv]],
                                [np.linspace(P_bdg_in, 140.e9, 21), [ppv, fper, cpv]]]:

        assemblage = burnman.Composite(phases)
        equality_constraints = [('P', pressures), ('S', S)]
        sols, prm = equilibrate(composition, assemblage,
                                equality_constraints)
        prm_list.append(prm)
        prp_list.append(np.array([sol.x for sol in sols if sol.success]).T)

    # Finally, do some plotting
    fig = plt.figure(figsize=(16, 8))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]
    phases = [bdg, fper, cpv, ppv]
    colors = ['red', 'green', 'blue', 'orange', 'purple']
    for j in range(3):
        prm = prm_list[j]
        p = prp_list[j]
        pressures, temperatures = p[[0, 1], :]
        ax[0].plot(pressures/1.e9, temperatures, color='black')
        ax[0].set_xlabel('Pressure (GPa)')
        ax[0].set_ylabel('Temperature (K)')
        for i, phase in enumerate(phases):
            try:
                idx = prm.parameter_names.index(f'x({phase.name})')
                x_phase = p[idx, :]
                ax[i+1].plot(pressures/1.e9, x_phase, color=colors[i])
                ax[i+1].set_ylim(0, 2)
                ax[i+1].set_xlim(0, 140)
                ax[i+1].set_xlabel('Pressure (GPa)')
                ax[i+1].set_ylabel(f'n moles of {phase.name}')
            except ValueError:
                pass

        if ((f'x({phases[0].name})' in prm.parameter_names
             and f'x({phases[1].name})' in prm.parameter_names)):
            pv_idx = prm.parameter_names.index(f'x({phases[0].name})')
            per_idx = prm.parameter_names.index(f'x({phases[1].name})')
            KD = (p[pv_idx+1, :] * (1. - p[per_idx+1, :])
                  / ((1. - p[pv_idx+1, :]
                      - p[pv_idx+2, :])*p[per_idx+1, :]))
            if j == 1:
                ax[5].plot(pressures/1.e9, KD, color='red', label='bdg K$_D$')
            else:
                ax[5].plot(pressures/1.e9, KD, color='red')

        if ((f'x({phases[1].name})' in prm.parameter_names
             and f'x({phases[3].name})' in prm.parameter_names)):
            per_idx = prm.parameter_names.index(f'x({phases[1].name})')
            ppv_idx = prm.parameter_names.index(f'x({phases[3].name})')
            KD = (p[ppv_idx+1, :] * (1. - p[per_idx+1, :])
                  / ((1. - p[ppv_idx+1, :]
                      - p[ppv_idx+2, :])*p[per_idx+1, :]))
            if j == 1:
                ax[5].plot(pressures/1.e9, KD, color='blue', label='ppv K$_D$')
            else:
                ax[5].plot(pressures/1.e9, KD, color='blue')

    ax[5].set_ylim(0., 1.)
    ax[5].set_xlabel('Pressure (GPa)')
    ax[5].set_ylabel('[FeSiO$_3$/MgSiO$_3$]/[FeO/MgO]')
    ax[5].legend()

    plt.show()

if __name__ == "__main__" and run_olivine_polymorphs:
    """
    Produces a P-T pseudosection for a fo90 composition.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    T0 = 573.15
    T1 = 1773.15

    ol.set_composition([0.93, 0.07])
    wad.set_composition([0.91, 0.09])
    rw.set_composition([0.93, 0.07])
    bdg.set_composition([0.86, 0.1, 0.04])
    fper.set_composition([0.9, 0.1])

    plt.text(1750., 10., 'olivine', horizontalalignment='right')
    plt.text(1750., 15., 'wadsleyite', horizontalalignment='right')
    plt.text(1750., 20., 'ringwoodite', horizontalalignment='right')
    plt.text(1750., 24., 'bridgmanite + ferropericlase',
             horizontalalignment='right')

    composition = {'Fe': 0.2, 'Mg': 1.8, 'Si': 1.0, 'O': 4.0}
    assemblage = burnman.Composite([ol, wad, rw], [1., 0., 0.])
    equality_constraints = [('phase_fraction', (ol, 0.0)),
                            ('phase_fraction', (rw, 0.0))]

    sol, prm = equilibrate(composition, assemblage, equality_constraints)
    Pinv1, Tinv1 = sol.x[0:2]

    if sol.code != 0:
        raise Exception("Couldn't find ol-wad-rw invariant")

    assemblage = burnman.Composite([ol, wad, rw])
    equality_constraints = [('phase_fraction', (wad, 0.0)),
                            ('phase_fraction', (rw, 0.0))]
    sol, prm = equilibrate(composition, assemblage, equality_constraints)
    Pinv2, Tinv2 = sol.x[0:2]

    for d in [[np.linspace(T0, Tinv1, 8), [ol, wad, rw], ol, 0.],
              [np.linspace(T0, Tinv2, 8), [ol, wad, rw], wad, 0.],
              [np.linspace(Tinv2, Tinv1, 8), [ol, wad, rw], rw, 0.],
              [np.linspace(T0, Tinv2, 8), [ol, rw], rw, 0.],
              [np.linspace(Tinv2, T1, 8), [ol, wad], wad, 0.],
              [np.linspace(Tinv1, T1, 8), [ol, wad], wad, 1.],
              [np.linspace(Tinv1, T1, 8), [wad, rw], rw, 0.]]:

        temperatures, phases, phase, fraction = d

        assemblage = burnman.Composite(phases)
        equality_constraints = [('T', temperatures),
                                ('phase_fraction', (phase, fraction))]
        sols, prm = equilibrate(composition, assemblage, equality_constraints)
        Ps = np.array([sol.x[0] for sol in sols if sol.success])
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        ax.plot(Ts, Ps/1.e9, color='k')

    assemblage.set_fractions([0., 1.])
    temperatures = np.linspace(T0, T1, 8)
    equality_constraints = [('T', temperatures),
                            ('phase_fraction', (rw, 1.0))]
    sols, prm = equilibrate(composition, assemblage, equality_constraints)
    Ps = np.array([sol.x[0] for sol in sols if sol.success])
    Ts = np.array([sol.x[1] for sol in sols if sol.success])
    ax.plot(Ts, Ps/1.e9, color='k')

    assemblage = burnman.Composite([rw, bdg, fper])
    assemblage = simplify_composite_with_composition(assemblage,
                                                     composition)

    assemblage.phases[1].set_composition([0.86, 0.14])

    for phase_fraction in [0., 1.]:
        equality_constraints = [('T', temperatures),
                                ('phase_fraction', (rw, phase_fraction))]
        sols, prm = equilibrate(composition, assemblage, equality_constraints)
        Ps = np.array([sol.x[0] for sol in sols if sol.success])
        Ts = np.array([sol.x[1] for sol in sols if sol.success])
        ax.plot(Ts, Ps/1.e9, color='k')


    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Pressure (GPa)')
    ax.set_ylim(5., 25.)
    plt.show()
