# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_composite
-----------------

This example demonstrates the functionalities of the burnman.Composite
class.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Mineral`
* :class:`burnman.Solution`
* :class:`burnman.Composite`

*Demonstrates:*

* How to initialize a composite object containing minerals and solutions
* How to set state and composition of composite objects
* How to interrogate composite objects for their
  compositional, thermodynamic and thermoelastic properties.
* How to use the stoichiometric and reaction affinity methods to solve
  simple thermodynamic equilibrium problems.
"""

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from scipy.optimize import brentq, fsolve
import matplotlib.pyplot as plt


import burnman
from burnman import minerals


if __name__ == "__main__":

    print('Part I: Making and interrogating composite materials')

    # In BurnMan, composite materials are those made up of a mechanical
    # mixture of other materials. Those materials could be minerals,
    # solutions or other composites.

    # Initialising a composite material is easy; first initialise all the
    # material objects you want to add to the composite, and then call the
    # constructor for the Composite class.
    # The following lines do this for a mixture of bridgmanite,
    # ferropericlase and Ca-perovskite.
    bdg = burnman.minerals.SLB_2011.mg_fe_bridgmanite()
    fper = burnman.minerals.SLB_2011.ferropericlase()
    cpv = burnman.minerals.SLB_2011.ca_perovskite()

    rock = burnman.Composite([bdg, fper, cpv],
                             fractions=[0.7, 0.2, 0.1],
                             fraction_type='molar',
                             name='lower mantle assemblage')

    # The fractions and fraction_type arguments are optional.
    # Some methods are immediately available to us. For example, we can
    # print the endmembers of all the phases.

    print('Names of endmembers in composite material:')
    for name in rock.endmember_names:
        print(name)

    # We can also print a list of the potential elements in the composite:
    print('Elements which might be in the composite:')
    print(rock.elements)

    # and a stoichiometric array
    table = [rock.elements]
    table.extend([list(s) for s in rock.stoichiometric_array])

    print('Stoichiometric array:')
    burnman.tools.misc.pretty_print_table(table)

    # Before we can interrogate our new composite fully, we must also set
    # the compositions of any solutions, and also set the state of the
    # composite.
    bdg.set_composition([0.88, 0.07, 0.05])
    fper.set_composition([0.85, 0.15])
    rock.set_state(30.e9, 2000.)

    # Now we can print the rock:
    print('Composite composition')
    els = burnman.tools.chemistry.sort_element_list_to_IUPAC_order(
        rock.formula.keys())
    for e in els:
        print(f'{e}: {rock.formula[e]:.3f}')

    # and also many properties of the rock:
    print(f'Density: {rock.density:.3f} kg/m^3')
    print(f'Entropy: {rock.molar_entropy:.3f} J/K/mol')
    print(f'P-wave velocity: {rock.p_wave_velocity/1.e3:.3f} km/s')

    # We can also evaluate these properties at many pressures and
    # temperatures using the evaluate method:
    pressures = np.linspace(30.e9, 130.e9, 101)
    temperatures = 2000.*np.ones_like(pressures)
    Vp, Vs, rho = rock.evaluate(['p_wave_velocity',
                                 'shear_wave_velocity',
                                 'density'],
                                pressures, temperatures)
    plt.plot(pressures/1.e9, Vp/1.e3, label='V$_p$ (km/s)')
    plt.plot(pressures/1.e9, Vs/1.e3, label='V$_s$ (km/s)')
    plt.plot(pressures/1.e9, rho/1.e3, label='Density (g/cm$^3$)')
    plt.legend()
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Velocities, Density')
    plt.show()

    print('Part II: Manual equilibration')
    print('Example 1: Reaction in an andalusite-sillimanite rock\n')
    andalusite = minerals.HP_2011_ds62.andalusite()
    sillimanite = minerals.HP_2011_ds62.sill()
    kyanite = minerals.HP_2011_ds62.ky()

    rock = burnman.Composite([andalusite, sillimanite],
                             name='and-sill assemblage')

    rock.set_state(1.e9, 600.)

    print('Isochemical reaction basis for rock')
    print(rock.reaction_basis_as_strings[0])

    print(f'Are we equilibrated at {rock.pressure/1.e9} GPa'
          f' and {rock.temperature} K? {rock.equilibrated}')

    def delta_affinity_two_polymorphs(P, T, rock):
        rock.set_state(P, T)
        return rock.reaction_affinities

    T = 600.
    print(f'Equilibrating the rock with brentq at {int(T)} K')
    brentq(delta_affinity_two_polymorphs, 1.e5, 1.e9, args=(T, rock))

    print(f'Are we equilibrated now? {rock.equilibrated}')
    print(f'Equilibrium pressure: {rock.pressure/1.e9:.3f} GPa')

    print('\nExample 2: Reaction in an andalusite-sillimanite-kyanite rock\n')
    rock = burnman.Composite([andalusite, sillimanite, kyanite],
                             name='and-sill-ky assemblage')

    def delta_affinity_three_polymorphs(x):
        P, T = x
        rock.set_state(P, T)
        return rock.reaction_affinities

    print('Equilibrating the rock with fsolve')
    fsolve(delta_affinity_three_polymorphs, [1.e9, 1600.])

    print(f'Are we equilibrated now? {rock.equilibrated}')
    print(f'Equilibrium pressure: {rock.pressure/1.e9:.3f} GPa')
    print(f'Equilibrium temperature: {rock.temperature:.3f} K')

    print('\nExample 3: Reaction in an olivine-wadsleyite rock\n')
    ol = burnman.minerals.SLB_2011.mg_fe_olivine()
    wad = burnman.minerals.SLB_2011.mg_fe_wadsleyite()

    rock = burnman.Composite([ol, wad], name='Olivine-wadsleyite assemblage')

    print(f'There are {rock.n_reactions} independent reactions in this rock.')
    print('One possible basis set is:')
    for string in rock.reaction_basis_as_strings:
        print(string)

    print('The reaction basis can also be returned as a numpy array:')
    print(rock.reaction_basis)

    print(f'If there are {rock.n_reactions} independent reactions, '
          f'we need to select {rock.n_reactions} unknowns.')
    print('In this example, we find the pressure and composition \n'
          'of wadsleyite which coexists with fo90 at 1600 K.\n')

    def delta_affinity_two_binaries_P_x_Fe_wad(x, T):
        P, x_fe_wad = x
        rock.set_state(P, T)
        wad.set_composition([1. - x_fe_wad, x_fe_wad])
        return rock.reaction_affinities

    print('Equilibrating the rock with fsolve')
    rock.set_fractions([1., 0.])
    ol.set_composition([0.9, 0.1])
    fsolve(delta_affinity_two_binaries_P_x_Fe_wad, [13.e9, 0.1], args=(1600.))
    print(f'Are we equilibrated now? {rock.equilibrated}')
    print(rock)

    print('\nWe could alternatively choose to solve for the compositions of '
          'olivine and wadsleyite at fixed pressure and temperature')

    # First, get the univariant pressure at the Mg-rich end of the binary
    # at a temperature of 1600 K.
    T = 1600.
    fo_mwd = burnman.Composite([burnman.minerals.SLB_2011.forsterite(),
                                burnman.minerals.SLB_2011.mg_wadsleyite()])
    brentq(delta_affinity_two_polymorphs, 12.e9, 16.e9, args=(T, fo_mwd))

    # Now solve for the equilibrium state at a range of pressures
    # less than the univariant pressure
    def delta_affinity_two_binaries_x_Fe_ol_wad(x, P, T):
        x_fe_ol, x_fe_wad = x
        rock.set_state(P, T)
        ol.set_composition([1. - x_fe_ol, x_fe_ol])
        wad.set_composition([1. - x_fe_wad, x_fe_wad])
        return rock.reaction_affinities

    pressures = np.linspace(10.e9, fo_mwd.pressure, 21)
    x_fe_ols = np.zeros_like(pressures)
    x_fe_wads = np.zeros_like(pressures)
    for i, P in enumerate(pressures[:-1]):
        sol = fsolve(delta_affinity_two_binaries_x_Fe_ol_wad,
                     [0.1, 0.1], args=(P, T))
        x_fe_ols[i], x_fe_wads[i] = sol

    # Pretty-print the results
    print(f'Olivine-wadsleyite equilibria at {T} K:')
    table = [['Pressure (GPa)', 'p(Fe2SiO4, ol)', 'p(Fe2SiO4, wad)']]
    for i in range(len(pressures)):
        table.append([f'{pressures[i]/1.e9:.3f}',
                      f'{x_fe_ols[i]:.3f}',
                      f'{x_fe_wads[i]:.3f}'])
    burnman.tools.misc.pretty_print_table(table)

    # Plot the results
    plt.plot(x_fe_ols, pressures/1.e9, label='ol')
    plt.plot(x_fe_wads, pressures/1.e9, label='wad')
    plt.legend()
    plt.xlim(0,)
    plt.xlabel('p(Fe$_2$SiO$_4$)')
    plt.ylabel('Pressure (GPa)')
    plt.show()
