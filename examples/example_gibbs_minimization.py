# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

'''
example_gibbs_minimization
--------------------

This example demonstrates how burnman may be used to calculate the
equilibrium phase proportions and compositions for an assemblage
of a fixed bulk composition.

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
from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from burnman.equilibriumassemblage import *

if __name__ == "__main__":
    # Let's solve for the composition of Martian crust at 1173 K (900 C) and 1 GPa
    P = 1.e9
    T = 1173.15
    atomic_masses = read_masses()
    wt_composition = { 'K': 0.3740, 'Na': 2.20, 'Ca': 4.95, 'Fe': 14.1, 'Mg': 5.46, 'Al': 5.56, 'Si': 23.0, 'Ti': 0.5880}
    oxidation_state = { 'K': 1., 'Na': 1., 'Ca': 2., 'Fe': 2., 'Mg': 2., 'Al': 3., 'Si': 4., 'Ti': 4.}
    moles = {key: value/formula_mass({key: 1.}, atomic_masses) for key, value in wt_composition.iteritems()}
    moles['Na'] += moles['K']
    del moles['K']
    moles['Si'] += moles['Ti']
    del moles['Ti']
    moles['O'] = np.sum([value*oxidation_state[key]/2. for key, value in moles.iteritems()])
    total = np.sum(moles.values())
    moles = {key: value / total for key, value in moles.iteritems()}
    mol_oxides = {'Na2O': moles['Na']/2., 'CaO': moles['Ca'], 'FeO': moles['Fe'],
                  'MgO': moles['Mg'], 'Al2O3': moles['Al']/2., 'SiO2': moles['Si']}
    total = np.sum(mol_oxides.values())
    mol_oxides = {key: value / total*100. for key, value in mol_oxides.iteritems()}
    print(mol_oxides) 
    
    ol = SLB_2011.mg_fe_olivine()
    opx = SLB_2011.orthopyroxene()
    pl = SLB_2011.plagioclase()
    cpx = SLB_2011.clinopyroxene()
    sp = SLB_2011.mg_fe_aluminous_spinel()
    gt = SLB_2011.garnet()
    q = SLB_2011.quartz()
    
    assemblage = burnman.Composite([sp, opx, pl, q])
    sol = gibbs_minimizer(moles, assemblage, [['P', P], ['T', T]])
    
    exit()

    
    # Example 1: The classic aluminosilicate diagram.
    # This is a T-P phase diagram with composition Al2SiO5,
    # which is the composition of andalusite, sillimanite and
    # kyanite.    
    andalusite = HP_2011_ds62.andalusite()
    sillimanite = HP_2011_ds62.sill()
    kyanite = HP_2011_ds62.ky()

    # Here we define the composition of interest
    composition = andalusite.params['formula']
    
    # First, we need to solve for the invariant point where the
    # three phases coexist. The function find_invariant takes a bulk composition
    # and then two variables which constitute the stable phases and the two phases
    # which become stable at that point.
    P_inv, T_inv = find_invariant(composition, [andalusite], [sillimanite, kyanite])[0:2]
    
    # Here we print the P-T position of the invariant point
    print('Example 1: The aluminosilicates')
    print('The and-sill-ky invariant can be found at {0:.2f} GPa and {1:.0f} K'.format(P_inv/1.e9, T_inv))
    print()

    # Now let's plot the full phase diagram
    # The kyanite-sillimanite reaction line appears at pressures above the invariant point
    lo_pressures = np.linspace(1.e5, P_inv, 21)
    hi_pressures = np.linspace(P_inv, 1.e9, 21)

    # For isochemical diagrams, we can make use of the function find_univariant,
    # which for a given composition and a pressure (temperature) range outputs the temperatures
    # (pressures) where one phase in an assemblage becomes unstable.
    and_ky_eqm = find_univariant(composition, [andalusite], kyanite, 'P', lo_pressures)
    and_sill_eqm = find_univariant(composition, [andalusite], sillimanite, 'P', lo_pressures)
    sill_ky_eqm = find_univariant(composition, [sillimanite], kyanite, 'P', hi_pressures)

    plt.plot(zip(*and_ky_eqm)[1], lo_pressures/1.e9, label='and-ky')
    plt.plot(zip(*and_sill_eqm)[1], lo_pressures/1.e9, label='and-sill')
    plt.plot(zip(*sill_ky_eqm)[1], hi_pressures/1.e9, label='sill-ky')
    
    plt.title('Aluminosilicate phase diagram')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (GPa)')
    plt.legend(loc='upper left')
    plt.show()
    
    
    # Example 2: The olivine phase diagram at 1400 K
    # Our aim here is to plot the X-P phase diagram
    # for Mg2SiO4-Fe2SiO4 olivine at a fixed temperature
    # For this example we illustrate the lower level functions
    # behind find_invariant and find_univariant
    
    # Here we define the temperature and assemblage
    T = 1673.
    
    ol = SLB_2011.mg_fe_olivine()
    wad = SLB_2011.mg_fe_wadsleyite()
    rw = SLB_2011.mg_fe_ringwoodite()
    pv = SLB_2011.mg_fe_bridgmanite()
    fper = SLB_2011.ferropericlase()
    stv = SLB_2011.stishovite()
    assemblage = burnman.Composite([ol, wad, rw])

    
    # We are interested in reactions at various compositions,
    # so here we define the compositions bounding the binary
    composition1 = { 'Mg': 2., 'Si': 1., 'O': 4.}
    composition2 = { 'Fe': 2., 'Si': 1., 'O': 4.}

    '''
    # The following illustrates some starting-composition-sensitive solves. 
    # It would be quite nice to make the solver more robust with respect to these guys...  
    composition = binary_composition(composition1, composition2, 0.3)
    pressure0 = find_univariant(composition, [fper, rw], stv, 'T', [T], [2.e10, T, 0., 0., 0.5, 1., 0.3])[0]
    print(pressure0/1.e9)

    guesses = [23.6e9, T, 0.0, 0.1, 0.6, 2., 0.3, 0., 0.1]
    pressure = find_univariant(composition, [stv, fper, rw], pv, 'T', [T], guesses)[0]
    print(pressure/1.e9)

    #pressure1 = find_univariant(composition, [stv, fper], rw, 'T', [T])[0]
    #print(pressure1/1.e9)
    #pressure = find_univariant(composition, [stv, fper], pv, 'T', [T])[0]
    #print(pressure/1.e9)
    
    pressure = find_univariant(composition, [pv, fper], stv, 'T', [T])[0]
    print(pressure/1.e9)

    pressures = np.linspace(pressure0 + 1.e6, pressure1 - 1.e6, 21)
    pwus0 = np.empty_like(pressures)
    pwus1 = np.empty_like(pressures)
    
    sfr = burnman.Composite([stv, fper, rw])
    sfp = burnman.Composite([stv, fper, pv])
    for i, P in enumerate(pressures):
        sfr_guesses = [P, T, 0., 0., 0.5, 1., 0.3]
        sfp_guesses = [P, T, 1., 2., 0.3, 0., 0.1]
        sol = gibbs_minimizer(composition, sfr, [['P', P], ['T', T]], sfr_guesses)
        pwus0[i] = sol['c'][sol['variables'].index('x(Stishovite)')]
        sol = gibbs_minimizer(composition, sfp, [['P', P], ['T', T]], sfp_guesses)
        pwus1[i] = sol['c'][sol['variables'].index('x(Stishovite)')]

    plt.plot(pressures/1.e9, pwus0)
    plt.plot(pressures/1.e9, pwus1)
    plt.show()
    exit()
    '''
    
    # First, we want to find the pressure at which olivine,
    # wadsleyite and ringwoodite coexist, and their compositions
    # There are many ways to do this, but here we fix the temperature
    # and the proportions of wadsleyite and ringwoodite to zero, and
    # find the bulk composition (which much therefore correspond to
    # the composition of olivine) at which all three phases coexist.

    # The gibbs minimizer takes compositional
    # constraints as a list with the first element being the string 'X'.
    # The second two elements are the multipliers for the compositional vector,
    # which has the form amount of phase 0, fraction of endmember 1 in phase 0,
    # fraction of endmember 2 in phase 0 ... amount of phase i ...
    # In this case, we are dealing with pure phases, so the compositional
    # vector is the amount of andalusite, sillimanite and kyanite.
    # These do not have to sum to one. The compositional contraints are
    # expressed in the form:
    # dot(constraint[1],cvector)/dot(constraint[2],cvector) = constraint[3]
    # For the invariant point, the two constraints are that the proportion
    # of two of the phases goes to zero.
    
    constraints = [['T', T],
                   ['X', [0., 0., 1., 0., 0., 0.], [1., 0., 1., 0., 1., 0.], 0.],
                   ['X', [0., 0., 0., 0., 1., 0.], [1., 0., 1., 0., 1., 0.], 0.]]


    # Here we call the gibbs bulk minimizer.
    sol = gibbs_bulk_minimizer(composition1, composition2, 0.25, assemblage, constraints)
    
    # The output of the bulk minimizer is a dictionary containing P, T, X
    # (where X is the distance along the compositional binary),
    # followed by the compositional vector
    P_inv = sol['P']
    x_ol_inv = sol['c'][sol['variables'].index('p(Fayalite)')]
    x_wad_inv = sol['c'][sol['variables'].index('p(Fe_Wadsleyite)')]
    x_rw_inv = sol['c'][sol['variables'].index('p(Fe_Ringwoodite)')]

       
    # Here we print the properties of the invariant point
    print('Example 2: The olivine polymorphs')
    print('The pressure of the olivine polymorph invariant at {0:.0f} K is {1:.1f} GPa'.format(T, P_inv/1.e9))
    print('The molar Fe2SiO4 contents of coexisting ol, wad and rw')
    print('at this point are {0:.2f}, {1:.2f} and {2:.2f}'.format(x_ol_inv, x_wad_inv, x_rw_inv))
    print()


    # Now we find the rest of the phase diagram.
    # The compositions of the phases at the invariant point
    # tell us the ranges of compositions of those phases
    # we should expect within the stable parts of each binary loop
    assemblages = [burnman.Composite([ol, wad]),
                   burnman.Composite([ol, rw]),
                   burnman.Composite([wad, rw])]

    c_bounds = [[0., x_ol_inv, 21],
                [x_ol_inv, 1., 21],
                [0., x_wad_inv, 21]]

    # For a given composition, find the pressure at which the
    # second phase of each binary loop becomes stable,
    # and the composition of the phase at that point.

    # Because each calculation is at a constant composition, we use the
    # function gibbs_minimizer.
    for i, assemblage in enumerate(assemblages):
        x1 = np.linspace(*c_bounds[i])
        x2 = np.empty_like(x1)
        pressures = np.empty_like(x1)
        guesses = None
        for i, x in enumerate(x1):
            composition = { 'Mg': 2.*(1. - x), 'Fe': 2.*x, 'O': 4., 'Si': 1.}
            sol = find_univariant(composition, [assemblage.phases[0]], assemblage.phases[1], 'T', [T], guesses)
            pressures[i] = sol[0][0]
            x2[i] = assemblage.phases[1].molar_fractions[1]
        plt.plot(x1, pressures/1.e9, label=assemblage.phases[1].name+'-in')
        plt.plot(x2, pressures/1.e9, label=assemblage.phases[0].name+'-in')

    # Here we finish off the plotting    
    plt.plot([x_ol_inv, x_rw_inv], [P_inv/1.e9, P_inv/1.e9])
    plt.title('Mg2SiO4-Fe2SiO4 phase diagram at '+str(T)+' K')
    plt.xlabel('p (Fe2SiO4)')
    plt.ylabel('Pressure (GPa)')
    plt.legend(loc='upper right')
    plt.show()
    


    # Example 3: Lower mantle phase relations
    # The lower mantle is dominated by bridgmanite and periclase,
    # which exchange Mg and Fe. We can use our minimization function
    # to find the compositions of the coexisting phases as a function
    # of temperature and pressure.
    
    # Let's solve the equations at a fixed temperature of 2000 K
    # and Mg-rich composition with a little bit of Al2O3
    T = 2000.
    composition = { 'Mg': 1.775, 'Fe': 0.2, 'Al': 0.05, 'Si': 0.975, 'O': 4.}
    bdg = SLB_2011.mg_fe_bridgmanite()
    fper = SLB_2011.ferropericlase()
    assemblage = burnman.Composite([bdg, fper])

    # The next few lines do all the work, looping over lower mantle pressures
    # and finding the equilibrium composition at each P-T point.
    pressures = np.linspace(22.e9, 130.e9, 21)
    x_bdg = np.empty_like(pressures)
    x_fper = np.empty_like(pressures)
    for i, P in enumerate(pressures):
        sol = gibbs_minimizer(composition, assemblage, [['P', P], ['T', T]])
        
        x_bdg[i] = sol['c'][1]
        x_fper[i] = sol['c'][4]

    # Let's print out a bit of information
    print('Example 3: Lower mantle phase relations in peridotite')
    print('Composition: {}'.format(composition))
    print('At {0:.1f} GPa and {1:.0f} K, ferropericlase contains {2:.1f} mole percent FeO'.format(pressures[0]/1.e9, T, 100.*x_fper[0]))
    print('At {0:.1f} GPa and {1:.0f} K, this has changed to {2:.1f} mole percent FeO'.format(pressures[-1]/1.e9, T, 100.*x_fper[-1]))
    print()
    
    # Finally, we plot the coexisting compositions 
    plt.plot(pressures, x_bdg, label='x(Fe) bdg')
    plt.plot(pressures, x_fper, label='x(Fe) fper')
    plt.legend(loc='upper right')
    plt.title('Bridgmanite-ferropericlase compositions at '+str(T)+' K')
    plt.show()
    
    
    # Example 4: Solvus in pyrope-grossular
    composition = { 'Ca': 1.5, 'Mg': 1.5, 'Al': 2., 'Si': 3., 'O': 12.}
    garnet0 = SLB_2011.garnet()
    garnet1 = SLB_2011.garnet()
    assemblage = burnman.Composite([garnet0, garnet1])
    P = 1.e5
    temperatures = []
    c_g1 = []
    c_g2 = []

    c1=0.021
    c2=0.979
    for T in np.linspace(300., 600., 21):
        constraints=[['P', P], ['T', T]]
        try:
            sol = gibbs_minimizer(composition, assemblage, constraints, guesses=[P, T, 0.5, c1, 0.5, c2])
            if np.abs(sol['c'][1] - sol['c'][3]) > 1.e-10:
                temperatures.append(T)
                c1=sol['c'][1]
                c2=sol['c'][3]
                c_g1.append(c1)
                c_g2.append(c2)
        except:
            print('No solution found at {0:.1f} K'.format(T))
            
    print('Example 4: Garnet immiscibility')
    print('At {0:.1f} GPa and {1:.0f} K, the pyrope-grossular immiscibility gap extends from {2:.3f} to {3:.3f} mole percent grossular'.format(P/1.e9, temperatures[0], c_g1[0], c_g2[0]))

            
    plt.plot(c_g1, temperatures)
    plt.plot(c_g2, temperatures)
    plt.xlabel(sol['variables'][1])
    plt.title('Pyrope-grossular immiscibility at '+str(P/1.e9)+' GPa')
    plt.show()
