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

    def three_phase_eqm(args, P, T, composition):
        n_ppv, x_al_ppv, x_fe_ppv, x_fe_bdg = args

        n_bdg = (composition['Si'] + 0.5*composition['Al']) - n_ppv
        n_per = composition['O'] - 3.*(n_ppv + n_bdg)

        x_mg_ppv = 1. - x_al_ppv - x_fe_ppv

        x_al_bdg = (0.5*composition['Al'] - n_ppv*x_al_ppv)/n_bdg
        x_mg_bdg = 1. - x_al_bdg - x_fe_bdg
        
        x_fe_per = (composition['Fe'] - n_ppv*x_fe_ppv - n_bdg*x_fe_bdg)/n_per
        x_mg_per = 1. - x_fe_per

        amounts = np.array([n_ppv, n_bdg, n_per])
        
        ppv.set_composition([x_mg_ppv, x_fe_ppv, x_al_ppv])
        bdg.set_composition([x_mg_bdg, x_fe_bdg, x_al_bdg])
        fper.set_composition([x_mg_per, x_fe_per])
        ppv.set_state(P, T)
        bdg.set_state(P, T)
        fper.set_state(P, T)

        # Mg amount
        # Al2O3 activities
        # FeSiO3 activities (bdg, ppv)
        # Mg-Fe exchange (bdg, per)
        return [ppv.partial_gibbs[0] - bdg.partial_gibbs[0],
                ppv.partial_gibbs[2] - bdg.partial_gibbs[2],
                (ppv.partial_gibbs[1] + fper.partial_gibbs[0]) - (ppv.partial_gibbs[0] + fper.partial_gibbs[1]),
                (bdg.partial_gibbs[1] + fper.partial_gibbs[0]) - (bdg.partial_gibbs[0] + fper.partial_gibbs[1])]


        
    
    # Let's solve the equations at a fixed temperature of 2000 K
    # and Mg-rich composition with a little bit of Al2O3
    T = 2000.
    composition = { 'Mg': 1.775, 'Fe': 0.2, 'Al': 0.05, 'Si': 0.975, 'O': 4.}
    bdg = SLB_2011.mg_fe_bridgmanite() # Mg, Fe, Al
    ppv = SLB_2011.post_perovskite() # Mg, Fe, Al
    fper = SLB_2011.ferropericlase() # Mg, Fe
    assemblage1 = burnman.Composite([bdg, fper])
    assemblage2 = burnman.Composite([ppv, bdg, fper])
    assemblage3 = burnman.Composite([ppv, fper])
    
    # The next few lines do all the work, looping over lower mantle pressures
    # and finding the equilibrium composition at each P-T point.
    pressures = np.linspace(20.e9, 140.e9, 101)

    ppv_in = find_univariant(composition, [bdg, fper], ppv, 'T', [T], [100.e9, T, 0.0, 0.005, 0.025, 1.0, 0.005, 0.025, 1., 0.19525])
    pv_out = find_univariant(composition, [ppv, fper], bdg, 'T', [T], [100.e9, T, 0.0, 0.005, 0.025, 1.0, 0.005, 0.025, 1., 0.19525])

    P0 = ppv_in[0][0]
    plt.plot([P0, P0, P0], [ppv_in[0][3], ppv_in[0][6], ppv_in[0][9]], marker='o', linestyle='None')
    P1 = pv_out[0][0]
    plt.plot([P1, P1, P1], [pv_out[0][3], pv_out[0][6], pv_out[0][9]], marker='o', linestyle='None')

    P_bdg = []
    x_bdg = []
    P_ppv = []
    x_ppv = []
    P_fper = []
    x_fper = []
    guess_three_phase = [0.1, 0.03, 0.03, 0.03]
    for i, P in enumerate(pressures):

        if P < P0:
            sol = gibbs_minimizer(composition, assemblage1, [['P', P], ['T', T]])

            P_bdg.append(P)
            x_bdg.append(sol['c'][1])
            P_fper.append(P)
            x_fper.append(sol['c'][4])
            
        elif P > P1:
            sol = gibbs_minimizer(composition, assemblage3, [['P', P], ['T', T]])

            P_ppv.append(P)
            x_ppv.append(sol['c'][1])
            P_fper.append(P)
            x_fper.append(sol['c'][4])
            
        else:
            sol = opt.fsolve(three_phase_eqm, guess_three_phase, args=(P, T, composition), full_output=True)
            guess_three_phase = sol[0]
            P_ppv.append(P)
            P_bdg.append(P)
            P_fper.append(P)
            zeros = three_phase_eqm(guess_three_phase, P, T, composition)
            x_ppv.append(ppv.molar_fractions[1])
            x_bdg.append(bdg.molar_fractions[1])
            x_fper.append(fper.molar_fractions[1])
            


    # Let's print out a bit of information
    print('Example 3: Lower mantle phase relations in peridotite')
    print('Composition: {}'.format(composition))
    print('At {0:.1f} GPa and {1:.0f} K, ferropericlase contains {2:.1f} mole percent FeO'.format(pressures[0]/1.e9, T, 100.*x_fper[0]))
    print('At {0:.1f} GPa and {1:.0f} K, this has changed to {2:.1f} mole percent FeO'.format(pressures[-1]/1.e9, T, 100.*x_fper[-1]))
    print()
    
    # Finally, we plot the coexisting compositions 
    plt.plot(P_bdg, x_bdg, label='x(Fe) bdg')
    plt.plot(P_ppv, x_ppv, label='x(Fe) ppv')
    plt.plot(P_fper, x_fper, label='x(Fe) fper')
    plt.legend(loc='upper right')
    plt.title('Bridgmanite-ferropericlase compositions at '+str(T)+' K')
    plt.show()
    
    
