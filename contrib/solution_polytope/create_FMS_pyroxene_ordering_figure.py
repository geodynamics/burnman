# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman import minerals

from pylab import rcParams
plt.rcParams['font.family'] = "sans-serif"

# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


def equilibrate_px(mineral, x_Fe, P, T):
    """
    Function to find the equilibrium state of order of a two site pyroxene
    with site distributions [Fe][Fe], [Mg][Mg] and [Fe][Mg]
    """

    # Try completely disordered as a starting guess
    mineral.set_state(P, T)

    if x_Fe < 0.5:
        aord = np.array([0., 1. - 2*x_Fe, 2.*x_Fe])
        bord = np.array([2.*x_Fe, 1., -2.*x_Fe])
    else:
        x_Mg = 1. - x_Fe
        aord = np.array([1. - 2*x_Mg, 0., 2.*x_Mg])
        bord = np.array([1., 2.*x_Mg, -2.*x_Mg])

    def diff_pot(Q):
        c = aord*Q + bord*(1. - Q)
        mineral.set_composition(c)

        return (mineral.partial_gibbs[0] + mineral.partial_gibbs[1])/2. - mineral.partial_gibbs[2]

    root = brentq(diff_pot, 0.01, 1.) # make slightly asymmetric
    return root



if __name__ == "__main__":

    """
    First, we create three FMS (clino)pyroxene solution models,
    one which is globally disordered, one which undergoes convergent ordering,
    and the one from Holland et al., 2018.
    """

    cen = minerals.HP_2011_ds62.en()
    cfs = minerals.HP_2011_ds62.fs()

    disordered = burnman.SolidSolution(name = 'ordered phases unstable',
                                  solution_type = 'symmetric',
                                  endmembers = [[cfs,
                                                 '[Fe][Fe]Si2O6'],
                                                [cen,
                                                 '[Mg][Mg]Si2O6'],
                                                [burnman.CombinedMineral([cen, cfs], [0.5, 0.5], [2.e3, 0., 0.]),
                                                 '[Fe][Mg]Si2O6']],
                                  energy_interaction = [[2.3e3, -0.85e3],
                                                        [-0.85e3]]) # convergent, but ordered phases unstable


    """
    The nonconvergent model is taken from Holland et al (2018)
    The solution model file used by thermocalc contains the following lines:

    W(fs,cen)     2.3  0   0
    W(fs,cfm)     3.5  0   0
    W(cen,cfm)      4  0   0

    and

    cfs
      make 1 fs 1
      DQF    2.1      -0.002      0.045

    cen
      make 1 en 1
      DQF  3.5         -0.002       0.048

    cfm
      make  2       en  1/2   fs  1/2
      DQF   -1.6     -0.002      0.0465


    These lines imply that the ordering reaction
    1/2 cen + 1/2 cfs -> cfm
    is associated with a reaction energy of -4.2 kJ/mol

    We can now make the solution model:
    """
    fm_HGP2018 = burnman.CombinedMineral([cen,
                                          cfs],
                                          [0.5, 0.5],
                                          [-4.2e3, 0., 0.])
    nonconvergent_HGP2018 = burnman.SolidSolution(name = 'nonconvergent ordering',
                                  solution_type = 'symmetric',
                                  endmembers = [[cfs,
                                                 '[Fe][Fe]Si2O6'],
                                                [cen,
                                                 '[Mg][Mg]Si2O6'],
                                                [burnman.CombinedMineral([cen, cfs], [0.5, 0.5], [-4.2e3, 0., 0.]),
                                                 '[Fe][Mg]Si2O6']],
                                  energy_interaction = [[2.3e3, 3.5e3],
                                                        [4.e3]])

    """
    To create the convergent model, we need to
    understand some relations between energies and interaction parameters

    W'
    W(cen, cfs) = W(cen, cfs)
    W(cen, cmf) = W(cfs, cfm)
    W(cfs, cmf) = W(cen, cfm)

    G'
    G(cmf) = G(cen) - G(cfm) + G(cfs) - W(cen, cfm) + W(cen, cfs) - W(cfs, cfm)
    i.e. DQF(cmf) = -DQF(cfm) - W(cen, cfm) + W(cen, cfs) - W(cfs, cfm)

    For the nonconvergent HGP2018 solution
    DQF(cfm) = -4.2
    W(cen, cfm) = 4
    W(cen, cfs) = 2.3
    W(cfs, cfm) = 3.5
    so DQF(cmf) = 4.2 - 4 + 2.3 - 3.5 = -1

    To make a convergent solution (DQF(cmf) = DQF(cfm)) with interactions
    of a similar order of magnitude, we use the expression
    DQF(cmf) = (- W(cen, cfm) + W(cen, cfs) - W(cfs, cfm)) / 2

    let DQF(cmf) = -4.2 kJ/mol
    then
    -4.2*2 kJ/mol - 2.3 kJ/mol = -10.7 kJ/mol = (-W(cen,cfm) - W(cfs,cfm))
    thus
    W(cen,cfm) = W(cfs,cfm) = 5.35 kJ/mol
    """

    convergent = burnman.SolidSolution(name = 'convergent ordering',
                                  solution_type = 'symmetric',
                                  endmembers = [[cfs,
                                                 '[Fe][Fe]Si2O6'],
                                                [cen,
                                                 '[Mg][Mg]Si2O6'],
                                                [burnman.CombinedMineral([cen, cfs], [0.5, 0.5], [-4.2e3, 0., 0.]),
                                                 '[Fe][Mg]Si2O6']],
                                  energy_interaction = [[2.3e3, 5.35e3],
                                                        [5.35e3]])


    """
    Now we create the desired figures using these models
    """

    pMg1s = np.linspace(0., 1., 201)
    pMg2s = np.linspace(0., 1., 201)

    cen.set_state(1.e5, 100.)
    cfs.set_state(1.e5, 100.)
    H_cfs = cfs.molar_enthalpy
    H_cen = cen.molar_enthalpy


    fig = plt.figure(figsize=(8,8))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]


    pyroxenes = [disordered, convergent, nonconvergent_HGP2018]
    Hex = [np.zeros((201, 201)) for px in pyroxenes]
    Sex = np.zeros((201, 201))
    cs = [0, 0, 0, 0]
    for i in range(201):
        print(f'{i:03d}/201\r', end='')
        for j in range(201):
            p_en = pMg1s[i]
            p_fs = (1. - pMg2s[j])
            p_fm = 1. - p_en - p_fs

            p_Mg = (pMg1s[i] + pMg2s[j])/2.
            H_base = H_cfs*(1. - p_Mg) + H_cen*p_Mg

            if True:
                for px_idx, px in enumerate(pyroxenes):
                    px.set_composition([p_fs, p_en, p_fm])
                    px.set_state(1.e5, 100.)

                    Hex[px_idx][i][j] = px.molar_enthalpy  - H_base

                    if px_idx == 0:
                        Sex[i][j] = px.excess_entropy

    pMg1, pMg2 = np.meshgrid(pMg1s, pMg2s)
    lvl1 = np.linspace(-4000., 2000., 13)
    cs[0] = ax[0].contour(pMg2, pMg1, Hex[0], lvl1, cmap='viridis')
    cs[1] = ax[1].contour(pMg2, pMg1, Hex[1], lvl1, cmap='viridis')
    cs[2] = ax[2].contour(pMg2, pMg1, Hex[2], lvl1, cmap='viridis')
    cs[3] = ax[3].contour(pMg2, pMg1, Sex, [0, 2, 4, 6, 8, 10, 12], cmap='viridis')

    lbl = [f'a) $\\mathcal{{G}}$* ({pyroxenes[0].name})',
           f'b) $\\mathcal{{G}}$* ({pyroxenes[1].name})',
           f'c) $\\mathcal{{G}}$* ({pyroxenes[2].name})',
           'd) $S_{conf}$ (all cases)']

    for i in range(4):
        disordered_line, = ax[i].plot([0., 1.], [0., 1.], color='black', zorder=100)
        disordered_line.set_label('complete disorder')

        for j in np.linspace(0.2, 0.8, 4):
            ax[i].plot([0, j], [j, 0], color='grey', linestyle=':')
            ax[i].plot([1, j], [j, 1], color='grey', linestyle=':')

        ax[i].plot([0, 1], [1, 0], color='grey', linestyle=':')


    # Plot labelled boxes in the corners of the subplots
    for i in [0,1,2,3]:
        bbox = ax[i].fill_between([0.01, 0.11], [0.01, 0.01], [0.09, 0.09], color='white', linestyle='-', zorder=110)
        bbox.set_edgecolor('k')
        ax[i].text(0.06, 0.045, 'cfs', horizontalalignment='center', verticalalignment='center', color='black', zorder=120)

        bbox = ax[i].fill_between([0.01, 0.11], [0.91, 0.91], [0.99, 0.99], color='white', linestyle='-', zorder=110)
        bbox.set_edgecolor('k')
        ax[i].text(0.06, 0.945, 'cfm', horizontalalignment='center', verticalalignment='center', color='black', zorder=120)


        bbox = ax[i].fill_between([0.89, 0.99], [0.01, 0.01], [0.09, 0.09], color='white', linestyle='-', zorder=110)
        bbox.set_edgecolor('k')
        ax[i].text(0.94, 0.045, 'cmf', horizontalalignment='center', verticalalignment='center', color='black', zorder=120)

        bbox = ax[i].fill_between([0.89, 0.99], [0.91, 0.91], [0.99, 0.99], color='white', linestyle='-', zorder=110)
        bbox.set_edgecolor('k')
        ax[i].text(0.94, 0.945, 'cen', horizontalalignment='center', verticalalignment='center', color='black', zorder=120)

        # Plot text
        ax[i].text(0.5, 0.5, 'complete disorder', fontsize=10,
                   rotation=45, rotation_mode='anchor',
                   verticalalignment='bottom',
                   horizontalalignment='center')

        # Recast levels to new class
        cs[i].levels = [nf(val) for val in cs[i].levels]
        ax[i].clabel(cs[i], inline=1, fmt = '%r', fontsize=10)
        #fig.colorbar(cs[i], ax=ax[i])
        ax[i].set_title(lbl[i], loc='left')
        ax[i].set_xlabel('p(Mg on Site 1)')
        ax[i].set_ylabel('p(Mg on Site 2)')


    x_Fes = np.linspace(0.001, 0.999, 501)
    pMg1s = np.empty_like(x_Fes)
    pMg2s = np.empty_like(x_Fes)


    # calculate the equilibrium distribution of site-species
    # for the convergent and nonconvergent models at different temperatures.
    ls = ['--', ':']
    cs = ['orange', 'red']
    for j_px in [0, 1, 2]: # disordered, convergent and non convergent cases
        for k, T in enumerate([400., 500.]): # temperatures
            for i, x_Fe in enumerate(x_Fes):
                equilibrate_px(pyroxenes[j_px], x_Fe=x_Fe, P=1.e8, T=T) # none of the models are pressure dependent.
                pMg1s[i] = pyroxenes[j_px].molar_fractions[0]
                pMg2s[i] = 1. - pyroxenes[j_px].molar_fractions[1]

                if j_px == 1 and pMg1s[i] > pMg2s[i]:
                    a = pMg1s[i]
                    pMg1s[i] = pMg2s[i]
                    pMg2s[i] = a

            if j_px == 1:
                ax[j_px].plot(pMg2s, pMg1s, linestyle=ls[k], c=cs[k], zorder=105)

            ax[j_px].plot(pMg1s, pMg2s, linestyle=ls[k], c=cs[k], zorder=105)


    fig.tight_layout()
    fig.savefig('energy_entropy_2_site_ordering.pdf')
    plt.show()
