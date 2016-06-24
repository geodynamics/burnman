# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_solid_solution
----------------------

This example shows how to create different solid solution models and output
thermodynamic and thermoelastic quantities.

There are four main types of solid solution currently implemented in
BurnMan:

1. Ideal solid solutions
2. Symmmetric solid solutions
3. Asymmetric solid solutions
4. Subregular solid solutions

These solid solutions can potentially deal with:

* Disordered endmembers (more than one element on a crystallographic site)
* Site vacancies
* More than one valence/spin state of the same element on a site

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.solidsolution.SolidSolution`
* :class:`burnman.solutionmodel.SolutionModel`


*Demonstrates:*

* Different ways to define a solid solution
* How to set composition and state
* How to output thermodynamic and thermoelastic properties

"""
from __future__ import absolute_import

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":
    '''
    First, let's pick a starting pressure and temperature
    '''
    P = 1.e5
    T = 1000.

    """
    Here we create an ideal solid solution for aluminous Mg,Fe garnets.
    Further endmembers can be added simply by adding more endmember lists.
    """
    class mg_fe_garnet(burnman.SolidSolution):

        def __init__(self, molar_fractions=None):
            self.name = 'Ideal pyrope-almandine garnet'
            self.type = 'ideal'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'], [
                               minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g1 = mg_fe_garnet()

    """
    Initialisation can optionally include setting the composition of the solid solution, i.e.
    """
    g1 = mg_fe_garnet([0.0, 1.0])
    g1.set_state(P, T)
    alm_gibbs = g1.gibbs

    """
    Alternatively, the composition can be set after initialisation
    """
    g1.set_composition([1.0, 0.0])
    g1.set_state(P, T)
    py_gibbs = g1.gibbs

    comp = np.linspace(0.001, 0.999, 100)
    g1_gibbs = np.empty_like(comp)
    g1_excess_gibbs = np.empty_like(comp)
    g1_pyrope_activity = np.empty_like(comp)
    g1_almandine_activity = np.empty_like(comp)
    g1_pyrope_gamma = np.empty_like(comp)
    g1_almandine_gamma = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, c]
        g1.set_composition(molar_fractions)
        g1.set_state(P, T)
        g1_gibbs[i] = g1.gibbs
        g1_excess_gibbs[i] = g1.excess_gibbs
        g1_pyrope_activity[i] = g1.activities[0]
        g1_almandine_activity[i] = g1.activities[1]
        g1_pyrope_gamma[i] = g1.activity_coefficients[0]
        g1_almandine_gamma[i] = g1.activity_coefficients[1]

    plt.plot([0., 1.], [py_gibbs / 1000., alm_gibbs / 1000.],
             'b-', linewidth=1., label='Mechanical mixing')
    plt.plot(comp, g1_gibbs / 1000., 'r-',
             linewidth=1., label='Ideal solid solution')
    plt.title("Ideal pyrope-almandine join")
    plt.ylabel("Gibbs free energy of solution (kJ/mol)")
    plt.xlabel("Almandine fraction")
    plt.legend(loc='lower right')
    plt.show()

    plt.plot(comp, g1_pyrope_activity, 'r-',
             linewidth=1., label='Pyrope activity')
    plt.plot(comp, g1_almandine_activity, 'b-',
             linewidth=1., label='Almandine activity')

    plt.plot(comp, g1_pyrope_gamma, 'r-.',
             linewidth=1., label='Pyrope activity coefficient')
    plt.plot(comp, g1_almandine_gamma, 'b--',
             linewidth=1., label='Almandine activity coefficient')
    plt.title("Ideal pyrope-almandine join")
    plt.ylim(0.0, 1.01)
    plt.ylabel("Activities")
    plt.xlabel("Almandine fraction")
    plt.legend(loc='lower right')
    plt.show()

    '''
    Not included in this example document are ways to create solid solutions with spin transitions,
    vacancies and mixed valence states. However, the formula parsing in BurnMan is comprehensive,
    so any model should be reproducable. Some formatted formulae examples follow:
    - [Fe]O: simple wuestite
    - [Fef2/3Vac1/3]O: a theoretical ferric endmember of wuestite. Note the distinct element name for Fe3+.
    - [Fe2/3Vac1/3]O: an alternative ferric endmember. Mixing between this and [Fe]O will produce a
      different configurational entropy to the previous model because there is no distinction
      between Fe2+ and Fe3+.
    - [Fels]O: Low spin wuestite. Another example illustrating the free-form approach to element definition.
    '''

    '''
    The solid solution corresponding to the pyrope-almandine join is actually not quite ideal.
    It can be well-approximated with a symmetric regular solid solution model
    '''
    class mg_fe_garnet_2(burnman.SolidSolution):

        def __init__(self, molar_fractions=None):
            self.name = 'Symmetric pyrope-almandine garnet'
            self.type = 'symmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'], [
                               minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]
            self.energy_interaction = [[2.5e3]]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g2 = mg_fe_garnet_2()
    g2_excess_gibbs = np.empty_like(comp)
    g2_pyrope_activity = np.empty_like(comp)
    g2_almandine_activity = np.empty_like(comp)
    g2_pyrope_gamma = np.empty_like(comp)
    g2_almandine_gamma = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1. - c, c]
        g2.set_composition(molar_fractions)
        g2.set_state(P, T)
        g2_excess_gibbs[i] = g2.excess_gibbs
        g2_pyrope_activity[i] = g2.activities[0]
        g2_almandine_activity[i] = g2.activities[1]
        g2_pyrope_gamma[i] = g2.activity_coefficients[0]
        g2_almandine_gamma[i] = g2.activity_coefficients[1]

    plt.plot(comp, g1_excess_gibbs, 'r-', linewidth=1., label='Ideal solution')
    plt.plot(comp, g2_excess_gibbs, 'g-', linewidth=1.,
             label='Symmetric solution, 2.5 kJ/mol')
    plt.title("Pyrope-almandine join (model comparison)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Almandine fraction")
    plt.legend(loc='upper left')
    plt.show()

    plt.plot(comp, g2_pyrope_activity, 'r-',
             linewidth=1., label='Pyrope activity')
    plt.plot(comp, g2_almandine_activity, 'b-',
             linewidth=1., label='Almandine activity')

    plt.plot(comp, g2_pyrope_gamma, 'r-.',
             linewidth=1., label='Pyrope activity coefficient')
    plt.plot(comp, g2_almandine_gamma, 'b--',
             linewidth=1., label='Almandine activity coefficient')
    plt.title("Non-ideal pyrope-almandine join")
    plt.ylabel("Activities")
    plt.xlabel("Almandine fraction")
    plt.legend(loc='lower right')
    plt.show()

    '''
    Adding more endmembers is very straightforward.
    Interaction terms must be added in the order [[12,13,14,...],[23,24,...],[34,...],...]
    Here, the new endmember is majorite, illustrating the addition of endmembers with
    multiple occupancy on a single site (here Mg and Si to replace Al)

    Here we also illustrate the addition of excess entropy and volume terms.
    V and S excess lists are optional parameters to the solid solution initialisation.
    In the initialisation, the order of interaction terms is H, V, S
    (S excesses are least commonly constrained).

    An important note: In Holland and Powell datasets, the "DQF" terms
    are equivalent to the excesses in burnman with the *exception* of the temperature correction,
    which is the negative of the entropy correction. i.e.
    DQF[0] = H_excess
    DQF[1] = -S_excess
    DQF[2] = V_excess
    '''
    class hp_mg_fe_garnet(burnman.SolidSolution):

        def __init__(self, molar_fractions=None):
            self.name = 'Symmetric pyrope-almandine-majorite garnet'
            self.type = 'symmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'], [
                               minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'], [minerals.HP_2011_ds62.maj(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
            self.energy_interaction = [[2.5e3, 0.0e3], [10.0e3]]
            self.entropy_interaction = [[0.0e3, 0.0e3], [0.0e3]]
            self.volume_interaction = [[0.0e3, 0.0e3], [0.0e3]]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g3 = hp_mg_fe_garnet()
    g3_configurational_entropy = np.empty_like(comp)
    g3_excess_entropy = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, 0., c]
        g3.set_composition(molar_fractions)
        g3.set_state(P, T)
        g3_configurational_entropy[
            i] = g3.solution_model._configurational_entropy(molar_fractions)
        g3_excess_entropy[i] = -np.dot(
            molar_fractions, g3.solution_model._ideal_excess_partial_gibbs(T, molar_fractions)) / T

    plt.plot(comp, g3_configurational_entropy, 'g-',
             linewidth=1., label='Configurational entropy')
    plt.plot(comp, g3_excess_entropy, 'r-',
             linewidth=1., label='Excess entropy')
    plt.title("Pyrope-majorite join")
    plt.ylabel("Entropy (J/K/mol)")
    plt.xlabel("Majorite fraction")
    plt.legend(loc='upper left')
    plt.show()

    '''
    The addition of grossular garnet (Ca-Al garnet) illustrates the use of asymmetric solution models
    This model is also found in the HP_2011_ds62 database
    '''

    class CFMASO_garnet(burnman.SolidSolution):

        def __init__(self, molar_fractions=None):
            self.name = 'garnet'
            self.endmembers = [[minerals.HP_2011_ds62.py(
            ), '[Mg]3[Al]2Si3O12'],
                [minerals.HP_2011_ds62.alm(
                ), '[Fe]3[Al]2Si3O12'],
                [minerals.HP_2011_ds62.gr(
                ), '[Ca]3[Al]2Si3O12'],
                [minerals.HP_2011_ds62.andr(), '[Ca]3[Fe]2Si3O12']]
            self.type = 'asymmetric'
            self.alphas = [1.0, 1.0, 2.7, 2.7]
            self.energy_interaction = [[2.5e3, 31.e3, 53.2e3],
                                        [5.e3, 37.24e3],
                                        [2.e3]]
            burnman.SolidSolution.__init__(self, molar_fractions)

    g4 = burnman.minerals.HP_2011_ds62.CFMASO_garnet()
    g4_excess_gibbs_400 = np.empty_like(comp)
    g4_excess_gibbs_800 = np.empty_like(comp)
    g4_excess_gibbs_1200 = np.empty_like(comp)

    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, 0., c, 0.]
        g4.set_composition(molar_fractions)
        g4.set_state(P, 400.)
        g4_excess_gibbs_400[i] = g4.excess_gibbs
        g4.set_state(P, 800.)
        g4_excess_gibbs_800[i] = g4.excess_gibbs
        g4.set_state(P, 1200.)
        g4_excess_gibbs_1200[i] = g4.excess_gibbs

    plt.plot(comp, g4_excess_gibbs_400, 'r-', linewidth=1., label='400 K')
    plt.plot(comp, g4_excess_gibbs_800, 'g-', linewidth=1., label='800 K')
    plt.plot(comp, g4_excess_gibbs_1200, 'b-', linewidth=1., label='1200 K')
    plt.title("Pyrope-grossular join (asymmetric model)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Grossular fraction")
    plt.legend(loc='lower left')
    plt.show()

    '''
    The subregular solution model (Helffrich and Wood, 1989)
    provides a more flexible way of constructing an asymmetric
    model
    '''
    class mg_fe_ca_garnet_Ganguly(burnman.SolidSolution):

        def __init__(self, molar_fractions=None):
            self.name = 'Subregular pyrope-almandine-grossular garnet'
            self.type = 'subregular'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'], [minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'], [
                               minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12'], [minerals.HP_2011_ds62.spss(), '[Mn]3[Al]2Si3O12']]
            self.energy_interaction = [
                [[2117., 695.], [9834., 21627.], [12083., 12083.]], [[6773., 873.], [539., 539.]], [[0., 0.]]]
            self.volume_interaction = [[[0.07e-5, 0.], [0.058e-5, 0.012e-5], [0.04e-5, 0.03e-5]], [
                [0.03e-5, 0.], [0.04e-5, 0.01e-5]], [[0., 0.]]]
            self.entropy_interaction = [
                [[0., 0.], [5.78, 5.78], [7.67, 7.67]], [[1.69, 1.69], [0., 0.]], [[0., 0.]]]

            # Published values are on a 4-oxygen (1-cation) basis
            for interaction in [self.energy_interaction, self.volume_interaction, self.entropy_interaction]:
                for i in range(len(interaction)):
                    for j in range(len(interaction[i])):
                        for k in range(len(interaction[i][j])):
                            interaction[i][j][k] *= 3.

            burnman.SolidSolution.__init__(self, molar_fractions)

    g5 = mg_fe_ca_garnet_Ganguly()
    g5_excess_gibbs = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1. - c, 0., c, 0.]
        g5.set_composition(molar_fractions)
        g5.set_state(1.e5, 298.15)
        g5_excess_gibbs[i] = g5.excess_gibbs

    plt.plot(comp, g5_excess_gibbs / 3., 'r-', linewidth=1.,
             label='Py-Gr excess gibbs (J/cation-mole)')

    plt.title("Asymmetric py-gr join (Ganguly et al., 1996; Figure 5)")
    plt.ylabel("Excess gibbs free energy of solution (J/cation-mol)")
    plt.xlabel("XCa")
    plt.legend(loc='lower left')
    plt.show()
