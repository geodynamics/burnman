# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

example_solid_solution
-------------------
    
This example shows how to create different solid solution models and output
thermodynamic and thermoelastic quantities.

There are three main types of solid solution currently implemented in 
BurnMan:

1. Ideal solid solutions
2. Symmmetric solid solutions
3. Asymmetric solid solutions

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

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":
    '''
    First, let's pick a starting pressure and temperature
    '''
    P=1.e5
    T=1000.

    """
    Here we create an ideal solid solution for aluminous Mg,Fe garnets.
    Further endmembers can be added simply by adding more endmember lists.
    """ 
    class mg_fe_garnet(burnman.SolidSolution):
        def __init__(self):
            self.name='Ideal pyrope-almandine garnet'
            endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]
            burnman.SolidSolution.__init__(self, endmembers, burnman.solutionmodel.IdealSolution(endmembers) )

    g1=mg_fe_garnet()

    comp = np.linspace(0.001, 0.999, 100)

    g1.set_composition([0.0, 1.0])
    g1.set_state(P,T)
    alm_gibbs = g1.gibbs
    g1.set_composition([1.0, 0.0])
    g1.set_state(P,T)
    py_gibbs = g1.gibbs

    g1_gibbs = np.empty_like(comp)
    g1_excess_gibbs = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, c]
        g1.set_composition(molar_fraction)
        g1.set_state(P,T)
        g1_gibbs[i] = g1.gibbs
        g1_excess_gibbs[i] = g1.excess_gibbs

    plt.plot( [0., 1.], [py_gibbs/1000., alm_gibbs/1000.], 'b-', linewidth=1., label='Mechanical mixing')
    plt.plot( comp, g1_gibbs/1000., 'r-', linewidth=1., label='Ideal solid solution')
    plt.title("Pyrope-almandine join")
    plt.ylabel("Gibbs free energy of solution (kJ/mol)")
    plt.xlabel("Pyrope fraction")
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
        def __init__(self):
            self.name='Symmetric pyrope-almandine garnet'
            endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]
            enthalpy_interaction=[[2.5e3]]
            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.SymmetricRegularSolution(endmembers, enthalpy_interaction) )

    g2=mg_fe_garnet_2()
    g2_excess_gibbs = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, c]
        g2.set_composition(molar_fraction)
        g2.set_state(P,T)
        g2_excess_gibbs[i] = g2.excess_gibbs

    plt.plot( comp, g1_excess_gibbs, 'r-', linewidth=1., label='Ideal solution')
    plt.plot( comp, g2_excess_gibbs, 'g-', linewidth=1., label='Symmetric solution, 2.5 kJ/mol')
    plt.title("Pyrope-almandine join (model comparison)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='upper left')
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
        def __init__(self):
            self.name='Symmetric pyrope-almandine-majorite garnet'
            endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.maj(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
            enthalpy_interaction=[[2.5e3,0.0e3],[10.0e3]]
            entropy_interaction=[[0.0e3,0.0e3],[0.0e3]]
            volume_interaction=[[0.0e3,0.0e3],[0.0e3]]

            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.SymmetricRegularSolution(endmembers, enthalpy_interaction, volume_interaction, entropy_interaction) )

    g3=hp_mg_fe_garnet()
    g3_configurational_entropy = np.empty_like(comp)
    g3_excess_entropy = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, 0., c]
        g3.set_composition(molar_fraction)
        g3.set_state(P,T)
        g3_configurational_entropy[i] = g3.solution_model._configurational_entropy( molar_fraction )
        g3_excess_entropy[i]= -np.dot(molar_fraction,g3.solution_model._ideal_excess_partial_gibbs(T, molar_fraction))/T

    plt.plot( comp, g3_configurational_entropy, 'g-', linewidth=1., label='Configurational entropy')
    plt.plot( comp, g3_excess_entropy, 'r-', linewidth=1., label='Excess entropy')
    plt.title("Pyrope-majorite join")
    plt.ylabel("Entropy (J/K/mol)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='upper left')
    plt.show()

    '''
    The addition of grossular garnet (Ca-Al garnet) illustrates the use of asymmetric solution models
    '''
    class mg_fe_ca_garnet(burnman.SolidSolution):
        def __init__(self):
            self.name='Asymmetric pyrope-almandine-grossular garnet'
            endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12']]
            alphas=[1.0, 1.0, 2.7]
            enthalpy_interaction=[[2.5e3,30.1e3],[1.0e3]]
            volume_interaction=[[0.,0.169e-5],[0.122e-5]]
            burnman.SolidSolution.__init__(self, endmembers, \
                          burnman.solutionmodel.AsymmetricRegularSolution(endmembers, alphas, enthalpy_interaction, volume_interaction) )


    g4=mg_fe_ca_garnet()
    g4_excess_gibbs_400 = np.empty_like(comp)
    g4_excess_gibbs_800 = np.empty_like(comp)
    g4_excess_gibbs_1200 = np.empty_like(comp)

    for i,c in enumerate(comp):
        molar_fraction=[1.0-c, 0., c]
        g4.set_composition(molar_fraction)
        g4.set_state(P,400.)
        g4_excess_gibbs_400[i] = g4.excess_gibbs
        g4.set_state(P,800.)
        g4_excess_gibbs_800[i] = g4.excess_gibbs
        g4.set_state(P,1200.)
        g4_excess_gibbs_1200[i] = g4.excess_gibbs

    plt.plot( comp, g4_excess_gibbs_400, 'r-', linewidth=1., label='400 K')
    plt.plot( comp, g4_excess_gibbs_800, 'g-', linewidth=1., label='800 K')
    plt.plot( comp, g4_excess_gibbs_1200, 'b-', linewidth=1., label='1200 K')
    plt.title("Pyrope-grossular join (asymmetric model)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='lower left')
    plt.show()

