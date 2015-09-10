# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

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
        def __init__(self, molar_fractions=None):
            self.name='Ideal pyrope-almandine garnet'
            self.type='ideal'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g1=mg_fe_garnet()

    """
    Initialisation can optionally include setting the composition of the solid solution, i.e.
    """
    g1=mg_fe_garnet([0.0, 1.0])
    g1.set_state(P,T)
    alm_gibbs = g1.gibbs

    """
    Alternatively, the composition can be set after initialisation
    """
    g1.set_composition([1.0, 0.0])
    g1.set_state(P,T)
    py_gibbs = g1.gibbs

    comp = np.linspace(0.001, 0.999, 100)
    g1_gibbs = np.empty_like(comp)
    g1_excess_gibbs = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, c]
        g1.set_composition(molar_fractions)
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
        def __init__(self, molar_fractions=None):
            self.name='Symmetric pyrope-almandine garnet'
            self.type='symmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12']]
            self.enthalpy_interaction=[[2.5e3]]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g2=mg_fe_garnet_2()
    g2_excess_gibbs = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, c]
        g2.set_composition(molar_fractions)
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
        def __init__(self, molar_fractions=None):
            self.name='Symmetric pyrope-almandine-majorite garnet'
            self.type='symmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.maj(), '[Mg]3[Mg1/2Si1/2]2Si3O12']]
            self.enthalpy_interaction=[[2.5e3,0.0e3],[10.0e3]]
            self.entropy_interaction=[[0.0e3,0.0e3],[0.0e3]]
            self.volume_interaction=[[0.0e3,0.0e3],[0.0e3]]

            burnman.SolidSolution.__init__(self, molar_fractions)

    g3=hp_mg_fe_garnet()
    g3_configurational_entropy = np.empty_like(comp)
    g3_excess_entropy = np.empty_like(comp)
    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, 0., c]
        g3.set_composition(molar_fractions)
        g3.set_state(P,T)
        g3_configurational_entropy[i] = g3.solution_model._configurational_entropy( molar_fractions )
        g3_excess_entropy[i]= -np.dot(molar_fractions,g3.solution_model._ideal_excess_partial_gibbs(T, molar_fractions))/T

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
        def __init__(self, molar_fractions=None):
            self.name='Asymmetric pyrope-almandine-grossular garnet'
            self.type='asymmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12']]
            self.alphas=[1.0, 1.0, 2.7]
            self.enthalpy_interaction=[[2.5e3,30.1e3],[1.0e3]]
            self.volume_interaction=[[0.,0.169e-5],[0.122e-5]]
            burnman.SolidSolution.__init__(self, molar_fractions)


    g4=mg_fe_ca_garnet()
    g4_excess_gibbs_400 = np.empty_like(comp)
    g4_excess_gibbs_800 = np.empty_like(comp)
    g4_excess_gibbs_1200 = np.empty_like(comp)

    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, 0., c]
        g4.set_composition(molar_fractions)
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


    '''
    The subregular solution model (Helffrich and Wood, 1989)
    provides a more flexible way of constructing an asymmetric
    model
    '''
    class mg_fe_ca_garnet_Ganguly(burnman.SolidSolution):
        def __init__(self, molar_fractions=None):
            self.name='Subregular pyrope-almandine-grossular garnet'
            self.type='subregular'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],[minerals.HP_2011_ds62.alm(), '[Fe]3[Al]2Si3O12'],[minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12'], [minerals.HP_2011_ds62.spss(), '[Mn]3[Al]2Si3O12']]
            self.enthalpy_interaction=[[[2117., 695.], [9834., 21627.], [12083., 12083.]],[[6773., 873.],[539., 539.]],[[0., 0.]]]
            self.volume_interaction=[[[0.07e-5, 0.], [0.058e-5, 0.012e-5], [0.04e-5, 0.03e-5]],[[0.03e-5, 0.],[0.04e-5, 0.01e-5]],[[0., 0.]]]
            self.entropy_interaction=[[[0., 0.], [5.78, 5.78], [7.67, 7.67]],[[1.69, 1.69],[0., 0.]],[[0., 0.]]]
            
            # Published values are on a 4-oxygen (1-cation) basis
            for interaction in [self.enthalpy_interaction, self.volume_interaction, self.entropy_interaction]:
                for i in range(len(interaction)):
                    for j in range(len(interaction[i])):
                        for k in range(len(interaction[i][j])):
                            interaction[i][j][k]*=3.

            burnman.SolidSolution.__init__(self, molar_fractions)

    g5=mg_fe_ca_garnet_Ganguly()
    g5_excess_enthalpy = np.empty_like(comp)

    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, 0., c, 0.]
        g5.set_composition(molar_fractions)
        g5.set_state(1.e5,298.15)
        g5_excess_enthalpy[i] = g5.excess_enthalpy

    plt.plot( comp, g5_excess_enthalpy/3., 'r-', linewidth=1., label='Py-Gr excess enthalpy (J/cation-mole)')
    plt.title("Asymmetric py-gr join (Ganguly et al., 1996; Figure 5)")
    plt.ylabel("Excess enthalpy of solution (J/cation-mol)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='lower left')
    plt.show()


    '''
    Finally, the full subregular model is extremely flexible.
    The compositional variability in excess properties at any given 
    P and T is given by (Helffrich and Wood, 1989), as in the
    standard subregular model.
    The difference is in the construction of the intermediates. 
    Instead of defining excess properties based on the endmembers and
    constant H_ex, S_ex and V_ex, each binary interaction term is 
    described via an intermediate compound
    '''


    # First let's define a simple symmetric garnet for comparison
    class symmetric_garnet(burnman.SolidSolution):
        def __init__(self, molar_fractions=None):
            self.name='Symmetric pyrope-grossular garnet'
            self.type='symmetric'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],
                               [minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12']]
            self.enthalpy_interaction=[[H_ex]]
            self.entropy_interaction=[[S_ex]]
            self.volume_interaction=[[V_ex]]

            burnman.SolidSolution.__init__(self, molar_fractions)


    # Now let's define a full subregular garnet model
    # First, we need to create an intermediate compound
    # We can do this by using experimental data augmented by some reasonable heuristics

    # Here are the excess properties we want to use
    H_ex = 15000.*3. # J/mol; making a symmetric version of Ganguly et al., 1996.
    S_ex = 5.78*3. # m^3/mol; making a symmetric version of Ganguly et al., 1996.
    V_ex = 4.e-6 # m^3/mol, see Du et al, 2015


    # Now we define properties relative to the binary
    py = minerals.HP_2011_ds62.py()
    gr = minerals.HP_2011_ds62.gr()

    V_pygr = (py.params['V_0'] + gr.params['V_0'])*0.5 + V_ex/4.

    # Heuristics
    KV_py = py.params['K_0']*py.params['V_0']
    KV_gr = gr.params['K_0']*gr.params['V_0']
    
    aK_py = py.params['K_0']*py.params['a_0']
    aK_gr = gr.params['K_0']*gr.params['a_0']
    
    KV_pygr = 0.5*(KV_py + KV_gr)
    aK_pygr = 0.5*(aK_py + aK_gr)
    
    K_pygr = KV_pygr / V_pygr
    a_pygr = aK_pygr / K_pygr
    
    Kprime_0 = V_pygr*2./(py.params['V_0']/(py.params['Kprime_0'] + 1.) \
                           + py.params['V_0']/(py.params['Kprime_0'] + 1.))
    
    Sconf = 2.*burnman.constants.gas_constant*0.5*3.*np.log(0.5)
    S_pygr = (py.params['S_0'] + gr.params['S_0'])*0.5 - Sconf + S_ex/4.
    H_pygr = (py.params['H_0'] + gr.params['H_0'])*0.5 + H_ex/4.

    # Cp_scaling = ((py.params['S_0'] + gr.params['S_0'])*0.5 + S_ex/4)/((py.params['S_0'] + gr.params['S_0'])*0.5)
    # Haselton and Westrum (1980) show that the excess entropy is primarily a result of a low temperature spike in Cp
    # Therefore, at >=298.15 K, Cp is well approximated by a linear combination of pyrope and grossular 
    Cp_scaling = 1. 

    Cp_pygr = [(py.params['Cp'][0] + gr.params['Cp'][0])*0.5*Cp_scaling,
               (py.params['Cp'][1] + gr.params['Cp'][1])*0.5*Cp_scaling,
               (py.params['Cp'][2] + gr.params['Cp'][2])*0.5*Cp_scaling,
               (py.params['Cp'][3] + gr.params['Cp'][3])*0.5*Cp_scaling] 
  
    # Overloaded heuristics
    K_pygr = 162.e9
    Kprime_pygr = 5.0
    a_pygr = 2.3e-5

    # Creating the intermediate endmember is just like creating any other endmember
    from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass
    atomic_masses=read_masses()

    class pygr (burnman.Mineral):
        def __init__(self):
            formula='Mg1.5Ca1.5Al2.0Si3.0O12.0'
            formula = dictionarize_formula(formula)
            self.params = {
                'name': 'py-gr intermediate',
                'formula': formula,
                'equation_of_state': 'hp_tmt',
                'H_0': H_pygr,
                'S_0': S_pygr ,
                'V_0': V_pygr,
                'Cp': Cp_pygr ,
                'a_0': a_pygr, 
                'K_0': K_pygr,
                'Kprime_0': Kprime_pygr ,
                'Kdprime_0': -Kprime_pygr/K_pygr ,
                'n': sum(formula.values()),
                'molar_mass': formula_mass(formula, atomic_masses)}
            burnman.Mineral.__init__(self)

    # Finally, here's the solid solution class
    class full_garnet(burnman.SolidSolution):
        def __init__(self, molar_fractions=None):
            self.name='Subregular pyrope-almandine-grossular garnet'
            self.type='full_subregular'
            self.endmembers = [[minerals.HP_2011_ds62.py(), '[Mg]3[Al]2Si3O12'],
                               [minerals.HP_2011_ds62.gr(), '[Ca]3[Al]2Si3O12']]
            self.intermediates=[[[pygr(), pygr()]]]

            burnman.SolidSolution.__init__(self, molar_fractions)

    # Now, let's see what the model looks like in practise
    g6=full_garnet()
    g7=symmetric_garnet()

    pressures = np.linspace(1.e5, 25.e9, 26)
    py_volumes = np.empty_like(pressures)
    gr_volumes = np.empty_like(pressures)
    g6_volumes = np.empty_like(pressures)
    g7_volumes = np.empty_like(pressures)
    g6_excess_volumes = np.empty_like(pressures)
    g7_excess_volumes = np.empty_like(pressures)

    g6_excess_gibbs = np.empty_like(pressures)
    g7_excess_gibbs = np.empty_like(pressures)

    T = 298.15
    molar_fractions=[0.5, 0.5]
    for i,P in enumerate(pressures):
        py.set_state(P,T)
        py_volumes[i] = py.V
        gr.set_state(P,T)
        gr_volumes[i] = gr.V

        g6.set_composition(molar_fractions)
        g6.set_state(P,T)
        g6_volumes[i] = g6.V
        g6_excess_volumes[i] = g6.excess_volume
        g6_excess_gibbs[i] = g6.excess_gibbs

        g7.set_composition(molar_fractions)
        g7.set_state(P,T)
        g7_volumes[i] = g7.V
        g7_excess_volumes[i] = g7.excess_volume
        g7_excess_gibbs[i] = g7.excess_gibbs


    plt.plot( pressures/1.e9, py_volumes*1.e6, 'r-', linewidth=1., label='Py volume')
    plt.plot( pressures/1.e9, gr_volumes*1.e6, 'r-', linewidth=1., label='Gr volume')
    plt.plot( pressures/1.e9, g6_volumes*1.e6, 'g-', linewidth=1., label='Py50Gr50 (full)')
    plt.plot( pressures/1.e9, g7_volumes*1.e6, 'b-', linewidth=1., label='Py50Gr50 (simple)')
    plt.title("Room temperature py-gr equations of state")
    plt.ylabel("Volume (cm^3/mole)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc='lower left')
    plt.show()

    plt.plot( pressures/1.e9, g6_excess_volumes*1.e6, 'g-', linewidth=1., label='Py50Gr50 (full)')
    plt.plot( pressures/1.e9, g7_excess_volumes*1.e6, 'b-', linewidth=1., label='Py50Gr50 (simple)')
    plt.title("Py-gr volume excesses")
    plt.ylabel("Excess volume (cm^3/mole)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc='lower left')
    plt.show()


    plt.plot( pressures/1.e9, g6_excess_gibbs*1.e-3, 'g-', linewidth=1., label='Py50Gr50 (full)')
    plt.plot( pressures/1.e9, g7_excess_gibbs*1.e-3, 'b-', linewidth=1., label='Py50Gr50 (simple)')
    plt.title("Py-gr gibbs excesses")
    plt.ylabel("Excess gibbs (kJ/mole)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc='lower left')
    plt.show()

    dT = 1.
    g6_K_T = np.empty_like(comp)
    g7_K_T = np.empty_like(comp)

    g6_Cp = np.empty_like(comp)
    g7_Cp = np.empty_like(comp)

    g6_a = np.empty_like(comp)
    g7_a = np.empty_like(comp)

    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, c]
        g6.set_composition(molar_fractions)
        g6.set_state(1.e5,298.15)
        g6_K_T[i] = g6.K_T
        g6_a[i] = g6.alpha
        g6_Cp[i] = g6.C_p

        g7.set_composition(molar_fractions)
        g7.set_state(1.e5,298.15)
        g7_K_T[i] = g7.K_T
        g7_a[i] = g7.alpha
        g7_Cp[i] = g7.C_p

    plt.plot( comp, g6_K_T, 'g-', linewidth=1., label='Py-Gr (full)')
    plt.plot( comp, g7_K_T, 'b-', linewidth=1., label='Py-Gr (simple)')
    plt.title("Bulk modulus along the py-gr join")
    plt.ylabel("Bulk modulus (GPa)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='lower left')
    plt.show()

    plt.plot( comp, g6_Cp, 'g-', linewidth=1., label='Py-Gr (full)')
    plt.plot( comp, g7_Cp, 'b-', linewidth=1., label='Py-Gr (simple)')
    plt.title("Isobaric heat capacity along the py-gr join")
    plt.ylabel("Cp")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='lower left')
    plt.show()

    plt.plot( comp, g6_a, 'g-', linewidth=1., label='Py-Gr (full)')
    plt.plot( comp, g7_a, 'b-', linewidth=1., label='Py-Gr (simple)')
    plt.title("Thermal expansion along the py-gr join")
    plt.ylabel("alpha (/K)")
    plt.xlabel("Pyrope fraction")
    plt.legend(loc='lower left')
    plt.show()
