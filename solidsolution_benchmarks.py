# Benchmarks for the solid solution class
import burnman
from burnman import minerals
from burnman.processchemistry import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

atomic_masses=read_masses()




'''
Solvus shapes (a proxy for Gibbs free energy checking
'''

# van Laar parameter
# Figure 2a of Holland and Powell, 2003

# Temperature dependence
# Figure 2b of Holland and Powell, 2003

# A specific solvus example: sanidine-high albite
# Includes asymmetry and pressure, temperature dependence
# Figure 3 of Holland and Powell, 2003


'''
Excess properties
'''
# Configurational entropy
# Navrotsky and Kleppa, 1967
class o_d_spinel(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[minerals.HP_2011.spinel(), '[Mg][Al]2O4'],[minerals.HP_2011.spinel(), '[Al][Mg1/2Al1/2]2O4']]

        # Interaction parameters
        enthalpy_interaction=[[0.0]]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )

comp = np.linspace(0.001, 0.999, 100)
sp=o_d_spinel()
sp.set_method('mtait')
sp_entropies = np.empty_like(comp)
sp_entropies_NK1967= np.empty_like(comp)
for i,c in enumerate(comp):
        molar_fractions=[1.0-c, c]
        sp.set_composition( np.array(molar_fractions) )
        sp.set_state( 1e5, 298.15 )
        sp_entropies[i] = sp.solution_model.configurational_entropy( molar_fractions )
        sp_entropies_NK1967[i] = -8.3145*(c*np.log(c) + (1.-c)*np.log(1.-c) + c*np.log(c/2.) + (2.-c)*np.log(1.-c/2.)) # eq. 7 in Navrotsky and Kleppa, 1967.

#fig1 = mpimg.imread('configurational_entropy.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
#plt.imshow(fig1, extent=[0.0, 1.0,0.,17.0], aspect='auto')
plt.plot( comp, sp_entropies_NK1967, 'b-', linewidth=3.)
plt.plot( comp, sp_entropies, 'r--', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(0.,17.0)
plt.ylabel("Configurational entropy of solution (J/K/mol)")
plt.xlabel("fraction inverse spinel")
plt.show()

# Configurational entropy
# Figure 3b of Stixrude and Lithgow-Bertelloni, 2011
class enstatite (burnman.Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -2913.e3,
            'V_0': 62.68e-6,
            'K_0': 107.e9,
            'Kprime_0': 7.0 ,
            'Debye_0': 812.0,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4,
            'G_0': 77.e9,
            'Gprime_0':  1.5,
            'eta_s_0': 2.5 }

class mg_tschermaks_molecule (burnman.Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -3003.e3 ,
            'V_0': 59.15e-6,
            'K_0': 107.e9,
            'Kprime_0': 7.0,
            'Debye_0': 784.,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4,
            'G_0': 97.e9,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 }

class orthopyroxene_red(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[enstatite(), 'Mg[Mg][Si]SiO6'],[mg_tschermaks_molecule(), 'Mg[Al][Al]SiO6'] ]

        # Interaction parameters
        enthalpy_interaction=[[0.0]]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )
class orthopyroxene_blue(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[enstatite(), 'Mg[Mg]Si2O6'],[mg_tschermaks_molecule(), 'Mg[Al]AlSiO6'] ]

        # Interaction parameters
        enthalpy_interaction=[[0.0]]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )

class orthopyroxene_long_dashed(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[enstatite(), 'Mg[Mg]Si2O6'],[mg_tschermaks_molecule(), '[Mg1/2Al1/2]2AlSiO6'] ]

        # Interaction parameters
        enthalpy_interaction=[[10.0e3]]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )

class orthopyroxene_short_dashed(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[enstatite(), 'Mg[Mg][Si]2O6'],[mg_tschermaks_molecule(), 'Mg[Al][Al1/2Si1/2]2O6'] ]

        # Interaction parameters
        enthalpy_interaction=[[0.0]]

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.SymmetricRegularSolution(base_material, enthalpy_interaction) )

comp = np.linspace(0, 1.0, 100)
opx_models=[orthopyroxene_red(), orthopyroxene_blue(), orthopyroxene_long_dashed(), orthopyroxene_short_dashed()]
opx_entropies = [ np.empty_like(comp) for model in opx_models ]
for idx, model in enumerate(opx_models):
    model.set_method('slb3')

    for i,c in enumerate(comp):
        molar_fractions=[1.0-c, c]
        model.set_composition( np.array(molar_fractions) )
        model.set_state( 0.0, 0.0 )
        opx_entropies[idx][i] = model.solution_model.configurational_entropy( molar_fractions )



fig1 = mpimg.imread('configurational_entropy.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
plt.imshow(fig1, extent=[0.0, 1.0,0.,17.0], aspect='auto')
plt.plot( comp, opx_entropies[0], 'r--', linewidth=3.)
plt.plot( comp, opx_entropies[1], 'b--', linewidth=3.)
plt.plot( comp, opx_entropies[2], 'g--', linewidth=3.)
plt.plot( comp, opx_entropies[3], 'g-.', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(0.,17.0)
plt.ylabel("Configurational entropy of solution (J/K/mol)")
plt.xlabel("cats fraction")
plt.show()

# Excess volume of solution

# Excess enthalpy of solution
# Figure 5 of Stixrude and Lithgow-Bertelloni, 2011

class ca_tschermaks_molecule (burnman.Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -3120.0e3,
            'V_0': 63.57e-6,
            'K_0': 112.e9,
            'Kprime_0': 5.2 ,
            'Debye_0': 804.,
            'grueneisen_0': 0.78 ,
            'q_0': 1.5,
            'G_0': 76.e9,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.0 }

class diopside (burnman.Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -3030.0e3,
            'V_0': 66.04e-6,
            'K_0': 112.e9,
            'Kprime_0': 5.2 ,
            'Debye_0': 782.,
            'grueneisen_0': 0.96 ,
            'q_0': 1.5,
            'G_0': 67.e9,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.6 }


class clinopyroxene(burnman.SolidSolution):
    def __init__(self):
        # Name
        self.name='clinopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[diopside(), '[Ca][Mg][Si]2O6'],[ca_tschermaks_molecule(), '[Ca][Al][Si1/2Al1/2]2O6'] ]

        # Interaction parameters
        enthalpy_interaction=[[26.e3]]
        alphas = np.array( [1.0, 3.5] ) 

        burnman.SolidSolution.__init__(self, base_material, \
                          burnman.solutionmodel.AsymmetricRegularSolution(base_material, alphas, enthalpy_interaction) )



cpx = clinopyroxene()
cpx.set_method('slb3')

comp = np.linspace(0, 1.0, 100)
gibbs = np.empty_like(comp)

for i,c in enumerate(comp):
   cpx.set_composition( np.array([1.0-c, c]) )
   cpx.set_state( 0.0, 0.0 )
   gibbs[i] = cpx.excess_gibbs


fig1 = mpimg.imread('dicats.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
plt.imshow(fig1, extent=[0.0, 1.0,-2.,8.0], aspect='auto')

plt.plot( comp, gibbs/1000., 'b--', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(-2.,8.0)
plt.ylabel("Excess enthalpy of solution (kJ/mol)")
plt.xlabel("cats fraction")
plt.show()

# Check endmember excess Gibbs goes to zero... 

opx = orthopyroxene_long_dashed()
opx.set_method('slb3')

comp = np.linspace(0, 1.0, 100)
gibbs = np.empty_like(comp)

for i,c in enumerate(comp):
   opx.set_composition( np.array([1.0-c, c]) )
   opx.set_state( 0.0, 2000.0 )
   gibbs[i] = opx.excess_gibbs


plt.plot( comp, gibbs/1000., 'b--', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(-2.,8.0)
plt.ylabel("Excess enthalpy of solution (kJ/mol)")
plt.xlabel("mgts fraction")
plt.show()
