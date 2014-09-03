# Benchmarks for the solid solution class
import burnman
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

opx_red = orthopyroxene_red()
opx_red.set_method('slb3')
opx_blue = orthopyroxene_blue()
opx_blue.set_method('slb3')

comp = np.linspace(0, 1.0, 100)
entropy_red = np.empty_like(comp)
entropy_blue = np.empty_like(comp)

for i,c in enumerate(comp):
   opx_red.set_composition( np.array([1.0-c, c]) )
   opx_red.set_state( 0.0, 0.0 )
   entropy_red[i] = opx_red.solution_model.configurational_entropy( [1.0-c, c] )

   opx_blue.set_composition( np.array([1.0-c, c]) )
   opx_blue.set_state( 0.0, 0.0 )
   entropy_blue[i] = opx_blue.solution_model.configurational_entropy( [1.0-c, c] )


#fig1 = mpimg.imread('configurational_entropy.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
#plt.imshow(fig1, extent=[0.0, 1.0,0.,17.0], aspect='auto')

plt.plot( comp, entropy_red, 'r--', linewidth=3.)
plt.plot( comp, entropy_blue, 'b--', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(0.,17.0)
plt.ylabel("Configurational enthalpy of solution")
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


#fig1 = mpimg.imread('dicats.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
#plt.imshow(fig1, extent=[0.0, 1.0,-2.,8.0], aspect='auto')

plt.plot( comp, gibbs/1000., 'b--', linewidth=3.)
plt.xlim(0.0,1.0)
plt.ylim(-2.,8.0)
plt.ylabel("Excess enthalpy of solution")
plt.xlabel("cats fraction")
plt.show()

# Excess volume of solution
