# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""

SLB_2011
^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2011 and references therein
Solid solutions from inv251010 of HeFESTo

"""

from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution
from burnman.solutionmodel import *
from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass

atomic_masses=read_masses()

'''
SOLID SOLUTIONS
'''

class c2c_pyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='C2/c pyroxene'

        # Endmembers (C2/c is symmetric)
        base_material = [[hp_clinoenstatite(), '[Mg]2Si2O6'],[hp_clinoferrosilite(), '[Fe]2Si2O6']]

        SolidSolution.__init__(self, base_material, IdealSolution(base_material))



class ca_ferrite_structured_phase(SolidSolution):
    def __init__(self):
        # Name
        self.name='calcium ferrite structured phase'

        # Endmembers (CF is symmetric)
        base_material = [[mg_ca_ferrite(), '[Mg]Al[Al]O4'],[fe_ca_ferrite(), '[Fe]Al[Al]O4'],[na_ca_ferrite(), '[Na]Al[Si]O4']]

        SolidSolution.__init__(self, base_material, IdealSolution(base_material))

class clinopyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='clinopyroxene'

        # Endmembers (cpx is symmetric)
        base_material = [[diopside(), '[Ca][Mg][Si]2O6'],[hedenbergite(), '[Ca][Fe][Si]2O6'],[clinoenstatite(), '[Mg][Mg][Si]2O6'],[ca_tschermaks(), '[Ca][Al][Si1/2Al1/2]2O6'],[jadeite(), '[Na][Al][Si]2O6']]

        # Interaction parameters
        enthalpy_interaction=[[0., 24.74e3, 26.e3, 24.3e3],[24.74e3, 0., 0.e3], [60.53136e3, 0.0], [10.e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )

class garnet(SolidSolution):
    def __init__(self):
        # Name
        self.name='garnet'

        # Endmembers (garnet is symmetric)
        base_material = [[pyrope(), '[Mg]3[Al][Al]Si3O12'],[almandine(), '[Fe]3[Al][Al]Si3O12'],[grossular(), '[Ca]3[Al][Al]Si3O12'],[mg_majorite(), '[Mg]3[Mg][Si]Si3O12'],[jd_majorite(), '[Na2/3Al1/3]3[Al][Si]Si3O12']]
        # Interaction parameters
        enthalply_interaction=[[0.0, 30.e3, 21.20278e3, 0.0],[0.0,0.0,0.0],[57.77596e3, 0.0],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )


class akimotoite(SolidSolution):
    def __init__(self):
        # Name
        self.name='akimotoite/ilmenite'

        # Endmembers (ilmenite/akimotoite is symmetric)
        base_material = [[mg_akimotoite(), '[Mg][Si]O3'],[fe_akimotoite(), '[Fe][Si]O3'],[corundum(), '[Al][Al]O3']]
        # Interaction parameters
        enthalpy_interaction=[[0.0, 66.e3],[66.e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )

class ferropericlase(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesiowustite/ferropericlase'

        # Endmembers (ferropericlase is symmetric)
        base_material = [[periclase(), '[Mg]O'],[wuestite(), '[Fe]O']]
        # Interaction parameters
        enthalpy_interaction=[[13.e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )

class mg_fe_olivine(SolidSolution):
    def __init__(self):
        # Name
        self.name='olivine'

        # Endmembers (olivine is symmetric)
        base_material = [[forsterite(), '[Mg]2SiO4'],[fayalite(), '[Fe]2SiO4']]
        # Interaction parameters
        enthalpy_interaction=[[7.81322e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )

class orthopyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (orthopyroxene is symmetric)
        base_material = [[enstatite(), '[Mg][Mg][Si]SiO6'],[ferrosilite(), '[Fe][Fe][Si]SiO6'],[mg_tschermaks(), '[Mg][Al][Al]SiO6'],[ortho_diopside(), '[Ca][Mg][Si]SiO6']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 0.0, 32.11352e3],[0.0, 0.0],[48.35316e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class plagioclase(SolidSolution):
    def __init__(self):
        # Name
        self.name='plagioclase'

        # Endmembers (plagioclase is symmetric)
        base_material = [[anorthite(), '[Ca][Al]2Si2O8'],[albite(), '[Na][Al1/2Si1/2]2Si2O8']]
        # Interaction parameters
        enthalpy_interaction=[[26.0e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class post_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='post-perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_post_perovskite(), '[Mg][Si]O3'],[fe_post_perovskite(), '[Fe][Si]O3'],[al_post_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 60.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class mg_fe_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesium silicate perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_perovskite(), '[Mg][Si]O3'],[fe_perovskite(), '[Fe][Si]O3'],[al_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 116.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class mg_fe_ringwoodite(SolidSolution):
    def __init__(self):
        # Name
        self.name='ringwoodite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_ringwoodite(), '[Mg]2SiO4'],[fe_ringwoodite(), '[Fe]2SiO4']]

        # Interaction parameters
        enthalpy_interaction=[[9.34084e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class mg_fe_aluminous_spinel(SolidSolution):
    def __init__(self):
        # Name
        self.name='spinel-hercynite binary, fixed order'

        # Endmembers (post perovskite is symmetric)
        base_material = [[spinel(), '[Mg3/4Al1/4]4[Al7/8Mg1/8]8O16'],[hercynite(), '[Fe3/4Al1/4]4[Al7/8Fe1/8]8O16']]

        # Interaction parameters
        enthalpy_interaction=[[5.87646e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction))

class mg_fe_wadsleyite(SolidSolution):
    def __init__(self):
        # Name
        self.name='wadsleyite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_wadsleyite(), '[Mg]2SiO4'],[fe_wadsleyite(), '[Fe]2SiO4']]

        # Interaction parameters
        enthalpy_interaction=[[16.74718e3]]

        SolidSolution.__init__(self, base_material, SymmetricRegularSolution(base_material, enthalpy_interaction) )

'''
ENDMEMBERS
'''

class anorthite (Mineral):
    def __init__(self):
       formula='CaAl2Si2O8'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Anorthite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -4015000.0 ,
            'V_0': 0.00010061 ,
            'K_0': 84000000000.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 752.0 ,
            'grueneisen_0': 0.39 ,
            'q_0': 1.0 ,
            'G_0': 40000000000.0 ,
            'Gprime_0': 1.1 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 3000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class albite (Mineral):
    def __init__(self):
       formula='NaAlSi3O8'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Albite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3719000.0 ,
            'V_0': 0.00010045 ,
            'K_0': 60000000000.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 716.0 ,
            'grueneisen_0': 0.57 ,
            'q_0': 1.0 ,
            'G_0': 36000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 5000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class spinel (Mineral):
    def __init__(self):
       formula='Mg4Al8O16'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Spinel',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -8668000.0 ,
            'V_0': 0.00015905 ,
            'K_0': 1.97e+11 ,
            'Kprime_0': 5.7 ,
            'Debye_0': 843.0 ,
            'grueneisen_0': 1.02 ,
            'q_0': 2.7 ,
            'G_0': 1.08e+11 ,
            'Gprime_0': 0.4 ,
            'eta_s_0': 2.7 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 32000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 33.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.6 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 0.6 }

class hercynite (Mineral):
    def __init__(self):
       formula='Fe4Al8O16'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Hercynite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -7324000.0 ,
            'V_0': 0.00016337 ,
            'K_0': 2.09e+11 ,
            'Kprime_0': 5.7 ,
            'Debye_0': 763.0 ,
            'grueneisen_0': 1.22 ,
            'q_0': 2.7 ,
            'G_0': 84000000000.0 ,
            'Gprime_0': 0.4 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 35000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 32.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 13000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class forsterite (Mineral):
    def __init__(self):
       formula='Mg2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Forsterite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2055000.0 ,
            'V_0': 4.36e-05 ,
            'K_0': 1.28e+11 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 809.0 ,
            'grueneisen_0': 0.99 ,
            'q_0': 2.1 ,
            'G_0': 82000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 1.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 0.2 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.1 }

class fayalite (Mineral):
    def __init__(self):
       formula='Fe2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fayalite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1371000.0 ,
            'V_0': 4.629e-05 ,
            'K_0': 1.35e+11 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 619.0 ,
            'grueneisen_0': 1.06 ,
            'q_0': 3.6 ,
            'G_0': 51000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 0.6 }

class mg_wadsleyite (Mineral):
    def __init__(self):
       formula='Mg2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Wadsleyite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2028000.0 ,
            'V_0': 4.052e-05 ,
            'K_0': 1.69e+11 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 844.0 ,
            'grueneisen_0': 1.21 ,
            'q_0': 2.0 ,
            'G_0': 1.12e+11 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.09 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.2 ,
            'err_eta_s_0': 0.4 }

class fe_wadsleyite (Mineral):
    def __init__(self):
       formula='Fe2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Wadsleyite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1365000.0 ,
            'V_0': 4.28e-05 ,
            'K_0': 1.69e+11 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 665.0 ,
            'grueneisen_0': 1.21 ,
            'q_0': 2.0 ,
            'G_0': 72000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 7000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 13000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 21.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 12000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class mg_ringwoodite (Mineral):
    def __init__(self):
       formula='Mg2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Ringwoodite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2017000.0 ,
            'V_0': 3.949e-05 ,
            'K_0': 1.85e+11 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 878.0 ,
            'grueneisen_0': 1.11 ,
            'q_0': 2.4 ,
            'G_0': 1.23e+11 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.1 ,
            'err_q_0': 0.4 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.5 }

class fe_ringwoodite (Mineral):
    def __init__(self):
       formula='Fe2SiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Ringwoodite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1363000.0 ,
            'V_0': 4.186e-05 ,
            'K_0': 2.13e+11 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 679.0 ,
            'grueneisen_0': 1.27 ,
            'q_0': 2.4 ,
            'G_0': 92000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 7000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.23 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class enstatite (Mineral):
    def __init__(self):
       formula='Mg2Si2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Enstatite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2913000.0 ,
            'V_0': 6.268e-05 ,
            'K_0': 1.07e+11 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 812.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 77000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.4 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.4 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.1 }

class ferrosilite (Mineral):
    def __init__(self):
       formula='Fe2Si2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Ferrosilite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2226000.0 ,
            'V_0': 6.594e-05 ,
            'K_0': 1.01e+11 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 674.0 ,
            'grueneisen_0': 0.72 ,
            'q_0': 3.4 ,
            'G_0': 52000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4000000000.0 ,
            'err_K_prime_0': 0.5 ,
            'err_Debye_0': 10.0 ,
            'err_grueneisen_0': 0.08 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class mg_tschermaks (Mineral):
    def __init__(self):
       formula='MgAl2SiO6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Tschermaks',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3003000.0 ,
            'V_0': 5.914e-05 ,
            'K_0': 1.07e+11 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 784.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 97000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 9000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 24.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class ortho_diopside (Mineral):
    def __init__(self):
       formula='CaMgSi2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Ortho_Diopside',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3016000.0 ,
            'V_0': 6.805e-05 ,
            'K_0': 1.07e+11 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 745.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 60000000000.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class diopside (Mineral):
    def __init__(self):
       formula='CaMgSi2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Diopside',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3030000.0 ,
            'V_0': 6.604e-05 ,
            'K_0': 1.12e+11 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 782.0 ,
            'grueneisen_0': 0.96 ,
            'q_0': 1.5 ,
            'G_0': 67000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5000000000.0 ,
            'err_K_prime_0': 1.8 ,
            'err_Debye_0': 3.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 2.0 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class hedenbergite (Mineral):
    def __init__(self):
       formula='CaFeSi2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Hedenbergite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2677000.0 ,
            'V_0': 6.787e-05 ,
            'K_0': 1.19e+11 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 702.0 ,
            'grueneisen_0': 0.94 ,
            'q_0': 1.5 ,
            'G_0': 61000000000.0 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 45000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class clinoenstatite (Mineral):
    def __init__(self):
       formula='Mg2Si2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Clinoenstatite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2906000.0 ,
            'V_0': 6.25e-05 ,
            'K_0': 1.12e+11 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 805.0 ,
            'grueneisen_0': 0.96 ,
            'q_0': 1.5 ,
            'G_0': 81000000000.0 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 1.7 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 10.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class ca_tschermaks (Mineral):
    def __init__(self):
       formula='CaAl2SiO6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Ca_Tschermaks',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3120000.0 ,
            'V_0': 6.357e-05 ,
            'K_0': 1.12e+11 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 804.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 1.5 ,
            'G_0': 76000000000.0 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 5000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.0 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class jadeite (Mineral):
    def __init__(self):
       formula='NaAlSi2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Jadeite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2855000.0 ,
            'V_0': 6.051e-05 ,
            'K_0': 1.42e+11 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 821.0 ,
            'grueneisen_0': 0.9 ,
            'q_0': 0.4 ,
            'G_0': 85000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.2 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 12.0 ,
            'err_grueneisen_0': 0.08 ,
            'err_q_0': 1.4 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class hp_clinoenstatite (Mineral):
    def __init__(self):
       formula='Mg2Si2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'HP_Clinoenstatite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2905000.0 ,
            'V_0': 6.076e-05 ,
            'K_0': 1.16e+11 ,
            'Kprime_0': 6.2 ,
            'Debye_0': 824.0 ,
            'grueneisen_0': 1.12 ,
            'q_0': 0.2 ,
            'G_0': 88000000000.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.5 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.5 }

class hp_clinoferrosilite (Mineral):
    def __init__(self):
       formula='Fe2Si2O6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'HP_Clinoferrosilite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2222000.0 ,
            'V_0': 6.385e-05 ,
            'K_0': 1.16e+11 ,
            'Kprime_0': 6.2 ,
            'Debye_0': 692.0 ,
            'grueneisen_0': 1.12 ,
            'q_0': 0.2 ,
            'G_0': 71000000000.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 0.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 11.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class ca_perovskite (Mineral):
    def __init__(self):
       formula='CaSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Ca_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1463000.0 ,
            'V_0': 2.745e-05 ,
            'K_0': 2.36e+11 ,
            'Kprime_0': 3.9 ,
            'Debye_0': 796.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 0.9 ,
            'G_0': 1.57e+11 ,
            'Gprime_0': 2.2 ,
            'eta_s_0': 1.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 8000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 44.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.6 ,
            'err_G_0': 12000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class mg_akimotoite (Mineral):
    def __init__(self):
       formula='MgSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Akimotoite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1410000.0 ,
            'V_0': 2.635e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 5.6 ,
            'Debye_0': 934.0 ,
            'grueneisen_0': 1.19 ,
            'q_0': 2.3 ,
            'G_0': 1.32e+11 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4000000000.0 ,
            'err_K_prime_0': 0.8 ,
            'err_Debye_0': 12.0 ,
            'err_grueneisen_0': 0.13 ,
            'err_q_0': 0.8 ,
            'err_G_0': 8000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class fe_akimotoite (Mineral):
    def __init__(self):
       formula='FeSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Akimotoite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1068000.0 ,
            'V_0': 2.685e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 5.6 ,
            'Debye_0': 888.0 ,
            'grueneisen_0': 1.19 ,
            'q_0': 2.3 ,
            'G_0': 1.5e+11 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 3.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 21000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 120.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class corundum (Mineral):
    def __init__(self):
       formula='AlAlO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Corundum',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1582000.0 ,
            'V_0': 2.558e-05 ,
            'K_0': 2.53e+11 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 933.0 ,
            'grueneisen_0': 1.32 ,
            'q_0': 1.3 ,
            'G_0': 1.63e+11 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 3.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.2 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.2 }

class pyrope (Mineral):
    def __init__(self):
       formula='Mg3Al2Si3O12'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Pyrope',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -5936000.0 ,
            'V_0': 0.00011308 ,
            'K_0': 1.7e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 823.0 ,
            'grueneisen_0': 1.01 ,
            'q_0': 1.4 ,
            'G_0': 94000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 10000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 0.5 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.2 ,
            'err_eta_s_0': 0.3 }

class almandine (Mineral):
    def __init__(self):
       formula='Fe3Al2Si3O12'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Almandine',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -4935000.0 ,
            'V_0': 0.00011543 ,
            'K_0': 1.74e+11 ,
            'Kprime_0': 4.9 ,
            'Debye_0': 741.0 ,
            'grueneisen_0': 1.06 ,
            'q_0': 1.4 ,
            'G_0': 96000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 29000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class grossular (Mineral):
    def __init__(self):
       formula='Ca3Al2Si3O12'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Grossular',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -6278000.0 ,
            'V_0': 0.00012512 ,
            'K_0': 1.67e+11 ,
            'Kprime_0': 3.9 ,
            'Debye_0': 823.0 ,
            'grueneisen_0': 1.05 ,
            'q_0': 1.9 ,
            'G_0': 1.09e+11 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 11000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 0.2 ,
            'err_G_0': 4000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.1 }

class mg_majorite (Mineral):
    def __init__(self):
       formula='Mg4Si4O12'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Majorite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -5691000.0 ,
            'V_0': 0.00011432 ,
            'K_0': 1.65e+11 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 822.0 ,
            'grueneisen_0': 0.98 ,
            'q_0': 1.5 ,
            'G_0': 85000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 10000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3000000000.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 0.5 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.2 ,
            'err_eta_s_0': 0.3 }

class jd_majorite (Mineral):
    def __init__(self):
       formula='Na2Al2Si4O12'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Jd_Majorite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -5519000.0 ,
            'V_0': 0.00011094 ,
            'K_0': 1.77e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 896.0 ,
            'grueneisen_0': 1.01 ,
            'q_0': 1.4 ,
            'G_0': 1.25e+11 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 3.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 14000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 7000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 18.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 4000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class quartz (Mineral):
    def __init__(self):
       formula='SiO2'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Quartz',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -859000.0 ,
            'V_0': 2.367e-05 ,
            'K_0': 50000000000.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 816.0 ,
            'grueneisen_0': 0.0 ,
            'q_0': 1.0 ,
            'G_0': 45000000000.0 ,
            'Gprime_0': 1.0 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 31.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class coesite (Mineral):
    def __init__(self):
       formula='SiO2'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Coesite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -855000.0 ,
            'V_0': 2.066e-05 ,
            'K_0': 1.14e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 857.0 ,
            'grueneisen_0': 0.39 ,
            'q_0': 1.0 ,
            'G_0': 62000000000.0 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class stishovite (Mineral):
    def __init__(self):
       formula='SiO2'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Stishovite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -819000.0 ,
            'V_0': 1.402e-05 ,
            'K_0': 3.14e+11 ,
            'Kprime_0': 3.8 ,
            'Debye_0': 1108.0 ,
            'grueneisen_0': 1.37 ,
            'q_0': 2.8 ,
            'G_0': 2.2e+11 ,
            'Gprime_0': 1.9 ,
            'eta_s_0': 4.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 8000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.17 ,
            'err_q_0': 2.2 ,
            'err_G_0': 12000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class seifertite (Mineral):
    def __init__(self):
       formula='SiO2'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Seifertite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -794000.0 ,
            'V_0': 1.367e-05 ,
            'K_0': 3.28e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 1141.0 ,
            'grueneisen_0': 1.37 ,
            'q_0': 2.8 ,
            'G_0': 2.27e+11 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 5.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 16.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class mg_perovskite (Mineral):
    def __init__(self):
       formula='MgSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1368000.0 ,
            'V_0': 2.445e-05 ,
            'K_0': 2.51e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 905.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 1.73e+11 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 2.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.3 ,
            'err_G_0': 2000000000.0 ,
            'err_Gprime_0': 0.0 ,
            'err_eta_s_0': 0.3 }

class fe_perovskite (Mineral):
    def __init__(self):
       formula='FeSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1041000.0 ,
            'V_0': 2.549e-05 ,
            'K_0': 2.72e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 871.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 1.33e+11 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 6000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 40000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 26.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 40000000000.0 ,
            'err_Gprime_0': 0.0 ,
            'err_eta_s_0': 1.0 }

class al_perovskite (Mineral):
    def __init__(self):
       formula='AlAlO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Al_perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1534000.0 ,
            'V_0': 2.494e-05 ,
            'K_0': 2.58e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 886.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 1.71e+11 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 0.5 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.5 }

class mg_post_perovskite (Mineral):
    def __init__(self):
       formula='MgSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Post_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1348000.0 ,
            'V_0': 2.442e-05 ,
            'K_0': 2.31e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 855.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 1.1 ,
            'G_0': 1.5e+11 ,
            'Gprime_0': 2.0 ,
            'eta_s_0': 1.2 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 0.1 ,
            'err_G_0': 4000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.2 }

class fe_post_perovskite (Mineral):
    def __init__(self):
       formula='FeSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Post_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -982000.0 ,
            'V_0': 2.546e-05 ,
            'K_0': 2.31e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 782.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 1.1 ,
            'G_0': 1.29e+11 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 21000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 52.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class al_post_perovskite (Mineral):
    def __init__(self):
       formula='AlAlO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Al_Post_Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1378000.0 ,
            'V_0': 2.385e-05 ,
            'K_0': 2.49e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 762.0 ,
            'grueneisen_0': 1.65 ,
            'q_0': 1.1 ,
            'G_0': 92000000000.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 20000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.02 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.2 }

class periclase (Mineral):
    def __init__(self):
       formula='MgO'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Periclase',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -569000.0 ,
            'V_0': 1.124e-05 ,
            'K_0': 1.61e+11 ,
            'Kprime_0': 3.8 ,
            'Debye_0': 767.0 ,
            'grueneisen_0': 1.36 ,
            'q_0': 1.7 ,
            'G_0': 1.31e+11 ,
            'Gprime_0': 2.1 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 0.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.2 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 0.2 }

class wuestite (Mineral):
    def __init__(self):
       formula='FeO'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Wuestite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -242000.0 ,
            'V_0': 1.226e-05 ,
            'K_0': 1.79e+11 ,
            'Kprime_0': 4.9 ,
            'Debye_0': 454.0 ,
            'grueneisen_0': 1.53 ,
            'q_0': 1.7 ,
            'G_0': 59000000000.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': -0.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 21.0 ,
            'err_grueneisen_0': 0.13 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class mg_ca_ferrite (Mineral):
    def __init__(self):
       formula='MgAl2O4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Ca_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2122000.0 ,
            'V_0': 3.618e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 838.0 ,
            'grueneisen_0': 1.31 ,
            'q_0': 1.0 ,
            'G_0': 1.3e+11 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 16.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class fe_ca_ferrite (Mineral):
    def __init__(self):
       formula='FeAl2O4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Ca_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1790000.0 ,
            'V_0': 3.726e-05 ,
            'K_0': 2.11e+11 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 804.0 ,
            'grueneisen_0': 1.31 ,
            'q_0': 1.0 ,
            'G_0': 1.52e+11 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 3.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 25000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 69.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class na_ca_ferrite (Mineral):
    def __init__(self):
       formula='NaAlSiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Na_Ca_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1851000.0 ,
            'V_0': 3.627e-05 ,
            'K_0': 1.58e+11 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 812.0 ,
            'grueneisen_0': 1.17 ,
            'q_0': 1.0 ,
            'G_0': 1.21e+11 ,
            'Gprime_0': 2.1 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 11000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 51.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class kyanite (Mineral):
    def __init__(self):
       formula='Al2SiO5'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Kyanite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2446000.0 ,
            'V_0': 4.423e-05 ,
            'K_0': 1.6e+11 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 943.0 ,
            'grueneisen_0': 0.93 ,
            'q_0': 1.0 ,
            'G_0': 1.21e+11 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 3.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 0.0 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class nepheline (Mineral):
    def __init__(self):
       formula='NaAlSiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Nepheline',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1993000.0 ,
            'V_0': 5.467e-05 ,
            'K_0': 53000000000.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 701.0 ,
            'grueneisen_0': 0.69 ,
            'q_0': 1.0 ,
            'G_0': 31000000000.0 ,
            'Gprime_0': 1.3 ,
            'eta_s_0': 0.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3000.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1000000000.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1000000000.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

'''
Mineral aliases
'''

# Feldspars
ab = albite
an = anorthite

# LP Spinels
sp = spinel
hc = hercynite

# Olivine polymorphs
fo = forsterite
fa = fayalite
mgwa = mg_wadsleyite
fewa = fe_wadsleyite
mgri = mg_ringwoodite
feri = fe_ringwoodite

# Orthopyroxenes
en = enstatite
fs = ferrosilite
mgts = mg_tschermaks
odi = ortho_diopside

# Clinopyroxenes
di = diopside
he = hedenbergite
cen = clinoenstatite
cats = ca_tschermaks
jd = jadeite
mgc2 = hp_clinoenstatite
fec2 = hp_clinoferrosilite
hpcen = hp_clinoenstatite
hpcfs = hp_clinoferrosilite

# Perovskites
mgpv = mg_perovskite
mg_bridgmanite = mg_perovskite
fepv = fe_perovskite
fe_bridgmanite = fe_perovskite
alpv = al_perovskite
capv = ca_perovskite

# Ilmenite group
mgil = mg_akimotoite
feil = fe_akimotoite
co = corundum

# Garnet group
py = pyrope
al = almandine
gr = grossular
mgmj = mg_majorite
jdmj = jd_majorite

# Quartz polymorphs
qtz = quartz
coes = coesite
st = stishovite
seif = seifertite

# Post perovskites
mppv = mg_post_perovskite
fppv = fe_post_perovskite
appv = al_post_perovskite

# Magnesiowuestite
pe = periclase
wu = wuestite

# Calcium ferrite structured phases
mgcf = mg_ca_ferrite
fecf = fe_ca_ferrite
nacf = na_ca_ferrite

# Al2SiO5 polymorphs
ky = kyanite

# Nepheline group
neph = nepheline



# Solid solution aliases
c2c = c2c_pyroxene
cf = ca_ferrite_structured_phase
cpx = clinopyroxene
gt = garnet
il = akimotoite
ilmenite_group = akimotoite
mw = ferropericlase
magnesiowuestite = ferropericlase
ol = mg_fe_olivine
opx = orthopyroxene
plag = plagioclase
ppv = post_perovskite
pv = mg_fe_perovskite
mg_fe_bridgmanite = mg_fe_perovskite
mg_fe_silicate_perovskite = mg_fe_perovskite
ri = mg_fe_ringwoodite
spinel_group=mg_fe_aluminous_spinel
wa = mg_fe_wadsleyite
spinelloid_III = mg_fe_wadsleyite
