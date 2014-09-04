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

atomic_masses=read_masses('data/input_masses/atomic_masses.dat')

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
        base_material = [[mg_calcium_ferrite(), '[Mg]Al[Al]O4'],[fe_calcium_ferrite(), '[Fe]Al[Al]O4'],[na_calcium_ferrite(), '[Na]Al[Si]O4']]

        SolidSolution.__init__(self, base_material, IdealSolution(base_material))

class ca_rich_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='calcium perovskite'

        # Endmembers (cpv is symmetric)
        base_material = [[ca_perovskite(), '[Ca][Ca]Si2O6'],[jd_perovskite(), '[Na][Al]Si2O6']]

        SolidSolution.__init__(self, base_material, IdealSolution(base_materia))

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
        base_material = [[enstatite(), '[Mg][Mg][Si]SiO6'],[ferrosilite(), '[Fe][Fe][Si]SiO6'],[mg_tschermaks(), '[Mg][Al][Al]SiO6'],[orthodiopside(), '[Ca][Mg][Si]SiO6']]

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
            'F_0': -4015.0 ,
            'V_0': 100.61 ,
            'K_0': 84.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 752.0 ,
            'grueneisen_0': 0.39 ,
            'q_0': 1.0 ,
            'G_0': 40.0 ,
            'Gprime_0': 1.1 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 3.0 ,
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
            'F_0': -3719.0 ,
            'V_0': 100.45 ,
            'K_0': 60.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 716.0 ,
            'grueneisen_0': 0.57 ,
            'q_0': 1.0 ,
            'G_0': 36.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 5.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5.0 ,
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
            'F_0': -8668.0 ,
            'V_0': 159.05 ,
            'K_0': 197.0 ,
            'Kprime_0': 5.7 ,
            'Debye_0': 843.0 ,
            'grueneisen_0': 1.02 ,
            'q_0': 2.7 ,
            'G_0': 108.0 ,
            'Gprime_0': 0.4 ,
            'eta_s_0': 2.7 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 32.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 33.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.6 ,
            'err_G_0': 10.0 ,
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
            'F_0': -7324.0 ,
            'V_0': 163.37 ,
            'K_0': 209.0 ,
            'Kprime_0': 5.7 ,
            'Debye_0': 763.0 ,
            'grueneisen_0': 1.22 ,
            'q_0': 2.7 ,
            'G_0': 84.0 ,
            'Gprime_0': 0.4 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 35.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 32.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 13.0 ,
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
            'F_0': -2055.0 ,
            'V_0': 43.6 ,
            'K_0': 128.0 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 809.0 ,
            'grueneisen_0': 0.99 ,
            'q_0': 2.1 ,
            'G_0': 82.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 1.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 0.2 ,
            'err_G_0': 2.0 ,
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
            'F_0': -1371.0 ,
            'V_0': 46.29 ,
            'K_0': 135.0 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 619.0 ,
            'grueneisen_0': 1.06 ,
            'q_0': 3.6 ,
            'G_0': 51.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2.0 ,
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
            'F_0': -2028.0 ,
            'V_0': 40.52 ,
            'K_0': 169.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 844.0 ,
            'grueneisen_0': 1.21 ,
            'q_0': 2.0 ,
            'G_0': 112.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.09 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2.0 ,
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
            'F_0': -1365.0 ,
            'V_0': 42.8 ,
            'K_0': 169.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 665.0 ,
            'grueneisen_0': 1.21 ,
            'q_0': 2.0 ,
            'G_0': 72.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 7.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 13.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 21.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 12.0 ,
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
            'F_0': -2017.0 ,
            'V_0': 39.49 ,
            'K_0': 185.0 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 878.0 ,
            'grueneisen_0': 1.11 ,
            'q_0': 2.4 ,
            'G_0': 123.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.1 ,
            'err_q_0': 0.4 ,
            'err_G_0': 2.0 ,
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
            'F_0': -1363.0 ,
            'V_0': 41.86 ,
            'K_0': 213.0 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 679.0 ,
            'grueneisen_0': 1.27 ,
            'q_0': 2.4 ,
            'G_0': 92.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 7.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.23 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -2913.0 ,
            'V_0': 62.68 ,
            'K_0': 107.0 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 812.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 77.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.4 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.4 ,
            'err_G_0': 1.0 ,
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
            'F_0': -2226.0 ,
            'V_0': 65.94 ,
            'K_0': 101.0 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 674.0 ,
            'grueneisen_0': 0.72 ,
            'q_0': 3.4 ,
            'G_0': 52.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4.0 ,
            'err_K_prime_0': 0.5 ,
            'err_Debye_0': 10.0 ,
            'err_grueneisen_0': 0.08 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5.0 ,
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
            'F_0': -3003.0 ,
            'V_0': 59.14 ,
            'K_0': 107.0 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 784.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 97.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 9.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 24.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -3016.0 ,
            'V_0': 68.05 ,
            'K_0': 107.0 ,
            'Kprime_0': 7.0 ,
            'Debye_0': 745.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 3.4 ,
            'G_0': 60.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 1.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -3030.0 ,
            'V_0': 66.04 ,
            'K_0': 112.0 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 782.0 ,
            'grueneisen_0': 0.96 ,
            'q_0': 1.5 ,
            'G_0': 67.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5.0 ,
            'err_K_prime_0': 1.8 ,
            'err_Debye_0': 3.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 2.0 ,
            'err_G_0': 2.0 ,
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
            'F_0': -2677.0 ,
            'V_0': 67.87 ,
            'K_0': 119.0 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 702.0 ,
            'grueneisen_0': 0.94 ,
            'q_0': 1.5 ,
            'G_0': 61.0 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 45.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
            'F_0': -2906.0 ,
            'V_0': 62.5 ,
            'K_0': 112.0 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 805.0 ,
            'grueneisen_0': 0.96 ,
            'q_0': 1.5 ,
            'G_0': 81.0 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 1.7 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 10.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class ca-tschermaks (Mineral):
    def __init__(self):
       formula='CaAl2SiO6'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Ca-Tschermaks',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -3120.0 ,
            'V_0': 63.57 ,
            'K_0': 112.0 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 804.0 ,
            'grueneisen_0': 0.78 ,
            'q_0': 1.5 ,
            'G_0': 76.0 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 5.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.0 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -2855.0 ,
            'V_0': 60.51 ,
            'K_0': 142.0 ,
            'Kprime_0': 5.2 ,
            'Debye_0': 821.0 ,
            'grueneisen_0': 0.9 ,
            'q_0': 0.4 ,
            'G_0': 85.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.2 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 12.0 ,
            'err_grueneisen_0': 0.08 ,
            'err_q_0': 1.4 ,
            'err_G_0': 2.0 ,
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
            'F_0': -2905.0 ,
            'V_0': 60.76 ,
            'K_0': 116.0 ,
            'Kprime_0': 6.2 ,
            'Debye_0': 824.0 ,
            'grueneisen_0': 1.12 ,
            'q_0': 0.2 ,
            'G_0': 88.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.5 ,
            'err_G_0': 1.0 ,
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
            'F_0': -2222.0 ,
            'V_0': 63.85 ,
            'K_0': 116.0 ,
            'Kprime_0': 6.2 ,
            'Debye_0': 692.0 ,
            'grueneisen_0': 1.12 ,
            'q_0': 0.2 ,
            'G_0': 71.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 0.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 11.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -1463.0 ,
            'V_0': 27.45 ,
            'K_0': 236.0 ,
            'Kprime_0': 3.9 ,
            'Debye_0': 796.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 0.9 ,
            'G_0': 157.0 ,
            'Gprime_0': 2.2 ,
            'eta_s_0': 1.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 8.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 44.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.6 ,
            'err_G_0': 12.0 ,
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
            'F_0': -1410.0 ,
            'V_0': 26.35 ,
            'K_0': 211.0 ,
            'Kprime_0': 5.6 ,
            'Debye_0': 934.0 ,
            'grueneisen_0': 1.19 ,
            'q_0': 2.3 ,
            'G_0': 132.0 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 4.0 ,
            'err_K_prime_0': 0.8 ,
            'err_Debye_0': 12.0 ,
            'err_grueneisen_0': 0.13 ,
            'err_q_0': 0.8 ,
            'err_G_0': 8.0 ,
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
            'F_0': -1068.0 ,
            'V_0': 26.85 ,
            'K_0': 211.0 ,
            'Kprime_0': 5.6 ,
            'Debye_0': 888.0 ,
            'grueneisen_0': 1.19 ,
            'q_0': 2.3 ,
            'G_0': 150.0 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 3.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 21.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 120.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -1582.0 ,
            'V_0': 25.58 ,
            'K_0': 253.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 933.0 ,
            'grueneisen_0': 1.32 ,
            'q_0': 1.3 ,
            'G_0': 163.0 ,
            'Gprime_0': 1.6 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 5.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 3.0 ,
            'err_grueneisen_0': 0.04 ,
            'err_q_0': 0.2 ,
            'err_G_0': 2.0 ,
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
            'F_0': -5936.0 ,
            'V_0': 113.08 ,
            'K_0': 170.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 823.0 ,
            'grueneisen_0': 1.01 ,
            'q_0': 1.4 ,
            'G_0': 94.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 10.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 0.5 ,
            'err_G_0': 2.0 ,
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
            'F_0': -4935.0 ,
            'V_0': 115.43 ,
            'K_0': 174.0 ,
            'Kprime_0': 4.9 ,
            'Debye_0': 741.0 ,
            'grueneisen_0': 1.06 ,
            'q_0': 1.4 ,
            'G_0': 96.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 29.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
            'F_0': -6278.0 ,
            'V_0': 125.12 ,
            'K_0': 167.0 ,
            'Kprime_0': 3.9 ,
            'Debye_0': 823.0 ,
            'grueneisen_0': 1.05 ,
            'q_0': 1.9 ,
            'G_0': 109.0 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 11.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 2.0 ,
            'err_grueneisen_0': 0.06 ,
            'err_q_0': 0.2 ,
            'err_G_0': 4.0 ,
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
            'F_0': -5691.0 ,
            'V_0': 114.32 ,
            'K_0': 165.0 ,
            'Kprime_0': 4.2 ,
            'Debye_0': 822.0 ,
            'grueneisen_0': 0.98 ,
            'q_0': 1.5 ,
            'G_0': 85.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 10.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3.0 ,
            'err_K_prime_0': 0.3 ,
            'err_Debye_0': 4.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 0.5 ,
            'err_G_0': 2.0 ,
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
            'F_0': -5519.0 ,
            'V_0': 110.94 ,
            'K_0': 177.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 896.0 ,
            'grueneisen_0': 1.01 ,
            'q_0': 1.4 ,
            'G_0': 125.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 3.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 14.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 7.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 18.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 4.0 ,
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
            'F_0': -859.0 ,
            'V_0': 23.67 ,
            'K_0': 50.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 816.0 ,
            'grueneisen_0': 0.0 ,
            'q_0': 1.0 ,
            'G_0': 45.0 ,
            'Gprime_0': 1.0 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 31.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
            'F_0': -855.0 ,
            'V_0': 20.66 ,
            'K_0': 114.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 857.0 ,
            'grueneisen_0': 0.39 ,
            'q_0': 1.0 ,
            'G_0': 62.0 ,
            'Gprime_0': 1.2 ,
            'eta_s_0': 2.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
            'F_0': -819.0 ,
            'V_0': 14.02 ,
            'K_0': 314.0 ,
            'Kprime_0': 3.8 ,
            'Debye_0': 1108.0 ,
            'grueneisen_0': 1.37 ,
            'q_0': 2.8 ,
            'G_0': 220.0 ,
            'Gprime_0': 1.9 ,
            'eta_s_0': 4.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 8.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.17 ,
            'err_q_0': 2.2 ,
            'err_G_0': 12.0 ,
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
            'F_0': -794.0 ,
            'V_0': 13.67 ,
            'K_0': 328.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 1141.0 ,
            'grueneisen_0': 1.37 ,
            'q_0': 2.8 ,
            'G_0': 227.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 5.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 2.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 16.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 2.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class mg_perovskite (Mineral):
    def __init__(self):
       formula='MgSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg-Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1368.0 ,
            'V_0': 24.45 ,
            'K_0': 251.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 905.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 173.0 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 2.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 5.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.3 ,
            'err_G_0': 2.0 ,
            'err_Gprime_0': 0.0 ,
            'err_eta_s_0': 0.3 }

class fe_perovskite (Mineral):
    def __init__(self):
       formula='FeSiO3'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe-Perovskite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1041.0 ,
            'V_0': 25.49 ,
            'K_0': 272.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 871.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 133.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 2.3 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 6.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 40.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 26.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 40.0 ,
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
            'F_0': -1534.0 ,
            'V_0': 24.94 ,
            'K_0': 258.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 886.0 ,
            'grueneisen_0': 1.57 ,
            'q_0': 1.1 ,
            'G_0': 171.0 ,
            'Gprime_0': 1.5 ,
            'eta_s_0': 2.5 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 2.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 0.5 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -1348.0 ,
            'V_0': 24.42 ,
            'K_0': 231.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 855.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 1.1 ,
            'G_0': 150.0 ,
            'Gprime_0': 2.0 ,
            'eta_s_0': 1.2 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 7.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 0.1 ,
            'err_G_0': 4.0 ,
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
            'F_0': -982.0 ,
            'V_0': 25.46 ,
            'K_0': 231.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 782.0 ,
            'grueneisen_0': 1.89 ,
            'q_0': 1.1 ,
            'G_0': 129.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': 1.4 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 21.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 52.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 5.0 ,
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
            'F_0': -1378.0 ,
            'V_0': 23.85 ,
            'K_0': 249.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 762.0 ,
            'grueneisen_0': 1.65 ,
            'q_0': 1.1 ,
            'G_0': 92.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 20.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.02 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -569.0 ,
            'V_0': 11.24 ,
            'K_0': 161.0 ,
            'Kprime_0': 3.8 ,
            'Debye_0': 767.0 ,
            'grueneisen_0': 1.36 ,
            'q_0': 1.7 ,
            'G_0': 131.0 ,
            'Gprime_0': 2.1 ,
            'eta_s_0': 2.8 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 0.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 3.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 9.0 ,
            'err_grueneisen_0': 0.05 ,
            'err_q_0': 0.2 ,
            'err_G_0': 1.0 ,
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
            'F_0': -242.0 ,
            'V_0': 12.26 ,
            'K_0': 179.0 ,
            'Kprime_0': 4.9 ,
            'Debye_0': 454.0 ,
            'grueneisen_0': 1.53 ,
            'q_0': 1.7 ,
            'G_0': 59.0 ,
            'Gprime_0': 1.4 ,
            'eta_s_0': -0.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 1.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.2 ,
            'err_Debye_0': 21.0 ,
            'err_grueneisen_0': 0.13 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class mg_calcium_ferrite (Mineral):
    def __init__(self):
       formula='MgAl2O4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mg_Calcium_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -2122.0 ,
            'V_0': 36.18 ,
            'K_0': 211.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 838.0 ,
            'grueneisen_0': 1.31 ,
            'q_0': 1.0 ,
            'G_0': 130.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 2.1 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 16.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
            'err_Gprime_0': 0.1 ,
            'err_eta_s_0': 1.0 }

class fe_calcium_ferrite (Mineral):
    def __init__(self):
       formula='FeAl2O4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Fe_Calcium_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1790.0 ,
            'V_0': 37.26 ,
            'K_0': 211.0 ,
            'Kprime_0': 4.1 ,
            'Debye_0': 804.0 ,
            'grueneisen_0': 1.31 ,
            'q_0': 1.0 ,
            'G_0': 152.0 ,
            'Gprime_0': 1.8 ,
            'eta_s_0': 3.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 25.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 10.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 69.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
            'err_Gprime_0': 0.5 ,
            'err_eta_s_0': 1.0 }

class na_calcium_ferrite (Mineral):
    def __init__(self):
       formula='NaAlSiO4'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Na_Calcium_Ferrite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'F_0': -1851.0 ,
            'V_0': 36.27 ,
            'K_0': 158.0 ,
            'Kprime_0': 4.3 ,
            'Debye_0': 812.0 ,
            'grueneisen_0': 1.17 ,
            'q_0': 1.0 ,
            'G_0': 121.0 ,
            'Gprime_0': 2.1 ,
            'eta_s_0': 1.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 11.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.1 ,
            'err_Debye_0': 51.0 ,
            'err_grueneisen_0': 0.3 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
            'F_0': -2446.0 ,
            'V_0': 44.23 ,
            'K_0': 160.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 943.0 ,
            'grueneisen_0': 0.93 ,
            'q_0': 1.0 ,
            'G_0': 121.0 ,
            'Gprime_0': 1.7 ,
            'eta_s_0': 3.0 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 4.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 0.0 ,
            'err_Debye_0': 8.0 ,
            'err_grueneisen_0': 0.07 ,
            'err_q_0': 1.0 ,
            'err_G_0': 10.0 ,
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
            'F_0': -1993.0 ,
            'V_0': 54.67 ,
            'K_0': 53.0 ,
            'Kprime_0': 4.0 ,
            'Debye_0': 701.0 ,
            'grueneisen_0': 0.69 ,
            'q_0': 1.0 ,
            'G_0': 31.0 ,
            'Gprime_0': 1.3 ,
            'eta_s_0': 0.6 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

       self.uncertainties = {
            'err_F_0': 3.0 ,
            'err_V_0': 0.0 ,
            'err_K_0': 1.0 ,
            'err_K_prime_0': 1.0 ,
            'err_Debye_0': 13.0 ,
            'err_grueneisen_0': 0.03 ,
            'err_q_0': 1.0 ,
            'err_G_0': 1.0 ,
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
mgts = mg_tschermaks_molecule
odi = orthodiopside

# Clinopyroxenes
di = diopside
he = hedenbergite
cen = clinoenstatite
cats = ca_tschermaks_molecule
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
capv = calcium_perovskite
jdpv = jd_perovskite

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
mgcf = mg_calcium_ferrite
fecf = fe_calcium_ferrite
nacf = na_calcium_ferrite

# Al2SiO5 polymorphs
ky = kyanite

# Nepheline group
neph = nepheline



# Solid solution aliases
c2c = c2c_pyroxene
cf = calcium_ferrite_structured_phase
cpv = ca_rich_perovskite
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
