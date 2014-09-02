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
        base_material = [[diopside(), '[Ca][Mg][Si]2O6'],[hedenbergite(), '[Ca][Fe][Si]2O6'],[clinoenstatite(), '[Mg][Mg][Si]2O6'],[ca_tschermaks_molecule(), '[Ca][Al][Si1/2Al1/2]2O6'],[jadeite(), '[Na][Al][Si]2O6']]

        # Interaction parameters
        enthalpy_interaction=[[0., 24.74e3, 26.e3, 24.3e3],[24.74e3, 0., 0.e3], [60.53136e3, 0.0], [10.e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )

class garnet(SolidSolution):
    def __init__(self):
        # Name
        self.name='garnet'

        # Endmembers (garnet is symmetric)
        base_material = [[pyrope(), '[Mg]3[Al][Al]Si3O12'],[almandine(), '[Fe]3[Al][Al]Si3O12'],[grossular(), '[Ca]3[Al][Al]Si3O12'],[mg_majorite(), '[Mg]3[Mg][Si]Si3O12'],[jd_majorite(), '[Na2/3Al1/3]3[Al][Si]Si3O12']]
        # Interaction parameters
        enthalply_interaction=[[0.0, 30.e3, 21.20278e3, 0.0],[0.0,0.0,0.0],[57.77596e3, 0.0],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )


class akimotoite(SolidSolution):
    def __init__(self):
        # Name
        self.name='akimotoite/ilmenite'

        # Endmembers (ilmenite/akimotoite is symmetric)
        base_material = [[mg_akimotoite(), '[Mg][Si]O3'],[fe_akimotoite(), '[Fe][Si]O3'],[corundum(), '[Al][Al]O3']]
        # Interaction parameters
        enthalpy_interaction=[[0.0, 66.e3],[66.e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )

class ferropericlase(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesiowustite/ferropericlase'

        # Endmembers (ferropericlase is symmetric)
        base_material = [[periclase(), '[Mg]O'],[wuestite(), '[Fe]O']]
        # Interaction parameters
        enthalpy_interaction=[[13.e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )

class mg_fe_olivine(SolidSolution):
    def __init__(self):
        # Name
        self.name='olivine'

        # Endmembers (olivine is symmetric)
        base_material = [[forsterite(), '[Mg]2SiO4'],[fayalite(), '[Fe]2SiO4']]
        # Interaction parameters
        enthalpy_interaction=[[7.81322e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )

class orthopyroxene(SolidSolution):
    def __init__(self):
        # Name
        self.name='orthopyroxene'

        # Endmembers (orthopyroxene is symmetric)
        base_material = [[enstatite(), '[Mg][Mg][Si]SiO6'],[ferrosilite(), '[Fe][Fe][Si]SiO6'],[mg_tschermaks_molecule(), '[Mg][Al][Al]SiO6'],[orthodiopside(), '[Ca][Mg][Si]SiO6']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 0.0, 32.11352e3],[0.0, 0.0],[48.35316e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class plagioclase(SolidSolution):
    def __init__(self):
        # Name
        self.name='plagioclase'

        # Endmembers (plagioclase is symmetric)
        base_material = [[anorthite(), '[Ca][Al]2Si2O8'],[albite(), '[Na][Al1/2Si1/2]2Si2O8']]
        # Interaction parameters
        enthalpy_interaction=[[26.0e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class post_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='post-perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_post_perovskite(), '[Mg][Si]O3'],[fe_post_perovskite(), '[Fe][Si]O3'],[al_post_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 60.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class mg_fe_perovskite(SolidSolution):
    def __init__(self):
        # Name
        self.name='magnesium silicate perovskite/bridgmanite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_perovskite(), '[Mg][Si]O3'],[fe_perovskite(), '[Fe][Si]O3'],[al_perovskite(), '[Al][Al]O3']]

        # Interaction parameters
        enthalpy_interaction=[[0.0, 116.0e3],[0.0]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class mg_fe_ringwoodite(SolidSolution):
    def __init__(self):
        # Name
        self.name='ringwoodite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_ringwoodite(), '[Mg]2SiO4'],[fe_ringwoodite(), '[Fe]2SiO4']]

        # Interaction parameters
        enthalpy_interaction=[[9.34084e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class mg_fe_aluminous_spinel(SolidSolution):
    def __init__(self):
        # Name
        self.name='spinel-hercynite binary, fixed order'

        # Endmembers (post perovskite is symmetric)
        base_material = [[spinel(), '[Mg3/4Al1/4]4[Al7/8Mg1/8]8O16'],[hercynite(), '[Fe3/4Al1/4]4[Al7/8Fe1/8]8O16']]

        # Interaction parameters
        enthalpy_interaction=[[5.87646e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction))

class mg_fe_wadsleyite(SolidSolution):
    def __init__(self):
        # Name
        self.name='wadsleyite'

        # Endmembers (post perovskite is symmetric)
        base_material = [[mg_wadsleyite(), '[Mg]2SiO4'],[fe_wadsleyite(), '[Fe]2SiO4']]

        # Interaction parameters
        enthalpy_interaction=[[16.74718e3]]

        SolidSolution.__init__(self, base_material, SymmetricVanLaar(base_material, enthalpy_interaction) )

'''
ENDMEMBERS
'''
class anorthite (Mineral):
    def __init__(self):
       formula='CaAl2Si2O8'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'anorthite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -4015.e3,
            'V_0': 100.61e-6,
            'K_0': 84.e9,
            'Kprime_0': 4.0,
            'Debye_0': 752.,
            'grueneisen_0': 0.39,
            'q_0': 1.0,
            'G_0': 40e9,
            'Gprime_0': 1.1,
            'eta_s_0': 1.6}

       self.uncertainties = {
            'err_F_0': 4.e3,
            'err_V_0': 0.,
            'err_K_0': 5.e9,
            'err_Kprime_0': 1.0,
            'err_Debye_0': 2.,
            'err_grueneisen_0': 0.05,
            'err_q_0': 1.0,
            'err_G_0': 3.e9,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 0.1
            }

class albite (Mineral):
    def __init__(self):
       formula='NaAlSi3O8'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'albite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -3719.e3,
            'V_0': 100.45e-6,
            'K_0': 60.e9,
            'Kprime_0': 4.0,
            'Debye_0': 716.,
            'grueneisen_0': 0.57,
            'q_0': 1.0,
            'G_0': 36.e9,
            'Gprime_0': 1.4,
            'eta_s_0': 1.0}

       self.uncertainties = {
            'err_F_0': 5.e3,
            'err_V_0': 0.0,
            'err_K_0': 5.e9,
            'err_Kprime_0': 1.0,
            'err_Debye_0': 13,
            'err_grueneisen_0': 0.03,
            'err_q_0': 1.0,
            'err_G_0': 5.e9,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0
            }

class spinel (Mineral):
    def __init__(self):
       formula='Mg4Al8O16'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'spinel',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -8668.0e3,
            'V_0': 159.05e-6,
            'K_0': 197.e9,
            'Kprime_0': 5.7,
            'Debye_0': 843.,
            'grueneisen_0': 1.02,
            'q_0': 2.7,
            'G_0': 108.e9,
            'Gprime_0': 0.4,
            'eta_s_0': 2.7}

       self.uncertainties = {
            'err_F_0': 32.e3,
            'err_V_0': 0.0,
            'err_K_0': 1.e9,
            'err_Kprime_0': 0.2,
            'err_Debye_0': 33.,
            'err_grueneisen_0': 0.04,
            'err_q_0': 0.6,
            'err_G_0': 10.e9,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 0.6
            }

class hercynite (Mineral):
    def __init__(self):
       formula='Fe4Al8O16'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'hercynite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -7324.e3,
            'V_0': 163.37e-6,
            'K_0': 209.e9,
            'Kprime_0': 5.7,
            'Debye_0': 763.,
            'grueneisen_0': 1.22,
            'q_0': 2.7,
            'G_0': 84.e9,
            'Gprime_0': 0.4,
            'eta_s_0': 2.8}

       self.uncertainties = {
            'err_F_0': 35.e3,
            'err_V_0': 0.0,
            'err_K_0': 2.e9,
            'err_Kprime_0': 1.0,
            'err_Debye_0': 32.,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.0,
            'err_G_0': 13.e9,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 1.0
            }

class forsterite (Mineral):
    def __init__(self):
       formula='Mg2SiO4'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'forsterite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -2055.e3,
            'V_0': 43.60e-6,
            'K_0': 128.e9,
            'Kprime_0': 4.2,
            'Debye_0': 809.,
            'grueneisen_0': 0.99,
            'q_0': 2.1,
            'G_0': 82.e9,
            'Gprime_0': 1.5,
            'eta_s_0': 2.3}

       self.uncertainties = {
            'err_F_0': 2.e3,
            'err_V_0': 0.0,
            'err_K_0': 2.e9,
            'err_Kprime_0': 0.2,
            'err_Debye_0': 1.,
            'err_grueneisen_0': 0.03,
            'err_q_0': 0.2,
            'err_G_0': 2.e9,
            'err_Gprime_0': 0.1,
            'err_eta_s_0': 0.1
            }

class fayalite (Mineral):
    def __init__(self):
       formula='Fe2SiO4'
       formula=dictionarize_formula(formula)
       self.params = {
            'name': 'fayalite',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': -1371.e3,
            'V_0': 46.29e-6,
            'K_0': 135.e9,
            'Kprime_0': 4.2,
            'Debye_0': 619.,
            'grueneisen_0': 1.06,
            'q_0': 3.6,
            'G_0': 51.e9,
            'Gprime_0': 1.5,
            'eta_s_0': 1.0}

       self.uncertainties = {
            'err_F_0': 1.e3,
            'err_V_0': 0.0,
            'err_K_0': 2.e9,
            'err_Kprime_0': 1.0,
            'err_Debye_0': 2.,
            'err_grueneisen_0': 0.07,
            'err_q_0': 1.0,
            'err_G_0': 2.e9,
            'err_Gprime_0': 0.5,
            'err_eta_s_0': 0.6
            }

class mg_wadsleyite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_wadsleyite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_ringwoodite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_ringwoodite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class enstatite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class ferrosilite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_tschermaks_molecule (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class orthodiopside (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class diopside (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class hedenbergite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class clinoenstatite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class ca_tschermaks_molecule (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class jadeite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class hp_clinoenstatite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class hp_clinoferrosilite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class ca_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_akimotoite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_akimotoite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class corundum (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class pyrope (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class almandine (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class grossular (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_majorite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class jd_majorite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class quartz (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class coesite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class stishovite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class seifertite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class al_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_post_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_post_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class al_post_perovskite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class periclase (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class wuestite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class mg_calcium_ferrite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class fe_calcium_ferrite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class na_calcium_ferrite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class kyanite (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }

class nepheline (Mineral):
    def __init__(self):
       formula=''
       formula=dictionarize_formula(formula)
       self.params = {
            'name': '',
            'formula': formula,
            'equation_of_state': 'slb3',
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'F_0': ,
            'V_0': ,
            'K_0': ,
            'Kprime_0': ,
            'Debye_0': ,
            'grueneisen_0': ,
            'q_0': ,
            'G_0': ,
            'Gprime_0': ,
            'eta_s_0': }

        self.uncertainties = {
            'err_F_0': ,
            'err_V_0': ,
            'err_K_0': ,
            'err_Kprime_0': ,
            'err_Debye_0': ,
            'err_grueneisen_0': ,
            'err_q_0': ,
            'err_G_0': ,
            'err_Gprime_0': ,
            'err_eta_s_0': 
            }


# Mineral aliases
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
