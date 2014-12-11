from burnman.mineral import Mineral
from burnman.solidsolution import SolidSolution
from burnman.processchemistry import read_masses, dictionarize_formula, formula_mass
from burnman.solutionmodel import *

atomic_masses=read_masses()

class Mo (Mineral):
    def __init__(self):
       formula='Ni1.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Mo',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 28.59 ,
            'V_0': 9.391e-06 ,
            'Cp': [33.9, 0.006276, 38859.7, -12.0] ,
            'a_0': 1.44e-05 ,
            'K_0': 2.608e+11 ,
            'Kprime_0': 4.46 ,
            'Kdprime_0': -1.71e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class MoO2 (Mineral):
    def __init__(self):
       formula='Ni1.0O2.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'MoO2',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -591500.0 ,
            'S_0': 50.016 ,
            'V_0': 1.9799e-05 ,
            'Cp': [56.1, 0.02559, -17.6, 18.9] ,
            'a_0': 4.4e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.25e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class Re (Mineral):
    def __init__(self):
       formula='Ni1.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Re',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 36.53 ,
            'V_0': 8.862e-06 ,
            'Cp': [23.7, 0.005448, 68.0, 0.0] ,
            'a_0': 1.9e-05 ,
            'K_0': 3.6e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -1.1e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class ReO2 (Mineral):
    def __init__(self):
       formula='Ni1.0O2.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'ReO2',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -445140.0 ,
            'S_0': 47.82 ,
            'V_0': 1.8779e-05 ,
            'Cp': [76.89, 0.00993, -1207130.0, -208.0] ,
            'a_0': 4.4e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.25e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class ReO3 (Mineral):
    def __init__(self):
       formula='Ni1.0O3.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'ReO3',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -589107.0 ,
            'S_0': 69.26 ,
            'V_0': 3.3844e-05 ,
            'Cp': [152.1, -0.00322, -191287.0, -1288.0] ,
            'a_0': 4.4e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.25e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class Re2O7 (Mineral):
    def __init__(self):
       formula='Ni2.0O7.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Re2O7',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': -1263150.0 ,
            'S_0': 207.28 ,
            'V_0': 7.9372e-05 ,
            'Cp': [303.0, 0.0, 284633.0, -2430.0] ,
            'a_0': 4.4e-05 ,
            'K_0': 1.8e+11 ,
            'Kprime_0': 4.05 ,
            'Kdprime_0': -2.25e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class Pt (Mineral):
    def __init__(self):
       formula='Ni1.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'Pt',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 0.0 ,
            'S_0': 41.53 ,
            'V_0': 9.107e-06 ,
            'Cp': [23519.2, 6.28644, 0.0, 0.0] ,
            'a_0': 2.64e-05 ,
            'K_0': 2.77715e+11 ,
            'Kprime_0': 4.82 ,
            'Kdprime_0': -1.38e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}

class PtO2 (Mineral):
    def __init__(self):
       formula='Ni1.0O2.0'
       formula = dictionarize_formula(formula)
       self.params = {
            'name': 'PtO2',
            'formula': formula,
            'equation_of_state': 'mtait',
            'H_0': 171500.0 ,
            'S_0': 0.0 ,
            'V_0': 1.9244e-05 ,
            'Cp': [0.0, 0.0, 0.0, 0.0] ,
            'a_0': 1.44e-05 ,
            'K_0': 2.46e+11 ,
            'Kprime_0': 4.0 ,
            'Kdprime_0': -1.63e-11 ,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses)}


