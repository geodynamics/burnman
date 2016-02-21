# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Sundman 1991 models for iron
Written in the Holland and Powell, 2011 form

EoS terms for bcc are from HP_2011_ds62 for iron
EoS terms for fcc are from an unpublished calibration
(R. Myhill, 02/2015)
"""
from __future__ import absolute_import

from ..mineral import Mineral
from ..solidsolution import SolidSolution
from ..solutionmodel import *
from ..processchemistry import read_masses, dictionarize_formula, formula_mass

atomic_masses = read_masses()


"""
ENDMEMBERS
"""


class bcc_iron (Mineral):

    def __init__(self):
        formula = 'Fe1.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'BCC iron',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': 9149.0,
            'S_0': 36.868,
            'V_0': 7.09e-06,
            'Cp': [21.09, 0.0101455, -221508., 47.1947],
            'a_0': 3.56e-05,
            'K_0': 1.64e+11,
            'Kprime_0': 5.16,
            'Kdprime_0': -3.1e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'curie_temperature': [1043., 0.0],
            'magnetic_moment': [2.22, 0.0],
            'magnetic_structural_parameter': 0.4}
        Mineral.__init__(self)


class fcc_iron (Mineral):

    def __init__(self):
        formula = 'Fe1.0'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'FCC iron',
            'formula': formula,
            'equation_of_state': 'hp_tmt',
            'H_0': 7973.0,
            'S_0': 35.907,
            'V_0': 6.93863394593e-06,
            'Cp': [22.24, 0.0088656, -221517., 47.1998],
            'a_0': 5.13e-05,
            'K_0': 1.539e+11,
            'Kprime_0': 5.2,
            'Kdprime_0': -3.37e-11,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula, atomic_masses),
            'curie_temperature': [201., 0.0],
            'magnetic_moment': [2.1, 0.0],
            'magnetic_structural_parameter': 0.28}
        Mineral.__init__(self)
