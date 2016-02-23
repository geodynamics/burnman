# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

HP_2011_fluids
^^^^^^^^

Fluids from Holland and Powell 2011 and references therein.
CORK parameters:
CHO gases from Holland and Powell, 1991. ["CO2",304.2,0.0738],["CH4",190.6,0.0460],["H2",41.2,0.0211],["CO",132.9,0.0350]
H2O and S2 from Wikipedia, 2012/10/23. ["H2O",647.096,0.22060],["S2",1314.00,0.21000]
H2S from ancyclopedia.airliquide.com, 2012/10/23. ["H2S",373.15,0.08937]

NB: Units for cork[i] in Holland and Powell datasets are
a = kJ^2/kbar*K^(1/2)/mol^2 -> multiply by 1e-2
b = kJ/kbar/mol -> multiply by 1e-5
c = kJ/kbar^1.5/mol -> multiply by 1e-9
d = kJ/kbar^2/mol -> multiply by 1e-13

Individual terms are divided through by P, P, P^1.5, P^2, so
[0][j] -> multiply by 1e6
[1][j] -> multiply by 1e3
[2][j] -> multiply by 1e3
[3][j] -> multiply by 1e3

cork_P: kbar -> multiply by 1e8
"""
from __future__ import absolute_import

from ..mineral import Mineral
from ..processchemistry import read_masses, dictionarize_formula, formula_mass

atomic_masses = read_masses()


class CO2 (Mineral):

    def __init__(self):
        formula = 'CO2'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'carbon dioxide',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 304.2,
            'cork_P': 0.0738e8,
            'H_0': -393.51e3,
            'S_0': 213.7,
            'Cp': [87.8, -2.644e-3, 706.4e3, -998.9]}
        Mineral.__init__(self)


class CH4 (Mineral):

    def __init__(self):
        formula = 'CH4'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'methane',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 190.6,
            'cork_P': 0.0460e8,
            'H_0': -74.81e3,
            'S_0': 186.26,
            'Cp': [150.1, 0.002063, 3427700., -2650.4]}
        Mineral.__init__(self)


class O2 (Mineral):

    def __init__(self):
        formula = 'O2'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'oxygen gas',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 0.,
            'cork_P': 1.0e5,
            'H_0': 0.,
            'S_0': 205.2,
            'Cp': [48.3, -0.000691, 499200., -420.7]}
        Mineral.__init__(self)


class H2 (Mineral):

    def __init__(self):
        formula = 'H2'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'hydrogen gas',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 41.2,
            'cork_P': 0.0211e8,
            'H_0': 0.,
            'S_0': 130.7,
            'Cp': [23.3, 0.004627, 0.0, 76.3]}
        Mineral.__init__(self)


class S2 (Mineral):

    def __init__(self):
        formula = 'S2'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'sulfur gas',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 1314.00,
            'cork_P': 0.21000e8,
            'H_0': 128.54e3,
            'S_0': 231.0,
            'Cp': [37.1, 0.002398, -161000.0, -65.0]}
        Mineral.__init__(self)


class H2S (Mineral):

    def __init__(self):
        formula = 'H2S'
        formula = dictionarize_formula(formula)
        self.params = {
            'name': 'hydrogen sulfide',
            'formula': formula,
            'equation_of_state': 'cork',
            'cork_params': [[5.45963e1, -8.63920e0], [9.18301e-1], [-3.30558e-2, 2.30524e-3], [6.93054e-4, -8.38293e-5]],
            'cork_T': 373.15,
            'cork_P': 0.08937e8,
            'H_0': 128.54e3,
            'S_0': 231.0,
            'Cp': [47.4, 0.010240, 615900., -397.8]}
        Mineral.__init__(self)
