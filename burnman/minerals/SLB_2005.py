# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
SLB_2005
^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2005 and references therein

"""
from __future__ import absolute_import
from ..classes import mineral_helpers as helpers
from ..classes.mineral import Mineral


class stishovite (Mineral):

    def __init__(self):
        self.params = {
            'formula': {'Si': 1., 'O': 2.},
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'Kprime_0': 4.4,
            'G_0': 220.0e9,
            'Gprime_0': 1.6,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1044.,
            'grueneisen_0': 1.34,
            'q_0': 2.4,
            'eta_s_0': 5.0}

        Mineral.__init__(self)


class periclase (Mineral):

    def __init__(self):
        self.params = {
            'formula': {'Mg': 1., 'O': 1.},
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.8}
        Mineral.__init__(self)


class wuestite (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'formula': {'Fe': 1., 'O': 1.},
            'V_0': 12.06e-6,
            'K_0': 152.0e9,
            'Kprime_0': 4.9,
            'G_0': 47.0e9,
            'Gprime_0': 0.7,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 455.,
            'grueneisen_0': 1.28,
            'q_0': 1.5,
            'eta_s_0': 0.8}

        Mineral.__init__(self)


class mg_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'formula': {'Mg': 1., 'Si': 1., 'O': 3.},
            'equation_of_state': 'slb3',
            'V_0': 24.45e-6,
            'K_0': 251.0e9,
            'Kprime_0': 4.1,
            'G_0': 175.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 1070.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.6}

        Mineral.__init__(self)


class fe_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'formula': {'Fe': 1., 'Si': 1., 'O': 3.},
            'equation_of_state': 'slb3',
            'V_0': 25.48e-6,
            'K_0': 281.0e9,
            'Kprime_0': 4.1,
            'G_0': 138.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 841.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.1}

        Mineral.__init__(self)

mg_bridgmanite = mg_perovskite
fe_bridgmanite = fe_perovskite
