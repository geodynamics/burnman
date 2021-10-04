# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
SLB_2011_ZSB_2013
^^^^^^^^^^^^^^^^^

Minerals from Stixrude & Lithgow-Bertelloni 2011, Zhang, Stixrude & Brodholt 2013, and references therein.

"""
from __future__ import absolute_import

from ..classes import mineral_helpers as helpers
from ..classes.mineral import Mineral


class stishovite (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 14.02e-6,
            'K_0': 314.0e9,
            'Kprime_0': 3.8,
            'G_0': 220.0e9,
            'Gprime_0': 1.9,
            'molar_mass': .0601,
            'n': 3,
            'Debye_0': 1108.,
            'grueneisen_0': 1.37,
            'q_0': 2.8,
            'eta_s_0': 4.6}

        self.uncertainties = {
            'err_K_0': 8.e9,
            'err_Kprime_0': 0.1,
            'err_G_0': 12.e9,
            'err_Gprime_0': 0.1,
            'err_Debye_0': 13.,
            'err_grueneisen_0': .17,
            'err_q_0': 2.2,
            'err_eta_s_0': 1.0
        }

        Mineral.__init__(self)


class periclase (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.8,
            'G_0': 131.0e9,
            'Gprime_0': 2.1,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 767.,
            'grueneisen_0': 1.36,
            'q_0': 1.7,  # 1.7
            'eta_s_0': 2.8}  # 2.8

        self.uncertainties = {
            'err_K_0': 3.e9,
            'err_Kprime_0': .2,
            'err_G_0': 1.0e9,
            'err_Gprime_0': .1,
            'err_Debye_0': 9.,
            'err_grueneisen_0': .05,
            'err_q_0': .2,
            'err_eta_s_0': .2}

        Mineral.__init__(self)


class wuestite (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 12.26e-6,
            'K_0': 179.0e9,
            'Kprime_0': 4.9,
            'G_0': 59.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 454.,
            'grueneisen_0': 1.53,
            'q_0': 1.7,  # 1.7
            'eta_s_0': -0.1}

        self.uncertainties = {
            'err_K_0': 1.e9,
            'err_Kprime_0': .2,
            'err_G_0': 1.e9,
            'err_Gprime_0': .1,
            'err_Debye_0': 21.,
            'err_grueneisen_0': .13,
            'err_q_0': 1.0,
            'err_eta_s_0': 1.0}

        Mineral.__init__(self)


class mg_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 24.45e-6,
            'K_0': 250.5e9,
            'Kprime_0': 4.01,
            'G_0': 172.9e9,
            'Gprime_0': 1.74,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 905.9,
            'grueneisen_0': 1.44,
            'q_0': 1.09,
            'eta_s_0': 2.13}  # 2.6

        self.uncertainties = {
            'err_K_0': 3.e9,
            'err_Kprime_0': 0.1,
            'err_G_0': 2.e9,
            'err_Gprime_0': 0.0,
            'err_Debye_0': 5.,
            'err_grueneisen_0': .05,
            'err_q_0': .3,
            'err_eta_s_0': .3}

        Mineral.__init__(self)


class fe_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'slb3',
            'V_0': 25.49e-6,
            'K_0': 272.0e9,
            'Kprime_0': 4.1,
            'G_0': 133.0e9,
            'Gprime_0': 1.4,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 871.,
            'grueneisen_0': 1.57,
            'q_0': 1.1,
            'eta_s_0': 2.3}  # 2.3

        self.uncertainties = {
            'err_K_0': 40e9,
            'err_Kprime_0': 1.,
            'err_G_0': 40e9,
            'err_Gprime_0': 0.0,
            'err_Debye_0': 26.,
            'err_grueneisen_0': .3,
            'err_q_0': 1.0,
            'err_eta_s_0': 1.0}

        Mineral.__init__(self)

mg_bridgmanite = mg_perovskite
fe_bridgmanite = fe_perovskite
