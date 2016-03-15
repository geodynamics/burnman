# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Matas_etal_2007
^^^^^^^^^^^^^^^

Minerals from Matas et al. 2007 and references therein. See Table 1 and 2.
"""
from __future__ import absolute_import

from ..mineral import Mineral


class mg_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 24.43e-6,
            'K_0': 250.0e9,
            'Kprime_0': 4.0,
            'G_0': 175.0e9,
            'Gprime_0': 1.8,
            'molar_mass': .1020,
            'n': 5,
            'Debye_0': 1070.,
            'grueneisen_0': 1.48,
            'q_0': 1.4}

        Mineral.__init__(self)


class fe_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 25.34e-6,
            'K_0': 250.0e9,
            'Kprime_0': 4.0,
            'G_0': 135.0e9,
            'Gprime_0': 1.3,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 841.,
            'grueneisen_0': 1.48,
            'q_0': 1.4}

        Mineral.__init__(self)


class al_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 24.58e-6,
            'K_0': 249.0e9,
            'Kprime_0': 4.0,
            'G_0': 165.0e9,
            'Gprime_0': 1.8,
            'molar_mass': .1005,
            'n': 5,
            'Debye_0': 1021.,
            'grueneisen_0': 1.48,
            'q_0': 1.4}

        Mineral.__init__(self)


class ca_perovskite(Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 27.45e-6,
            'K_0': 236.0e9,
            'Kprime_0': 3.9,
            'G_0': 165.0e9,
            'Gprime_0': 2.46,
            'molar_mass': .11616,
            'n': 5,
            'Debye_0': 984.,
            'grueneisen_0': 1.53,
            'q_0': 1.6}

        Mineral.__init__(self)


class periclase (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 11.25e-6,
            'K_0': 160.1e9,
            'Kprime_0': 3.83,
            'G_0': 130.0e9,
            'Gprime_0': 2.2,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 673.,
            'grueneisen_0': 1.41,
            'q_0': 1.3}

        Mineral.__init__(self)


class wuestite (Mineral):

    def __init__(self):
        self.params = {
            'equation_of_state': 'mgd2',
            'V_0': 12.26e-6,
            'K_0': 160.1e9,
            'Kprime_0': 3.83,
            'G_0': 46.0e9,
            'Gprime_0':  0.6,
            'molar_mass': .0718,
            'n': 2,
            'Debye_0': 673.,
            'grueneisen_0': 1.41,
            'q_0': 1.3}

        Mineral.__init__(self)

ca_bridgmanite = ca_perovskite
mg_bridgmanite = mg_perovskite
fe_bridgmanite = fe_perovskite
al_bridgmanite = al_perovskite
