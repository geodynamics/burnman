# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


"""
Murakami_2013
^^^^^^^^^^^^^

Minerals from Murakami 2013 and references therein.

"""
from __future__ import absolute_import

from .. import mineral_helpers as helpers
from ..mineral import Mineral


class periclase (Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 11.24e-6,
            'K_0': 161.0e9,
            'Kprime_0': 3.9,
            'G_0': 130.9e9,
            'Gprime_0': 1.92,
            'molar_mass': .0403,
            'n': 2,
            'Debye_0': 773.,
            'grueneisen_0': 1.5,
            'q_0': 1.5,
            'eta_s_0': 2.3 }
        Mineral.__init__(self)


class wuestite (Mineral):
    """
    Murakami 2013 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
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
            'eta_s_0': 0.8 }

        Mineral.__init__(self)

class mg_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 24.45e-6,
            'K_0': 253.0e9,
            'Kprime_0': 4.1,
            'G_0': 172.9e9,
            'Gprime_0': 1.56,
            'molar_mass': .1000,
            'n': 5,
            'Debye_0': 1100.,
            'grueneisen_0': 1.4,
            'q_0': 1.4,
            'eta_s_0': 2.6 }

        Mineral.__init__(self)

class fe_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb2',
            'V_0': 25.49e-6,
            'K_0': 281.0e9,
            'Kprime_0': 4.1,
            'G_0': 138.0e9,
            'Gprime_0': 1.7,
            'molar_mass': .1319,
            'n': 5,
            'Debye_0': 841.,
            'grueneisen_0': 1.48,
            'q_0': 1.4,
            'eta_s_0': 2.1 }

        Mineral.__init__(self)

mg_bridgmanite = mg_perovskite
fe_bridgmanite = fe_perovskite
