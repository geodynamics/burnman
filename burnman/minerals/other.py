# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Other minerals
^^^^^^^^^^^^^^

"""

import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral


class Speziale_fe_periclase(bmb.HelperSpinTransition):
    def __init__(self):
        bmb.HelperSpinTransition.__init__(self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
        self.cite = 'Speziale et al. 2007'


class Speziale_fe_periclase_HS(Mineral):
    """
    Speziale et al. 2007, Mg#=83
    """
    def __init__(self):
            self.params = {
                        'equation_of_state': 'mgd3',
                        'V_0': 22.9e-6,
                        'K_0': 157.5e9,
                        'Kprime_0': 3.92,
                        'molar_mass': .04567,
                        'n': 2,
                        'Debye_0': 587,
                        'grueneisen_0': 1.46,
                        'q_0': 1.2 }


class Speziale_fe_periclase_LS(Mineral):
    """
    Speziale et al. 2007, Mg#=83
    """
    def __init__(self):
        self.params = {
                        'equation_of_state': 'mgd3',
                        'V_0': 21.49e-6,
                        'K_0': 186.0e9,
                        'Kprime_0': 4.6,
                        'molar_mass': .04567,
                        'n': 2,
                        'Debye_0': 587.,
                        'grueneisen_0': 1.46,
                        'q_0': 1.2  }

"""
Katsura_2009
^^^^^^^^^^^^^

Minerals from Katsura 2009 and references therein

"""

class Katsura_2009_wadsleyite (Mineral):
    """
    Katsura 2009 and references therein
    """
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
            'V_0': 4.0499E-05,
            'K_0': 169.2e9,
            'Kprime_0': 4.1,
            'G_0': 113e9,
            'Gprime_0': 1.5,
            'molar_mass': .140695,
            'n': 7,
            'Debye_0': 814.,
            'grueneisen_0': 1.64,
            'q_0': 1.5,
            'eta_s_0': 2.3 }





