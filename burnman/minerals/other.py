# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


"""
Other minerals
^^^^^^^^^^^^^^

"""

import burnman.mineral_helpers as bmb
from burnman.mineral import Mineral

from SLB_2011 import periclase, wuestite, mg_perovskite, fe_perovskite

# This is ferropericlase with the deprecated solid solution setup, although it is still used in some of the /misc/paper* scripts
class ferropericlase(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        endmembers = [periclase(), wuestite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, endmembers, molar_fractions)

# this is mg_fe_perovskite with the depricated solid solution setup. Better not use...
class mg_fe_perovskite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        endmembers = [mg_perovskite(), fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, endmembers, molar_fractions)

# similar to ferropericlase, using the old solid solution setup. These values are based on Zhang, Stixrude and Brodholt 2013
class ZSB_2013_mg_fe_perovskite(bmb.HelperSolidSolution):
    def __init__(self, fe_num):
        endmembers = [ZSB_2013_mg_perovskite(), ZSB_2013_fe_perovskite()]
        molar_fractions = [1. - fe_num, 0.0 + fe_num] # keep the 0.0 +, otherwise it is an array sometimes
        bmb.HelperSolidSolution.__init__(self, endmembers, molar_fractions)

class ZSB_2013_mg_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
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
            'eta_s_0': 2.13 } #2.6

        self.uncertainties = {
            'err_K_0': 3.e9,
            'err_Kprime_0': 0.1,
            'err_G_0': 2.e9,
            'err_Gprime_0' : 0.0,
            'err_Debye_0': 5.,
            'err_grueneisen_0':.05,
            'err_q_0': .3,
            'err_eta_s_0':.3}

        Mineral.__init__(self)

class ZSB_2013_fe_perovskite(Mineral):
    def __init__(self):
        self.params = {
            'equation_of_state':'slb3',
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
            'eta_s_0': 2.3 } #2.3

        self.uncertainties = {
            'err_K_0':40e9,
            'err_Kprime_0':1.,
            'err_G_0':40e9,
            'err_Gprime_0':0.0,
            'err_Debye_0':26.,
            'err_grueneisen_0':.3,
            'err_q_0':1.0,
            'err_eta_s_0':1.0}

        Mineral.__init__(self)


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
        Mineral.__init__(self)

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
        Mineral.__init__(self)





