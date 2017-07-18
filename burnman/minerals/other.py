# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Other minerals
^^^^^^^^^^^^^^

"""
from __future__ import absolute_import

from .. import mineral_helpers as helpers
from ..mineral import Mineral
from ..processchemistry import dictionarize_formula, formula_mass

from .SLB_2011 import periclase, wuestite, mg_perovskite, fe_perovskite



class liquid_iron( Mineral ):
    """
    Liquid iron equation of state from Anderson and Ahrens (1994)
    """
    def __init__(self):
        formula='Fe1.0'
        formula = dictionarize_formula(formula)
        m = formula_mass(formula)
        rho_0 = 7019.
        V_0 = m/rho_0
        D = 7766.
        Lambda = 1146.
        self.params = {
            'name': 'liquid iron',
            'formula': formula,
            'equation_of_state': 'aa',
            'P_0': 1.e5,
            'T_0': 1811.,
            'S_0': 100., # to fit
            'molar_mass': m,
            'V_0': V_0,
            'E_0': 0.,
            'K_S': 109.7e9,
            'Kprime_S': 4.661,
            'Kprime_prime_S': -0.043e-9,
            'grueneisen_0': 1.735,
            'grueneisen_prime': -0.130/m*1.e-6,
            'grueneisen_n': -1.870,
            'a': [248.92*m, 289.48*m],
            'b': [0.4057*m, -1.1499*m],
            'Theta': [1747.3, 1.537],
            'theta': 5000.,
            'lmda': [-325.23*m, 302.07*m, 30.45*m],
            'xi_0': 282.67*m,
            'F': [D/rho_0, Lambda/rho_0],
            'n': sum(formula.values()),
            'molar_mass': m}
        Mineral.__init__(self)

        
class ZSB_2013_mg_perovskite(Mineral):

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


class ZSB_2013_fe_perovskite(Mineral):

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


class Speziale_fe_periclase(helpers.HelperSpinTransition):

    def __init__(self):
        helpers.HelperSpinTransition.__init__(
            self, 60.0e9, Speziale_fe_periclase_LS(), Speziale_fe_periclase_HS())
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
            'q_0': 1.2}
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
            'q_0': 1.2}
        Mineral.__init__(self)


class Liquid_Fe_Anderson(Mineral):

    """
    Anderson & Ahrens, 1994 JGR
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'bm4',
            'V_0': 7.95626e-6,
            'K_0': 109.7e9,
            'Kprime_0': 4.66,
            'Kprime_prime_0': -0.043e-9,
            'molar_mass': 0.055845,
        }
        Mineral.__init__(self)


class Fe_Dewaele(Mineral):

    """
    Dewaele et al., 2006, Physical Review Letters
    """

    def __init__(self):
        self.params = {
            'equation_of_state': 'vinet',
            'V_0': 6.75e-6,
            'K_0': 163.4e9,
            'Kprime_0': 5.38,
            'molar_mass': 0.055845,
            'n': 1}
        Mineral.__init__(self)
