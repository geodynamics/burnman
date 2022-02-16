# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
DKS_2013_solids
^^^^^^^^^^^^^^^

Solids from de Koker and Stixrude (2013) FPMD simulations
"""

from ..classes.mineral import Mineral
from ..utils.chemistry import formula_mass


class stishovite(Mineral):
    def __init__(self):
        p_fit = [0.2471304763E+05, 0.4793020138E+05]
        V_0 = 0.1513000000E+02 * 1e-6
        K_0 = p_fit[0] / (9. * V_0) * 1.e3
        a3 = 2. * p_fit[1] / (9. * K_0 * V_0) * 1.e3
        Kp_0 = 4. + ( a3 / 3. )
        Kdp_0 = ( -143./9. - Kp_0*(Kp_0 - 7.)) / K_0
        formula = {'Si': 1., 'O': 2.}
        self.params = {
            'name': 'stishovite',
            'formula': formula,
            'equation_of_state': 'dks_s',
            'V_0': V_0, # [m^3/mol]
            'T_0': 0.3000000000E+04, # [K]
            'E_0': -.2274840214E+04 * 1e3, # [J/mol]
            'S_0': 0.1668222552E+00 * 1e3, # [J/K/mol]
            'K_0': K_0,
            'Kprime_0': Kp_0,
            'Kdprime_0': Kdp_0,
            'n': 2., # called fsn in param file
            'Cv': 0.7794230433E-01 * 1e3, # [J/K/mol]
            'grueneisen_0': 0.1389501259E+01,
            'q_0': 0.1332025550E+01,
            'molar_mass': formula_mass(formula)
            }
        Mineral.__init__(self)


class perovskite(Mineral):
    def __init__(self):
        p_fit = [0.4067243956E+05, 0.1177159096E+05]
        V_0 = 0.2705000000E+02 * 1e-6
        K_0 = p_fit[0] / (9. * V_0) * 1.e3
        a3 = 2. * p_fit[1] / (9. * K_0 * V_0) * 1.e3
        Kp_0 = 4. + ( a3 / 3. )
        Kdp_0 = ( -143./9. - Kp_0*(Kp_0 - 7.)) / K_0
        formula = {'Mg': 1., 'Si': 1., 'O': 3.}
        self.params = {
            'name': 'perovskite',
            'formula': formula,
            'equation_of_state': 'dks_s',
            'V_0': V_0, # [m^3/mol]
            'T_0': 0.3000000000E+04, # [K]
            'E_0': -.3355012782E+04 * 1e3, # [J/mol]
            'S_0': 0.3384574347E+00 * 1e3, # [J/K/mol]
            'K_0': K_0,
            'Kprime_0': Kp_0,
            'Kdprime_0': Kdp_0,
            'n': 2., # called fsn in param file
            'Cv': 0.1338111589E+00 * 1e3, # [J/K/mol]
            'grueneisen_0': 0.1893754815E+01,
            'q_0': 0.1487809730E+01,
            'molar_mass': formula_mass(formula)
            }
        Mineral.__init__(self)


class periclase(Mineral):
    def __init__(self):
        p_fit = [0.1208938157E+05, 0.1133765229E+05]
        V_0 = 0.1223000000E+02 * 1e-6
        K_0 = p_fit[0] / (9. * V_0) * 1.e3
        a3 = 2. * p_fit[1] / (9. * K_0 * V_0) * 1.e3
        Kp_0 = 4. + ( a3 / 3. )
        Kdp_0 = ( -143./9. - Kp_0*(Kp_0 - 7.)) / K_0
        formula = {'Mg': 1., 'O': 1.}
        self.params = {
            'name': 'periclase',
            'formula': formula,
            'equation_of_state': 'dks_s',
            'V_0': V_0, # [m^3/mol]
            'T_0': 0.2000000000E+04, # [K]
            'E_0': -.1164949141E+04 * 1e3, # [J/mol]
            'S_0': 0.1198358648E+00 * 1e3, # [J/K/mol]
            'K_0': K_0,
            'Kprime_0': Kp_0,
            'Kdprime_0': Kdp_0,
            'n': 2., # called fsn in param file
            'Cv': 0.4904715075E-01 * 1e3, # [J/K/mol]
            'grueneisen_0': 0.1412003694E+01,
            'q_0': 0.6317609916E+00,
            'molar_mass': formula_mass(formula)
            }
        Mineral.__init__(self)
