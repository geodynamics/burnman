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
        p_fit = [0.2471304763e05, 0.4793020138e05]
        V_0 = 0.1513000000e02 * 1e-6
        K_0 = p_fit[0] / (9.0 * V_0) * 1.0e3
        a3 = 2.0 * p_fit[1] / (9.0 * K_0 * V_0) * 1.0e3
        Kp_0 = 4.0 + (a3 / 3.0)
        Kdp_0 = (-143.0 / 9.0 - Kp_0 * (Kp_0 - 7.0)) / K_0
        formula = {"Si": 1.0, "O": 2.0}
        self.params = {
            "name": "stishovite",
            "formula": formula,
            "equation_of_state": "dks_s",
            "V_0": V_0,  # [m^3/mol]
            "T_0": 0.3000000000e04,  # [K]
            "E_0": -0.2274840214e04 * 1e3,  # [J/mol]
            "S_0": 0.1668222552e00 * 1e3,  # [J/K/mol]
            "K_0": K_0,
            "Kprime_0": Kp_0,
            "Kdprime_0": Kdp_0,
            "n": 2.0,  # called fsn in param file
            "Cv": 0.7794230433e-01 * 1e3,  # [J/K/mol]
            "grueneisen_0": 0.1389501259e01,
            "q_0": 0.1332025550e01,
            "molar_mass": formula_mass(formula),
        }
        Mineral.__init__(self)


class perovskite(Mineral):
    def __init__(self):
        p_fit = [0.4067243956e05, 0.1177159096e05]
        V_0 = 0.2705000000e02 * 1e-6
        K_0 = p_fit[0] / (9.0 * V_0) * 1.0e3
        a3 = 2.0 * p_fit[1] / (9.0 * K_0 * V_0) * 1.0e3
        Kp_0 = 4.0 + (a3 / 3.0)
        Kdp_0 = (-143.0 / 9.0 - Kp_0 * (Kp_0 - 7.0)) / K_0
        formula = {"Mg": 1.0, "Si": 1.0, "O": 3.0}
        self.params = {
            "name": "perovskite",
            "formula": formula,
            "equation_of_state": "dks_s",
            "V_0": V_0,  # [m^3/mol]
            "T_0": 0.3000000000e04,  # [K]
            "E_0": -0.3355012782e04 * 1e3,  # [J/mol]
            "S_0": 0.3384574347e00 * 1e3,  # [J/K/mol]
            "K_0": K_0,
            "Kprime_0": Kp_0,
            "Kdprime_0": Kdp_0,
            "n": 2.0,  # called fsn in param file
            "Cv": 0.1338111589e00 * 1e3,  # [J/K/mol]
            "grueneisen_0": 0.1893754815e01,
            "q_0": 0.1487809730e01,
            "molar_mass": formula_mass(formula),
        }
        Mineral.__init__(self)


class periclase(Mineral):
    def __init__(self):
        p_fit = [0.1208938157e05, 0.1133765229e05]
        V_0 = 0.1223000000e02 * 1e-6
        K_0 = p_fit[0] / (9.0 * V_0) * 1.0e3
        a3 = 2.0 * p_fit[1] / (9.0 * K_0 * V_0) * 1.0e3
        Kp_0 = 4.0 + (a3 / 3.0)
        Kdp_0 = (-143.0 / 9.0 - Kp_0 * (Kp_0 - 7.0)) / K_0
        formula = {"Mg": 1.0, "O": 1.0}
        self.params = {
            "name": "periclase",
            "formula": formula,
            "equation_of_state": "dks_s",
            "V_0": V_0,  # [m^3/mol]
            "T_0": 0.2000000000e04,  # [K]
            "E_0": -0.1164949141e04 * 1e3,  # [J/mol]
            "S_0": 0.1198358648e00 * 1e3,  # [J/K/mol]
            "K_0": K_0,
            "Kprime_0": Kp_0,
            "Kdprime_0": Kdp_0,
            "n": 2.0,  # called fsn in param file
            "Cv": 0.4904715075e-01 * 1e3,  # [J/K/mol]
            "grueneisen_0": 0.1412003694e01,
            "q_0": 0.6317609916e00,
            "molar_mass": formula_mass(formula),
        }
        Mineral.__init__(self)
