# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
DKS_2013_liquids
^^^^^^^^^^^^^^^^

Liquids from de Koker and Stixrude (2013) FPMD simulations.
"""

import numpy as np
from ..classes.mineral import Mineral
from ..utils.chemistry import formula_mass


# Vector parsing for DKS liquid equation of state
def vector_to_array(a, Of, Otheta):
    array = np.empty([Of+1, Otheta+1])
    for i in range(Of+1):
        for j in range(Otheta+1):
            n = int((i+j)*((i+j)+1.)/2. + j)
            array[i][j] = a[n]
    return array


class SiO2_liquid(Mineral):
    def __init__(self):
        formula = {'Si': 1.0, 'O': 2.}
        self.params = {
            'name': 'SiO2_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 2.78e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 5,
            'm': 0.91,
            'a': [-1945.93156, -226.6835978, 455.0286309, 2015.65287, -200.585046, -216.6028187, 48369.72992, 441.5340414, 73.07765325, 0.0, -651587.652, 20701.69954, 892.12209, 0.0, 0.0, 4100181.286, -128258.7237, -1228.478753, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.0004266056389,
            'xi': 0.8639433047,
            'Tel_0': 5651.204964,
            'eta': -0.2783503528,
            'el_V_0': 1e-06,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class MgSiO3_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 1.0, 'Si': 1.0, 'O': 3.}
        self.params = {
            'name': 'MgSiO3_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 4.18e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.83,
            'a': [-2984.241297, -380.9839126, 601.8088234, 7307.69753, 7.626381912, -328.367174, 38737.46417, 6251.230413, 402.4716495, 0.0, 0.0, -23578.93569, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.008009960983,
            'xi': -0.08859010337,
            'Tel_0': 2194.563521,
            'eta': -0.775354875,
            'el_V_0': 3.89008e-05,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class MgSi2O5_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 1.0, 'Si': 2.0, 'O': 5.}
        self.params = {
            'name': 'MgSi2O5_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 6.75e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.79,
            'a': [-4958.560203, -607.6635229, 1089.553108, 9125.144702, -443.9654989, -603.1466364, 62485.19233, 10927.5085, 1425.929331, 0.0, 0.0, -27738.0811, -4055.024972, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.02219035084,
            'xi': 0.7754599642,
            'Tel_0': 1699.783718,
            'eta': -0.4712864331,
            'el_V_0': 0.0001080643234,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class MgSi3O7_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 1.0, 'Si': 3.0, 'O': 7.}
        self.params = {
            'name': 'MgSi3O7_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 9.35e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.86,
            'a': [-6925.370617, -832.0455172, 1439.840307, 12287.35224, -264.7754561, -780.6835127, 78558.47091, 8145.563063, 1444.65423, 0.0, 0.0, -13036.84144, -2360.783631, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.02273836197,
            'xi': 0.4506392324,
            'Tel_0': 2085.530204,
            'eta': -0.4804823168,
            'el_V_0': 0.0001080643234,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class MgSi5O11_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 1.0, 'Si': 5.0, 'O': 11.}
        self.params = {
            'name': 'MgSi5O11_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 0.000146,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 4,
            'm': 0.77,
            'a': [-10813.78126, -1297.292175, 2642.979479, 19993.66381, -1085.821183, -1226.314792, 49132.58238, 18886.40475, 1252.739819, 0.0, 528358.8509, -54556.11837, 2844.969895, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.03318133543,
            'xi': 0.4033708612,
            'Tel_0': 2037.798559,
            'eta': -0.6209203711,
            'el_V_0': 0.00015,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class Mg2SiO4_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 2.0, 'Si': 1.0, 'O': 4.}
        self.params = {
            'name': 'Mg2SiO4_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 5.84e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.75,
            'a': [-3944.769208, -531.7975964, 880.0460994, 11401.47398, 118.7409191, -456.3140461, 55778.07008, 12132.5261, 519.3612273, 0.0, 0.0, -48733.22459, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.01101820277,
            'xi': 1.175924196,
            'Tel_0': 2228.185561,
            'eta': -0.464192202,
            'el_V_0': 5.23613e-05,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class Mg3Si2O7_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 3.0, 'Si': 2.0, 'O': 7.}
        self.params = {
            'name': 'Mg3Si2O7_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 0.0001005,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.79,
            'a': [-6945.262972, -905.8656523, 1466.115121, 18498.28462, 260.3083362, -841.8330982, 89795.95729, 14752.47411, 1120.541194, 0.0, 0.0, -55594.62308, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.02603891379,
            'xi': 1.129966677,
            'Tel_0': 2230.685379,
            'eta': -0.3689626876,
            'el_V_0': 0.0001080643234,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class Mg5SiO7_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 5.0, 'Si': 1.0, 'O': 7.}
        self.params = {
            'name': 'Mg5SiO7_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 0.0001075,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.64,
            'a': [-6721.181931, -1008.922671, 1800.764267, 25856.0057, 2169.612789, -753.9019178, 103374.7345, 17933.40061, 127.5989699, 0.0, 0.0, -80394.40732, 570.0622605, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.0163355197,
            'xi': 0.2784006205,
            'Tel_0': 1662.606581,
            'eta': -0.9693629899,
            'el_V_0': 0.0001080643234,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)


class MgO_liquid(Mineral):
    def __init__(self):
        formula = {'Mg': 1.0, 'O': 1.}
        self.params = {
            'name': 'MgO_liquid',
            'formula': formula,
            'equation_of_state': 'dks_l',
            'V_0': 1.646e-05,
            'T_0': 3000.0,
            'O_theta': 2,
            'O_f': 3,
            'm': 0.63,
            'a': [-925.2677296, -155.3240992, 260.8211743, 5323.167667, 466.3722398, -88.30035696, 10473.87879, 1997.967054, 50.72520834, 0.0, 0.0, -9914.621337, 71.89989255, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            'zeta_0': 0.002194565772,
            'xi': 0.411459446,
            'Tel_0': 1620.106387,
            'eta': -0.986457555,
            'el_V_0': 1.620953559e-05,
            'n': sum(formula.values()),
            'molar_mass': formula_mass(formula)
            }
        self.params['a'] = vector_to_array(self.params['a'], self.params['O_f'], self.params['O_theta'])*1e3 # [J/mol]
        Mineral.__init__(self)
