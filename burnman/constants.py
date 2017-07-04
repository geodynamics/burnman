# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
import scipy.constants

"""
molar gas constant (R) in J mol^-1 K^-1
"""
gas_constant = scipy.constants.gas_constant


"""
Avogadro constant (N_A) in mol ^ -1
"""
Avogadro = scipy.constants.Avogadro


"""
Boltzmann constant (k_B) in J K^-1.

Note that we are not using scipy.constants.Boltzmann because it is not
available in older versions.
"""
Boltzmann = 1.3806488e-23


"""
Newtonian constant of gravitation (G) in m^3 kg^-1 s^-2
"""

G = scipy.constants.G


"""
Dirac constant (hbar, Planck constant / 2*pi) in J s^-1
"""
Dirac = 1.054571726e-34


"""
1 cm^-1 in J/mol
"""
invcm = 11.9627
