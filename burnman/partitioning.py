# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
import numpy as np
from . import constants

# TODO: add up weight percent and check <100 and tell them how much

molar_mass = {'Fe': 55.845 / 1000., 'Mg': 24.305 / 1000., 'O': 15.999 / 1000.,
              'Al': 26.982 / 1000., 'Ca': 40.078 / 1000., 'Si': 28.085 / 1000.}  # kg/mol
lower_mantle_mass = 4.043e24 * .75  # in kg


# convert weight percentage (amount, 1.00 = 100%) of a given element to
# molar mass
def weight_pct_to_mol(element, amount):

    return amount * lower_mantle_mass / molar_mass[element] * constants.Avogadro


def calculate_phase_percents(inp):
    """
    Converts given weight percentages into the requisite percent of each phase
    in mols and also returns the fraction of perovskite versus ferropericlase,
    assuming all of the silcon goes into the perovskite phase
    and with any remaining Fe or Mg going into the oxide phase.
    Input:
    inp={'Mg': ..., 'Fe': ..., ...} # in weight percent
    Returns:
    phase_per={'fp': ..., 'pv': ...} # as a fraction
    rel_mol_per={'MgO: ..., 'FeO': ..., ...} # in mols
    """
    names = {'Mg': 'MgO', 'Fe': 'FeO', 'Si': 'SiO2', 'Ca': 'Ca', 'Al': 'Al'}
    rel_mol_per = {}
    out = {}
    for a in inp:
        out[names[a]] = weight_pct_to_mol(a, inp[a])

    norm = out['MgO'] + out['FeO']
    for a in inp:
        rel_mol_per[names[a]] = out[names[a]] / norm

    frac_mol_SiO2 = rel_mol_per['SiO2']
    phase_per = {'fp': (1. - frac_mol_SiO2), 'pv': frac_mol_SiO2}
    return phase_per, rel_mol_per


def calculate_partition_coefficient(pressure, temperature, components, initial_distribution_coefficient):
    """ calculate the partition coefficient given [...] initial_distribution_coefficient is known as Kd_0 """

    frac_mol_FeO = components['FeO']
    frac_mol_SiO2 = components['SiO2']
    Kd_0 = initial_distribution_coefficient

    delV = 2.e-7  # in m^3/mol, average taken from Nakajima et al 2012, JGR

    # eq 5 Nakajima et al 2012, JGR. Solved for ln(K(P,T,X))
    rs = ((25.e9 - pressure) * (delV) / (constants.gas_constant * temperature)) + \
        np.log(Kd_0)

    # The exchange coefficent at P and T. K(P,T,X) in eq 5 Nakajima et al 2012
    K = np.exp(rs)

    # Solving equation 6 in Nakajima et al., 2012 for X_Fe_fp (we define as a) and X_Fe_pv (we define as b)
    # Solved using the definition of the distribution coefficient to define X_Fe_fp as a function of X_Fe_pv

    num_to_sqrt = (-4. * frac_mol_FeO * (K - 1.) * K * frac_mol_SiO2) + (
        pow(1. + (frac_mol_FeO * (K - 1)) + ((K - 1.) * frac_mol_SiO2), 2.))

    b = (-1. + frac_mol_FeO - (frac_mol_FeO * K) + frac_mol_SiO2 - (frac_mol_SiO2 * K) + np.sqrt(num_to_sqrt)) \
        / (2. * frac_mol_SiO2 * (1. - K))

    a = b / (((1. - b) * K) + b)

    return (a, b)  # a is partition coefficient array with P for fp, b is pcarray for pv
