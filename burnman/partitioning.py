# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
import numpy as np
from . import constants

def calculate_nakajima_fp_pv_partition_coefficient(pressure, temperature, bulk_composition_mol, initial_distribution_coefficient):
    """ 
    calculate the partition coefficient given [...] initial_distribution_coefficient is known as Kd_0 
    """

    norm = bulk_composition_mol['Mg'] + bulk_composition_mol['Fe']

    f_FeO = bulk_composition_mol['Fe']/norm
    f_SiO2 = bulk_composition_mol['Si']/norm

    Kd_0 = initial_distribution_coefficient

    delV = 2.e-7  # in m^3/mol, average taken from Nakajima et al 2012, JGR

    # eq 5 Nakajima et al 2012, JGR. Solved for ln(K(P,T,X))
    rs = ((25.e9 - pressure) * (delV) / (constants.gas_constant * temperature)) + \
        np.log(Kd_0)

    # The exchange coefficent at P and T. K(P,T,X) in eq 5 Nakajima et al 2012
    K = np.exp(rs)

    # Solving equation 6 in Nakajima et al., 2012 for X_Fe_fp (we define as a) and X_Fe_pv (we define as b)
    # Solved using the definition of the distribution coefficient to define X_Fe_fp as a function of X_Fe_pv

    num_to_sqrt = ((-4. * f_FeO * (K - 1.) * K * f_SiO2) + 
                   (pow(1. + (f_FeO * (K - 1)) + ((K - 1.) * f_SiO2), 2.)))

    b = ((-1. + f_FeO - (f_FeO * K) + f_SiO2 - (f_SiO2 * K) + np.sqrt(num_to_sqrt)) /
         (2. * f_SiO2 * (1. - K)))

    a = b / (((1. - b) * K) + b)

    return (a, b)  # a is partition coefficient array with P for fp, b is pcarray for pv
