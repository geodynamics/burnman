# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2019 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
from .. import constants


def calculate_nakajima_fp_pv_partition_coefficient(
    pressure, temperature, bulk_composition_mol, initial_distribution_coefficient
):
    """
    Calculate the partitioning of iron between periclase and bridgmanite as given
    by Nakajima et al., 2012.

    :param pressure: Equilibrium pressure [Pa]
    :type pressure: float
    :param temperature: Equilibrium temperature [K]
    :type temperature: float
    :param bulk_composition_mol: Bulk composition [mol].
        Only Mg, Fe, and Si are assumed to explicitly affect the partitioning
        with Al playing an implicit role.
    :type bulk_composition_mol: dict
    :param initial_distribution_coefficient: The distribution coefficient (Kd_0)
        at 25 GPa and 0 K.
    :type initial_distribution_coefficient: float

    :returns: The proportion of Fe in ferropericlase and perovskite, respectively
    :rtype: tuple
    """

    norm = bulk_composition_mol["Mg"] + bulk_composition_mol["Fe"]
    f_FeO = bulk_composition_mol["Fe"] / norm
    f_SiO2 = bulk_composition_mol["Si"] / norm

    Kd_0 = initial_distribution_coefficient
    delV = 2.0e-7  # in m^3/mol, average taken from Nakajima et al 2012, JGR

    # eq 5 Nakajima et al 2012, JGR. Solved for ln(K(P,T,X))
    rs = ((25.0e9 - pressure) * delV / (constants.gas_constant * temperature)) + np.log(
        Kd_0
    )

    # The exchange coefficent at P and T. K(P,T,X) in eq 5 Nakajima et al 2012
    K = np.exp(rs)

    # Solving equation 6 in Nakajima et al., 2012 for X_Fe_fp and X_Fe_pv
    # Solved using the definition of the distribution coefficient
    # to define X_Fe_fp as a function of X_Fe_pv

    num_to_sqrt = (-4.0 * f_FeO * (K - 1.0) * K * f_SiO2) + (
        pow(1.0 + (f_FeO * (K - 1)) + ((K - 1.0) * f_SiO2), 2.0)
    )

    X_Fe_pv = (
        -1.0 + f_FeO - (f_FeO * K) + f_SiO2 - (f_SiO2 * K) + np.sqrt(num_to_sqrt)
    ) / (2.0 * f_SiO2 * (1.0 - K))

    X_Fe_fp = X_Fe_pv / (((1.0 - X_Fe_pv) * K) + X_Fe_pv)

    return (X_Fe_fp, X_Fe_pv)
