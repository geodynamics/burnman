# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np


def attenuation_correction(v_p, v_s, v_phi, Qs, Qphi):
    """
    Applies the attenuation correction following :cite:`Matas2007`, page 4.
    This is simplified, and there is also currently no 1D Q model implemented.
    The correction, however, only slightly reduces the velocities,
    and can be ignored for our current applications.
    Arguably, it might not be as relevant when comparing computations
    to PREM for periods of 1s as is implemented here.
    Called from :func:`burnman.tools.seismic.attenuation_correction`

    :param v_p: P wave velocity in [m/s].
    :type v_p: float
    :param v_s: S wave velocitiy in [m/s].
    :type v_s: float
    :param v_phi: Bulk sound velocity in [m/s].
    :type v_phi: float
    :param Qs: Shear quality factor [dimensionless].
    :type Qs: float
    :param Qphi: Bulk quality factor [dimensionless].
    :type Qphi: float

    :returns: Corrected P wave, S wave and bulk sound velocities in [m/s].
    :rtype: tuple
    """
    beta = 0.3  # Matas et al. (2007) page 4
    Qp = 3.0 / 4.0 * pow((v_p / v_s), 2.0) * Qs  # Matas et al. (2007) page 4

    cot = 1.0 / np.tan(beta * np.pi / 2.0)
    v_p *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qp  # Matas et al. (2007) page 1
    v_s *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qs
    v_phi *= 1.0 - 1.0 / 2.0 * cot * 1.0 / Qphi
    return v_p, v_s, v_phi
