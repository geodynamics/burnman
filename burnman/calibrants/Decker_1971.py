# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
Decker_1971
^^^^^^^^^^^
"""

import numpy as np
from ..eos.birch_murnaghan import birch_murnaghan
from ..eos import MGD2
from ..eos import debye
from ..classes.calibrant import Calibrant


class NaCl_B1(Calibrant):
    """
    The NaCl (B1 structured) pressure standard reported by Decker (1971).

    Note: This calibrant is not exactly the same as that proposed by Decker.
    The cold compression curve has here been approximated by a
    Birch-Murnaghan EoS.

    TODO: Make the calibrant exactly match that published by Decker.
    """
    def __init__(self):

        def _pressure_Decker_NaCl(volume, temperature, params):
            p300 = birch_murnaghan(params['V_0'] / volume, params)
            grueneisen = MGD2._grueneisen_parameter(0.,
                                                    params['V_0'] / volume,
                                                    params)
            Debye_T = params['Debye_0'] * np.exp((params['grueneisen_0']
                                                  - grueneisen)
                                                 / params['q_0'])
            Eqh = debye.thermal_energy(temperature, Debye_T, params['n'])
            EqhR = debye.thermal_energy(params['refT'], Debye_T, params['n'])
            dpqh = (Eqh - EqhR) * (grueneisen / volume)
            pressure = p300 + dpqh
            return pressure

        _params_Decker_NaCl = {'V_0': 2.7015e-05,
                               'K_0': 24.02e9,
                               'Kprime_0': 4.7369,
                               'Debye_0': 279,
                               'grueneisen_0': 1.59,
                               'q_0': 0.93,
                               'n': 2.,
                               'refT': 298.15,
                               'P_0': 0,
                               'Z': 4.}

        Calibrant.__init__(self, _pressure_Decker_NaCl, 'pressure',
                           _params_Decker_NaCl)
