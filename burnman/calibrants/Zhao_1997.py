# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Zhao_1997
^^^^^^^^^
"""


class hBN(Calibrant):
    """
    The hBN pressure standard reported by
    Zhao et al. (1997; https://doi.org/10.1080/08957959708240481).
    """

    def __init__(self):
        def _pressure_Zhao_hBN(volume, temperature, params):

            # Modified BM3+Thermal (see Zhao et al. 1997)
            a0 = params["a0"]
            a1 = params["a1"]
            K0 = params["K_0"]
            Kprime0 = params["Kprime_0"]
            dK_dT = params["dK_dT"]
            dKprime_dT = params["dKprime_dT"]

            KT = (K0) + dK_dT * (temperature - 300.0)
            KprimeT = Kprime0 + dKprime_dT * (temperature - 300.0)
            V0T = params["V_0"] * np.exp(
                ((a0 * temperature) + (a1 * ((temperature**2.0) / 2.0)))
                - ((a0 * 300.0) + (a1 * ((300.0**2.0) / 2.0)))
            )
            V_V0T = volume / V0T
            f = 0.5 * ((V_V0T ** (-2.0 / 3.0)) - 1.0)

            P = (3.0 * KT * f * (1.0 + 2.0 * f) ** (5.0 / 2.0)) * (
                1.0 - 3.0 / 2.0 * (4.0 - KprimeT) * f
            )

            return P

        _params_Zhao_hBN = {
            "V_0": 1.08164e-05,
            "K_0": 17.6e9,
            "Kprime_0": 19.5,
            "a0": 4.38e-05,
            "a1": 1.75e-08,
            "dK_dT": -0.0069e9,
            "dKprime_dT": 0.0,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Zhao_hBN, "pressure", _params_Zhao_hBN)
