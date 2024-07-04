# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.classes.calibrant import Calibrant
import numpy as np


"""
LeGodec_2000
^^^^^^^^^^^^
"""


class hBN(Calibrant):
    """
    The hexagonal boron nitride pressure standard reported by
    LeGodec et al. (2000; http://dx.doi.org/10.1080/08957950008200304).
    """

    def __init__(self):
        def _pressure_LeGodec_hBN(volume, temperature, params):

            # Modified BM3+Thermal (see Zhao et al. 1997)
            a0 = params["a0"]
            a1 = params["a1"]
            K0 = params["K_0"]
            Kprime0 = params["Kprime_0"]
            dK_dT = params["dK_dT"]
            dKprime_dT = params["dKprime_dT"]

            KT = (K0) + dK_dT * (temperature - params["T_0"])
            KprimeT = Kprime0 + dKprime_dT * (temperature - params["T_0"])
            V0T = params["V_0"] * np.exp(
                ((a0 * temperature) + (a1 * ((temperature**2.0) / 2.0)))
                - ((a0 * params["T_0"]) + (a1 * ((params["T_0"] ** 2.0) / 2.0)))
            )
            V_V0T = volume / V0T
            f = 0.5 * ((V_V0T ** (-2.0 / 3.0)) - 1.0)

            P = (3.0 * KT * f * (1.0 + 2.0 * f) ** (5.0 / 2.0)) * (
                1.0 - 3.0 / 2.0 * (4.0 - KprimeT) * f
            )

            return P

        _params_LeGodec_hBN = {
            "V_0": 1.0891e-05,
            "K_0": 27.6e9,
            "Kprime_0": 10.5,
            "a0": 3.53e-05,
            "a1": 2.1e-09,
            "dK_dT": -0.0081e9,
            "dKprime_dT": 0.0016,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_LeGodec_hBN, "pressure", _params_LeGodec_hBN)
