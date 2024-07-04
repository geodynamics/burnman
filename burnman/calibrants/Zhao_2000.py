# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Zhao_2000
^^^^^^^^^
"""


class Mo(Calibrant):
    """
    The Mo pressure standard reported by
    Zhao et al. (2000; https://doi.org/10.1103/PhysRevB.62.8766).
    """

    def __init__(self):
        def _pressure_Zhao_Mo(volume, temperature, params):

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

        _params_Zhao_Mo = {
            "V_0": 9.37647e-06,
            "K_0": 268.0e9,
            "Kprime_0": 3.82,
            "a0": 1.27e-5,
            "a1": 1.12e-8,
            "dK_dT": -0.0213e9,
            "dKprime_dT": -1.41e-2,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Zhao_Mo, "pressure", _params_Zhao_Mo)
