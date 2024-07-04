# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.classes.calibrant import Calibrant
import numpy as np


"""
Dubrovinsky_1998
^^^^^^^^^^^^^^^^
"""


class Al2O3_corundum(Calibrant):
    """
    The Al2O3 pressure standard reported by
    Dubrovinsky (1998; https://dx.doi.org/10.1007/s002690050133).
    """

    def __init__(self):
        def _pressure_Dubrovinsky_Al2O3(volume, temperature, params):

            # Modified BM3+Thermal
            a0 = params["a"]
            a1 = params["b"]
            a2 = params["c"]
            K0 = params["K_0"]
            Kprime0 = params["Kprime_0"]
            dK_dT = params["dK_dT_P"]
            dKprime_dT = params["dKprime_dT"]

            KT = (K0) + dK_dT * (temperature - params["T_0"])
            KprimeT = Kprime0 + dKprime_dT * (temperature - params["T_0"])
            V0T = params["V_0"] * np.exp(
                (
                    (a0 * temperature)
                    + (a1 * ((temperature**2.0) / 2.0))
                    - (a2 / temperature)
                )
                - (
                    (a0 * params["T_0"])
                    + (a1 * ((params["T_0"] ** 2.0) / 2.0))
                    - (a2 / params["T_0"])
                )
            )
            V_V0T = volume / V0T
            f = 0.5 * ((V_V0T ** (-2.0 / 3.0)) - 1.0)

            P = (3.0 * KT * f * (1.0 + 2.0 * f) ** (5.0 / 2.0)) * (
                1.0 - 3.0 / 2.0 * (4.0 - KprimeT) * f
            )

            return P

        _params_Dubrovinsky_Al2O3 = {
            "V_0": 2.558e-05,
            "K_0": 253.0e9,
            "Kprime_0": 5.0,
            "a": 2.515e-5,
            "b": 1.958e-9,
            "c": -0.305,
            "dK_dT_P": -0.02e9,
            "dKprime_dT": 0.0,
            "n": 5.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 6.0,
        }

        Calibrant.__init__(
            self, _pressure_Dubrovinsky_Al2O3, "pressure", _params_Dubrovinsky_Al2O3
        )
