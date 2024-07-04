# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Ono_2022
^^^^^^^^
"""


class Re(Calibrant):
    """
    The Re pressure standard reported by
    Ono (2022; https://doi.org/10.1155/2022/7545777).
    """

    def __init__(self):
        def _pressure_Ono_Re(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth = (
                params["dK_dT"] + params["dK_dT_V"] * np.log(params["V_0"] / volume)
            ) * (temperature - params["T_0"])

            # Total pressure
            P = P0 + Pth

            return P

        _params_Ono_Re = {
            "V_0": 8.87798e-06,
            "K_0": 384.0e9,
            "Kprime_0": 3.26,
            "dK_dT": 0.0056e9,  # Pa/K
            "dK_dT_V": -0.00042e9,  # Pa/K
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Ono_Re, "pressure", _params_Ono_Re)
