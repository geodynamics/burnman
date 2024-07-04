# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Zha_2004
^^^^^^^^
"""


class Re(Calibrant):
    """
    The Re pressure standard reported by
    Zha et al. (2004; https://doi.org/10.1063/1.1765752).
    """

    def __init__(self):
        def _pressure_Zha_Re(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth = (
                params["aK(V,T)"] + params["dK_dT_V"] * np.log(params["V_0"] / volume)
            ) * (temperature - params["T_0"])

            # Total pressure
            P = P0 + Pth

            return P

        _params_Zha_Re = {
            "V_0": 8.85516e-06,
            "K_0": 360.0e9,
            "Kprime_0": 4.5,
            "aK(V,T)": 0.00776e9,
            "dK_dT_V": -0.00815e9,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Zha_Re, "pressure", _params_Zha_Re)
