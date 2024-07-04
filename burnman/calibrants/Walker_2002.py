# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Walker_2002
^^^^^^^^^^^
"""


class KCl_B2(Calibrant):
    """
    The B2 KCl pressure standard reported by
    Walker (2002; https://doi.org/10.2138/am-2002-0701).
    """

    def __init__(self):
        def _pressure_Walker_KCl(volume, temperature, params):

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

        _params_Walker_KCl = {
            "V_0": 3.22365e-05,
            "K_0": 23.7e9,
            "Kprime_0": 4.4,
            "aK(V,T)": 0.00284e9,
            "dK_dT_V": 0.00012e9,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(self, _pressure_Walker_KCl, "pressure", _params_Walker_KCl)
