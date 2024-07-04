# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
import numpy as np


"""
Dewaele_2012
^^^^^^^^^^^^
"""


class KBr_B2(Calibrant):
    """
    The B2 KBr pressure standard reported by
    Dewaele et al. (2012; https://doi.org/10.1103/PhysRevB.85.214105).
    """

    def __init__(self):
        def _pressure_Dewaele_KBr(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth = (
                params["aK(V,T)"] + params["dK_dT_V"] * np.log(params["V_0"] / volume)
            ) * (temperature - params["T_0"])

            # Total pressure
            P = P0 + Pth

            return P

        _params_Dewaele_KBr = {
            "V_0": 3.8180e-5,
            "K_0": 14.9e9,
            "Kprime_0": 5.81,
            "aK(V,T)": 0.00222e9,
            "dK_dT_V": 0.0,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(self, _pressure_Dewaele_KBr, "pressure", _params_Dewaele_KBr)


class KCl_B2(Calibrant):
    """
    The B2 KCl pressure standard reported by
    Dewaele et al. (2012; https://doi.org/10.1103/PhysRevB.85.214105).
    """

    def __init__(self):
        def _pressure_Dewaele_KCl(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth = (
                params["aK(V,T)"] + params["dK_dT_V"] * np.log(params["V_0"] / volume)
            ) * (temperature - params["T_0"])

            # Total pressure
            P = P0 + Pth

            return P

        _params_Dewaele_KCl = {
            "V_0": 3.28206e-5,
            "K_0": 17.2e9,
            "Kprime_0": 5.89,
            "aK(V,T)": 0.00224e9,
            "dK_dT_V": 0.0,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(self, _pressure_Dewaele_KCl, "pressure", _params_Dewaele_KCl)
