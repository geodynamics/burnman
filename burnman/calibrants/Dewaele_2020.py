# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
import numpy as np


"""
Dewaele_2020
^^^^^^^^^^^^
"""


class CsCl(Calibrant):
    """
    The CsCl pressure standard reported by
    Dewaele (2020; https://doi.org/10.1080/08957959.2020.1774754).
    The thermal model is from Decker (1971).
    """

    def __init__(self):
        def _pressure_Dewaele_CsCl(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Total pressure
            P = P0 + Pth - Pth0

            return P

        _params_Dewaele_CsCl = {
            "V_0": 4.2179e-5,
            "K_0": 16.74e9,
            "Kprime_0": 5.703,
            "Debye_0": 151.0,
            "grueneisen_0": 1.99,
            "q_0": 1.18,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(
            self, _pressure_Dewaele_CsCl, "pressure", _params_Dewaele_CsCl
        )
