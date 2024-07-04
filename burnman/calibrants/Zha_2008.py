# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
import numpy as np


"""
Zha_2008
^^^^^^^^
"""


class Pt(Calibrant):
    """
    The Pt pressure standard reported by
    Zha (2008; https://doi.org/10.1063/1.2844358).
    """

    def __init__(self):
        def _pressure_Zha_Pt(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Electronic pressure
            Pel = (
                1.1916e-15 * temperature**4.0
                - 1.4551e-11 * temperature**3.0
                + 1.6209e-07 * temperature**2.0
                + 1.8269e-4 * temperature
                - 0.069
            ) * 1.0e09

            # Total pressure
            P = P0 + Pth - Pth0 + Pel

            return P

        _params_Zha_Pt = {
            "V_0": 9.0904e-06,
            "K_0": 273.5e9,
            "Kprime_0": 4.7,
            "Debye_0": 230.0,  # 370-405
            "grueneisen_0": 2.75,
            "q_0": 0.25,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Zha_Pt, "pressure", _params_Zha_Pt)
