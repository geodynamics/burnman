# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Pigott_2015
^^^^^^^^^^^
"""


class Ni(Calibrant):
    """
    The Ni pressure standard reported by
    Pigott (2015; https://doi.org/10.1002/2015GL066577).
    """

    def __init__(self):
        def _pressure_Pigott_Ni(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Total pressure
            P = P0 + Pth - Pth0

            return P

        _params_Pigott_Ni = {
            "V_0": 6.5790e-06,
            "K_0": 201.0e9,
            "Kprime_0": 4.4,
            "Debye_0": 415.0,
            "grueneisen_0": 1.98,
            "q_0": 1.3,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Pigott_Ni, "pressure", _params_Pigott_Ni)
