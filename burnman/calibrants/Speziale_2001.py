# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Speziale_2001
^^^^^^^^^^^^^
"""


class MgO(Calibrant):
    """
    The MgO pressure standard reported by
    Speziale (2001; https://doi.org/10.1029/2000JB900318).
    """

    def __init__(self):
        def _pressure_Speziale_MgO(volume, temperature, params):

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

        _params_Speziale_MgO = {
            "V_0": 1.12463e-05,
            "K_0": 160.2e9,
            "Kprime_0": 3.99,
            "Debye_0": 773.0,
            "grueneisen_0": 1.524,
            "q_0": 1.65,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(
            self, _pressure_Speziale_MgO, "pressure", _params_Speziale_MgO
        )
