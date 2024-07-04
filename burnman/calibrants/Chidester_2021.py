# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Chidester_2021
^^^^^^^^^^^^^^
"""


class KCl_B2(Calibrant):
    """
    The B2 KCl pressure standard reported by
    Chidester et al. (2021; https://doi.org/10.1103/PhysRevB.104.094107).
    """

    def __init__(self):
        def _pressure_Chidester_KCl(volume, temperature, params):

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

        _params_Chidester_KCl = {
            "V_0": 3.1200e-05,
            "K_0": 24.0e9,
            "Kprime_0": 4.56,
            "Debye_0": 235.0,
            "grueneisen_0": 2.9,
            "q_0": 1.0,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(
            self, _pressure_Chidester_KCl, "pressure", _params_Chidester_KCl
        )
