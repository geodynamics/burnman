# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant

"""
Tateno_2019
^^^^^^^^^^^
"""


class KCl_B2(Calibrant):
    """
    The B2 KCl pressure standard reported by
    Tateno (2019; https://doi.org/10.2138/am-2019-6779).
    """

    def __init__(self):
        def _pressure_Tateno_KCl(volume, temperature, params):

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

        _params_Tateno_KCl = {
            "V_0": 3.28206e-05,
            "K_0": 17.4e9,
            "Kprime_0": 5.77,
            "Debye_0": 235.0,
            "grueneisen_0": 1.8,
            "q_0": 0.7,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 1.0,
        }

        Calibrant.__init__(self, _pressure_Tateno_KCl, "pressure", _params_Tateno_KCl)
