# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Armentrout_2015
^^^^^^^^^^^^^^^
"""


class Co_fcc(Calibrant):
    """
    The FCC Cobalt pressure standard reported by Armentrout
    (2015; https://doi.org/10.1063/1.4935087).
    """

    def __init__(self):
        def _pressure_Armentrout_fccCo(volume, temperature, params):

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

        _params_Armentrout_fccCo = {
            "V_0": 6.7529e-06,
            "K_0": 196.0e9,
            "Kprime_0": 4.7,
            "Debye_0": 385.0,
            "grueneisen_0": 2.0,
            "q_0": 1.3,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(
            self, _pressure_Armentrout_fccCo, "pressure", _params_Armentrout_fccCo
        )
