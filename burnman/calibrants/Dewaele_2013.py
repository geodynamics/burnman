# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
import numpy as np


"""
Dewaele_2013
^^^^^^^^^^^^
"""


class Al2O3_corundum(Calibrant):
    """
    The alpha-Al2O3 (corundum) pressure standard reported by
    Dewaele and Torrent (2013; https://doi.org/10.1103/PhysRevB.88.064107).
    The thermal model is from Shi (2022).
    """

    def __init__(self):
        def _pressure_Dewaele_Al2O3(volume, temperature, params):

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

        _params_Dewaele_Al2O3 = {
            "V_0": 2.5594e-5,
            "K_0": 254.1e9,
            "Kprime_0": 4.0,
            "Debye_0": 1100.0,
            "grueneisen_0": 1.32,
            "q_0": 0.8,
            "n": 5.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 6.0,
        }

        Calibrant.__init__(
            self, _pressure_Dewaele_Al2O3, "pressure", _params_Dewaele_Al2O3
        )
