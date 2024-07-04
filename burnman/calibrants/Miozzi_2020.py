# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Miozzi_2020
^^^^^^^^^^^
"""


class Fe_hcp(Calibrant):
    """
    The HCP Fe pressure standard reported by
    Miozzi (2020; https://doi.org/10.3390/min10020100).
    """

    def __init__(self):
        def _pressure_Miozzi_hcpFe(volume, temperature, params):

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

        _params_Miozzi_hcpFe = {
            "V_0": 6.86825e-06,
            "K_0": 129.0e9,
            "Kprime_0": 6.24,
            "Debye_0": 420.0,
            "grueneisen_0": 1.11,
            "q_0": 0.3,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(
            self, _pressure_Miozzi_hcpFe, "pressure", _params_Miozzi_hcpFe
        )
