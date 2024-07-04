# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant

"""
Dewaele_2008
^^^^^^^^^^^^
"""


class Ni(Calibrant):
    """
    The Ni pressure standard reported by
    Dewaele (2008; Table IV; https://doi.org/10.1103/PhysRevB.78.104102).
    """

    def __init__(self):
        def _pressure_Dewaele_Ni(volume, temperature, params):

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

        _params_Dewaele_Ni = {
            "V_0": 6.5792e-06,
            "K_0": 176.7e9,
            "Kprime_0": 5.23,
            "Debye_0": 415.0,
            "grueneisen_0": 1.98,
            "q_0": 1.3,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Dewaele_Ni, "pressure", _params_Dewaele_Ni)
