# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Campbell_2009
^^^^^^^^^^^^^
"""


class Ni_fcc(Calibrant):
    """
    The FCC Ni pressure standard reported by
    Campbell (2009; https://doi.org/10.1016/j.epsl.2009.07.022).
    """

    def __init__(self):
        def _pressure_Campbell_Ni(volume, temperature, params):

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

        _params_Campbell_Ni = {
            "V_0": 6.5870e-06,
            "K_0": 179.0e9,
            "Kprime_0": 4.3,
            "Debye_0": 415.0,
            "grueneisen_0": 2.5,
            "q_0": 1.0,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Campbell_Ni, "pressure", _params_Campbell_Ni)
