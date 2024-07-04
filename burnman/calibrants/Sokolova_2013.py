# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Sokolova_2013
^^^^^^^^^^^^^
"""


class Mo(Calibrant):
    """
    The Mo pressure standard reported by Sokolova et al. (2013),
    as reported by Huang et al. (2016; https://doi.org/10.1038/srep19923).
    """

    def __init__(self):
        def _pressure_Sokolova_Mo(volume, temperature, params):

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

        _params_Sokolova_Mo = {
            "V_0": 9.37647e-06,
            "K_0": 249.0e9,
            "Kprime_0": 4.47,
            "Debye_0": 470.0,
            "grueneisen_0": 1.98,
            "q_0": 1.99,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Sokolova_Mo, "pressure", _params_Sokolova_Mo)
