# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Matsui_2009
^^^^^^^^^^^
"""


class Pt(Calibrant):
    """
    The Pt pressure standard reported by
    Matsui (2009; https://doi.org/10.1063/1.3054331).
    """

    def __init__(self):
        def _pressure_Matsui_Pt(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Electronic pressure
            Pel = (
                1.1916e-15 * temperature**4.0
                - 1.4551e-11 * temperature**3.0
                + 1.6209e-07 * temperature**2.0
                + 1.8269e-4 * temperature
                - 0.069
            ) * 1.0e09

            # Total pressure
            P = P0 + Pth - Pth0 + Pel

            return P

        _params_Matsui_Pt = {
            "V_0": 9.0904e-06,
            "K_0": 273.0e9,
            "Kprime_0": 5.2,
            "Debye_0": 230.0,
            "grueneisen_0": 2.7,
            "q_0": 1.1,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Matsui_Pt, "pressure", _params_Matsui_Pt)
