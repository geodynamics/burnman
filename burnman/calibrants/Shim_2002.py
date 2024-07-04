# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Shim_2002
^^^^^^^^^
"""


class Au(Calibrant):
    """
    The Au pressure standard reported by
    Shim et al. (2002; https://doi.org/10.1016/S0012-821X(02)00917-2).
    """

    def __init__(self):
        def _pressure_Shim_Au(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Electronic pressure
            Pel = (
                1.4664e-16 * temperature**4.0
                - 8.0179e-13 * temperature**3.0
                + 1.6205e-8 * temperature**2.0
                - 5.4573e-6 * temperature
                - 7.8273e-5
            ) * 1.0e09

            # Total pressure
            P = P0 + Pth - Pth0 + Pel

            return P

        _params_Shim_Au = {
            "V_0": 1.0215e-05,
            "K_0": 167.0e9,
            "Kprime_0": 5.0,
            "Debye_0": 170.0,
            "grueneisen_0": 2.97,
            "q_0": 1.0,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Shim_Au, "pressure", _params_Shim_Au)
