# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Fei_2016
^^^^^^^^
"""


class Fe_hcp(Calibrant):
    """
    The HCP Fe pressure standard reported by
    Fei et al. (2016; https://doi.org/10.1002/2016GL069456).
    """

    def __init__(self):
        def _pressure_Fei_hcpFe(volume, temperature, params):

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

        _params_Fei_hcpFe = {
            "V_0": 6.8050e-6,
            "K_0": 172.7e9,
            "Kprime_0": 4.79,
            "Debye_0": 422.0,
            "grueneisen_0": 1.74,
            "q_0": 0.78,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Fei_hcpFe, "pressure", _params_Fei_hcpFe)
