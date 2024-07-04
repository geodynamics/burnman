# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Dorogokupets_2017
^^^^^^^^^^^^^^^^^
"""


class Fe_bcc(Calibrant):
    """
    The BCC Fe pressure standard reported by
    Dorogokupets (2017; https://doi.org/10.1038/srep41863).
    """

    def __init__(self):
        def _pressure_Dorogokupets_bccFe(volume, temperature, params):

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

        _params_Dorogokupets_bccFe = {
            "V_0": 7.09197e-06,
            "K_0": 164.0e9,
            "Kprime_0": 5.5,
            "Debye_0": 303.0,
            "grueneisen_0": 1.736,
            "q_0": 1.125,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(
            self, _pressure_Dorogokupets_bccFe, "pressure", _params_Dorogokupets_bccFe
        )


class Fe_fcc(Calibrant):
    """
    The FCC Fe pressure standard reported by
    Dorogokupets (2017; https://doi.org/10.1038/srep41863).
    """

    def __init__(self):
        def _pressure_Dorogokupets_fccFe(volume, temperature, params):

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

        _params_Dorogokupets_fccFe = {
            "V_0": 6.9285e-06,
            "K_0": 146.2e9,
            "Kprime_0": 4.67,
            "Debye_0": 222.5,
            "grueneisen_0": 2.203,
            "q_0": 0.01,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(
            self, _pressure_Dorogokupets_fccFe, "pressure", _params_Dorogokupets_fccFe
        )


class Fe_hcp(Calibrant):
    """
    The HCP Fe pressure standard reported by
    Dorogokupets (2017; https://doi.org/10.1038/srep41863).
    """

    def __init__(self):
        def _pressure_Dorogokupets_hcpFe(volume, temperature, params):

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

        _params_Dorogokupets_hcpFe = {
            "V_0": 6.81706e-6,
            "K_0": 148.0e9,
            "Kprime_0": 5.86,
            "Debye_0": 227.0,
            "grueneisen_0": 2.2,
            "q_0": 0.01,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(
            self, _pressure_Dorogokupets_hcpFe, "pressure", _params_Dorogokupets_hcpFe
        )
