# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.vinet import Vinet
from burnman.eos import debye
from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Tange_2009
^^^^^^^^^^
"""


class MgO_BM3(Calibrant):
    """
    The MgO pressure standard reported by
    Tange et al. (2009; https://doi.org/10.1029/2008JB005813).
    BM3 version.
    """

    def __init__(self):
        def _pressure_Tange_MgO(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            gra = params["grueneisen_0"] * (
                1.0 + params["a"] * ((volume / params["V_0"]) ** params["b"] - 1.0)
            )
            debye_Ta = (
                params["Debye_0"]
                * (volume / params["V_0"])
                ** (-params["grueneisen_0"] * (1.0 - params["a"]))
                * np.exp(-(gra - params["grueneisen_0"]) / params["b"])
            )

            E_th = debye.thermal_energy(temperature, debye_Ta, params["n"])
            E_th_ref = debye.thermal_energy(
                params["T_0"], debye_Ta, params["n"]
            )  # thermal energy at reference temperature

            Pth = gra * (E_th - E_th_ref) / volume

            # Total pressure
            P = P0 + Pth

            return P

        _params_Tange_MgO = {
            "V_0": 1.12463e-05,
            "K_0": 160.64e9,
            "Kprime_0": 4.221,
            "Debye_0": 761.0,
            "grueneisen_0": 1.431,
            "a": 0.29,
            "b": 3.5,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Tange_MgO, "pressure", _params_Tange_MgO)


class MgO_Vinet(Calibrant):
    """
    The MgO pressure standard reported by
    Tange et al. (2009; https://doi.org/10.1029/2008JB005813).
    Vinet version.
    """

    def __init__(self):
        def _pressure_Tange_MgO(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = Vinet()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            gra = params["grueneisen_0"] * (
                1.0 + params["a"] * ((volume / params["V_0"]) ** params["b"] - 1.0)
            )
            debye_Ta = (
                params["Debye_0"]
                * (volume / params["V_0"])
                ** (-params["grueneisen_0"] * (1.0 - params["a"]))
                * np.exp(-(gra - params["grueneisen_0"]) / params["b"])
            )

            E_th = debye.thermal_energy(temperature, debye_Ta, params["n"])
            E_th_ref = debye.thermal_energy(
                params["T_0"], debye_Ta, params["n"]
            )  # thermal energy at reference temperature

            Pth = gra * (E_th - E_th_ref) / volume

            # Total pressure
            P = P0 + Pth

            return P

        _params_Tange_MgO = {
            "V_0": 1.12463e-05,
            "K_0": 160.63e9,
            "Kprime_0": 4.367,
            "Debye_0": 761.0,
            "grueneisen_0": 1.442,
            "a": 0.138,
            "b": 5.4,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Tange_MgO, "pressure", _params_Tange_MgO)
