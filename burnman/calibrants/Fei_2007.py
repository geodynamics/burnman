# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume

"""
Fei_2007
^^^^^^^^
"""


class Pt(Calibrant):
    """
    The Pt pressure standard reported by
    Fei et al. (2007; https://doi.org/10.1073/pnas.0609013104).
    """

    def __init__(self):
        def _pressure_Fei_Pt(volume, temperature, params):

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

        Z = 4.0
        _params_Fei_Pt = {
            "V_0": molar_volume_from_unit_cell_volume(60.38, Z),
            "K_0": 277.0e9,
            "Kprime_0": 5.08,
            "Debye_0": 230.0,
            "grueneisen_0": 2.72,
            "q_0": 0.5,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": Z,
        }

        Calibrant.__init__(self, _pressure_Fei_Pt, "pressure", _params_Fei_Pt)


class Au(Calibrant):
    """
    The Au pressure standard reported by
    Fei et al. (2007; https://doi.org/10.1073/pnas.0609013104).
    """

    def __init__(self):
        def _pressure_Fei_Au(volume, temperature, params):

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

        Z = 4.0
        _params_Fei_Au = {
            "V_0": molar_volume_from_unit_cell_volume(67.850, Z),
            "K_0": 167.0e9,
            "Kprime_0": 6.0,
            "Debye_0": 170.0,
            "grueneisen_0": 2.97,
            "q_0": 0.6,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": Z,
        }

        Calibrant.__init__(self, _pressure_Fei_Au, "pressure", _params_Fei_Au)
