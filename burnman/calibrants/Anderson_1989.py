# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.classes.calibrant import Calibrant
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume

"""
Anderson_1989
^^^^^^^^^^^^^
"""


class Au(Calibrant):
    """
    The Au pressure standard reported by
    Anderson et al. (1989; https://doi.org/10.1063/1.342969).
    """

    def __init__(self):
        def _pressure(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth_rel = (7.14e6 + params["dKTdT"] * np.log(params["V_0"] / volume)) * (
                temperature - 300.0
            )

            return P0 + Pth_rel

        Z = 4.0
        _params = {
            "V_0": molar_volume_from_unit_cell_volume(67.850, Z),
            "K_0": 166.65e9,
            "Kprime_0": 5.4823,
            "dKTdT": -11.5e6,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": Z,
        }

        Calibrant.__init__(self, _pressure, "pressure", _params)
