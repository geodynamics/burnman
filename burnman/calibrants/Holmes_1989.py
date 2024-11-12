# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
Holmes_1989
^^^^^^^^^^^
"""

from burnman.classes.calibrant import Calibrant
import numpy as np
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume


class Pt(Calibrant):
    """
    The Pt pressure standard reported by
    Holmes et al. (1989; https://doi.org/10.1063/1.344177).
    """

    def __init__(self):
        def _pressure(volume, temperature, params):
            X = np.power(volume / params["V_0"], 1.0 / 3.0)
            P_300 = (
                3.0
                * params["beta_T"]
                * (1.0 - X)
                / (X * X)
                * np.exp(params["eta"] * (1.0 - X))
            )

            return P_300 + params["alpha_T"] * params["beta_T"] * (temperature - 300.0)

        Z = 4.0
        _params = {
            "V_0": molar_volume_from_unit_cell_volume(60.38, Z),
            "beta_T": 798.31e9 / 3.0,
            "eta": 7.2119,
            "beta_prime_T": (7.2119 / 1.5) + 1.0,
            "alpha_T": 2.61e-5,
            "Z": Z,
        }

        Calibrant.__init__(self, _pressure, "pressure", _params)
