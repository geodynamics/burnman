# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from burnman.classes.calibrant import Calibrant
from scipy.interpolate import RegularGridInterpolator

"""
Tsuchiya_2003
^^^^^^^^^^^^^
"""


class Au(Calibrant):
    """
    The Au pressure standard reported by
    Tsuchiya (2003; https://doi.org/10.1029/2003JB002446).
    """

    def __init__(self):

        grid_compressions = np.linspace(0.0, 0.34, 18)
        grid_temperatures = np.array([300.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0])
        grid_pressures = np.array(
            [
                [0.00, 1.52, 5.35, 9.19, 13.04, 16.88],
                [3.55, 5.04, 8.78, 12.54, 16.29, 20.05],
                [7.68, 9.13, 12.79, 16.45, 20.12, 23.79],
                [12.42, 13.83, 17.40, 20.98, 24.56, 28.14],
                [17.86, 19.23, 22.71, 26.20, 29.70, 33.19],
                [24.12, 25.46, 28.85, 32.25, 35.66, 39.07],
                [31.30, 32.60, 35.90, 39.22, 42.54, 45.86],
                [39.52, 40.78, 43.99, 47.22, 50.45, 53.68],
                [48.94, 50.17, 53.29, 56.43, 59.58, 62.72],
                [59.76, 60.95, 63.98, 67.03, 70.09, 73.15],
                [72.11, 73.26, 76.21, 79.18, 82.14, 85.11],
                [86.36, 87.48, 90.34, 93.22, 96.10, 98.98],
                [102.65, 103.73, 106.50, 109.29, 112.08, 114.88],
                [121.38, 122.42, 125.10, 127.80, 130.51, 133.21],
                [142.98, 143.99, 146.58, 149.19, 151.81, 154.43],
                [167.77, 168.74, 171.24, 173.77, 176.30, 178.83],
                [196.48, 197.41, 199.83, 202.26, 204.70, 207.15],
                [229.56, 230.45, 232.78, 235.13, 237.49, 239.84],
            ]
        )

        self.interpolate_pressure = RegularGridInterpolator(
            (grid_compressions, grid_temperatures),
            grid_pressures,
            bounds_error=False,
            fill_value=None,
            method="cubic",
        )

        def _pressure(volume, temperature, params):
            compression = 1.0 - volume / params["V_0"]
            return self.interpolate_pressure([compression, temperature])[0] * 1.0e9

        Calibrant.__init__(self, _pressure, "pressure", {"V_0": 10.207e-06})
