# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
Decker_1971
^^^^^^^^^^^
"""

from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


class NaCl_B1(Calibrant):
    """
    The NaCl (B1 structured) pressure standard reported by
    Decker (1971; https://doi.org/10.1063/1.1660714).

    .. note:: This calibrant is not exactly the same as that proposed by Decker.
        The cold compression curve has here been approximated by a 4th order
        Birch-Murnaghan EoS, as described in
        Matsui et al. (2012; https://doi.org/10.2138/am.2012.4136).
    """

    def __init__(self):
        def _pressure_Decker_NaCl(volume, temperature, params):

            # Isothermal pressure (GPa)
            a = (3 / 2) * (params["Kprime_0"] - 4)
            b = (
                9 * params["K_0"] * params["Kprime_prime_0"]
                + 9 * params["Kprime_0"] ** 2
                - 63 * params["Kprime_0"]
                + 143
            )
            f = 0.5 * ((volume / params["V_0"]) ** (-2 / 3) - 1)
            K_T = (
                params["K_0"]
                * ((1 + 2 * f) ** (5 / 2))
                * (1 + (7 + 2 * a) * f + (9 * a + 3 * b) * f**2 + 11 * b * f**3)
            )
            P0 = 3 * f * params["K_0"] * (1 + 2 * f) ** (5 / 2) * (1 + a * f + b * f**2)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            # Total pressure
            P = (P0 * 1e9) + Pth - Pth0

            return P

        _params_Decker_NaCl = {
            "V_0": 2.7013e-05,
            "K_0": 23.7,
            "Kprime_0": 5.04,  # 4.91 in Matsui (2012), however 5.04 is required to reproduce values in Table 4.
            "Kprime_prime_0": -0.267,
            "Debye_0": 279.0,
            "grueneisen_0": 1.59,
            "q_0": 0.93,
            "n": 2.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Decker_NaCl, "pressure", _params_Decker_NaCl)
