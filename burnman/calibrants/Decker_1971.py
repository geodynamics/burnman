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
from burnman.eos.birch_murnaghan_4th import birch_murnaghan_fourth


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
            P0 = birch_murnaghan_fourth(params["V_0"] / volume, params)

            # Thermal pressure
            thermal_model = MGDBase()
            Pth0 = thermal_model._thermal_pressure(params["T_0"], volume, params)
            Pth = thermal_model._thermal_pressure(temperature, volume, params)

            return P0 + Pth - Pth0

        _params_Decker_NaCl = {
            "V_0": 2.7013e-05,
            "K_0": 23.7e9,
            "Kprime_0": 4.91,
            "Kprime_prime_0": -0.267e-9,
            "Debye_0": 279.0,
            "grueneisen_0": 1.59,
            "q_0": 0.93,
            "n": 2.0,
            "T_0": 298.15,
            "P_0": 0.0,
            "Z": 4.0,
        }

        Calibrant.__init__(self, _pressure_Decker_NaCl, "pressure", _params_Decker_NaCl)
