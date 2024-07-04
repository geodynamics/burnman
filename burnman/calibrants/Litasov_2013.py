# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.vinet import Vinet
from burnman.eos.mie_grueneisen_debye import MGDBase
from burnman.classes.calibrant import Calibrant


"""
Litasov_2013
^^^^^^^^^^^^
"""


class Mo_bcc(Calibrant):
    """
    The BCC Mo pressure standard reported by
    Litasov et al. (2013; https://doi.org/10.1063/1.4794127).
    """

    def __init__(self):
        def _pressure_Litasov_Mo(volume, temperature, params):

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

        _params_Litasov_Mo = {
            "V_0": 9.37647e-06,
            "K_0": 260.0e9,
            "Kprime_0": 4.21,
            "Debye_0": 470.0,
            "grueneisen_0": 2.03,
            "q_0": 0.24,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Litasov_Mo, "pressure", _params_Litasov_Mo)


class W_bcc(Calibrant):
    """
    The BCC W pressure standard reported by
    Litasov et al. (2013; https://doi.org/10.1063/1.4799018).
    """

    def __init__(self):
        def _pressure_Litasov_W(volume, temperature, params):

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

        _params_Litasov_W = {
            "V_0": 9.5481e-06,
            "K_0": 308.0e9,
            "Kprime_0": 4.2,
            "Debye_0": 388.0,
            "grueneisen_0": 1.81,
            "q_0": 0.71,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Litasov_W, "pressure", _params_Litasov_W)
