# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from burnman.eos.birch_murnaghan import BirchMurnaghanBase as BM3
from burnman.classes.calibrant import Calibrant
import numpy as np

"""
Zhang_1999
^^^^^^^^^^
"""


class Fe_bcc(Calibrant):
    """
    The BCC Fe pressure standard reported by
    Zhang et al. (1999; https://dx.doi.org/10.1007/s002690050178).
    """

    def __init__(self):
        def _pressure_Zhang_bccFe(volume, temperature, params):

            # Isothermal pressure (GPa)
            pressure_model = BM3()
            P0 = pressure_model.pressure(params["T_0"], volume, params)

            # Thermal pressure
            Pth = (
                params["aK(V,T)"] + params["dK_dT_V"] * np.log(params["V_0"] / volume)
            ) * (temperature - params["T_0"])

            # Total pressure
            P = P0 + Pth

            return P

        _params_Zhang_bccFe = {
            "V_0": 7.08384e-06,
            "K_0": 155.0e9,
            "Kprime_0": 5.3,
            "aK(V,T)": 0.00648e9,
            "dK_dT_V": -0.022e9,
            "n": 1.0,
            "T_0": 300.0,
            "P_0": 0.0,
            "Z": 2.0,
        }

        Calibrant.__init__(self, _pressure_Zhang_bccFe, "pressure", _params_Zhang_bccFe)
