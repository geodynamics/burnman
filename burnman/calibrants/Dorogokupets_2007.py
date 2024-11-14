# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from burnman.classes.calibrant import Calibrant
from burnman.constants import gas_constant
from burnman.eos.vinet import Vinet

"""
Dorogokupets_2007
^^^^^^^^^^^^^^^^^
"""

_materials_data = {
    "Ag": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 10.272e-6,
        "K_0": 99.65e9,
        "Kprime_0": 6.11,
        "theta_B1": 130.6,
        "d_B1": 8.572,
        "m_B1": 0.121,
        "theta_B2": 103.6,
        "d_B2": 5.326,
        "m_B2": 0.449,
        "theta_E1": 111.9,
        "m_E1": 0.766,
        "theta_E2": 189.12,
        "m_E2": 1.664,
        "gamma_0": 2.376,
        "gamma_inf": 1.481,
        "beta": 2.507,
        "a": 6.70e-6,
        "m": 3.44,
        "e": 25.9e-6,
        "f": -1.0,
        "g": 0.666,
        "h": -2.0,
        "H": 15239,
        "S": 0.732,
    },
    "Al": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 9.999e-6,
        "K_0": 72.67e9,
        "Kprime_0": 4.62,
        "theta_B1": 245.8,
        "d_B1": 5.575,
        "m_B1": 0.987,
        "theta_B2": 300.0,  # Not used as m = 0
        "d_B2": 1.0,  # Not used as m = 0
        "m_B2": 0.0,
        "theta_E1": 240.2,
        "m_E1": 1.000,
        "theta_E2": 356.2,
        "m_E2": 1.013,
        "gamma_0": 2.144,
        "gamma_inf": 1.017,
        "beta": 3.942,
        "a": 5.14e-6,
        "m": 3.44,
        "e": 54.1e-6,
        "f": -1.0,
        "g": 1.8,
        "h": -2.0,
        "H": 8679,
        "S": 0.998,
    },
    "Au": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 10.215e-6,
        "K_0": 166.70e9,
        "Kprime_0": 6.00,
        "theta_B1": 95.7,
        "d_B1": 8.290,
        "m_B1": 0.681,
        "theta_B2": 106.4,
        "d_B2": 3.239,
        "m_B2": 0.417,
        "theta_E1": 170.6,
        "m_E1": 1.063,
        "theta_E2": 105.2,
        "m_E2": 0.839,
        "gamma_0": 2.965,
        "gamma_inf": 1.142,
        "beta": 3.030,
        "a": 25.33e-6,
        "m": 3.79,
        "e": 18.92e-6,
        "f": -1.0,
        "g": 0.66,
        "h": -2.0,
        "H": 11.69e3,  # note typo
        "S": 1.067,
    },
    "Cu": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 7.113e-6,
        "K_0": 133.41e9,
        "Kprime_0": 5.37,
        "theta_B1": 123.7,
        "d_B1": 3.776,
        "m_B1": 0.115,
        "theta_B2": 175.4,
        "d_B2": 10.372,
        "m_B2": 0.711,
        "theta_E1": 187.4,
        "m_E1": 0.756,
        "theta_E2": 286.9,
        "m_E2": 1.418,
        "gamma_0": 1.974,
        "gamma_inf": 1.554,
        "beta": 4.647,
        "a": 3.50e-6,
        "m": 3.46,
        "e": 27.698e-6,
        "f": -1.0,
        "g": 0.666,
        "h": -2.0,
        "H": 11687,
        "S": 1.407,
    },
    "Pt": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 9.091e-6,
        "K_0": 276.07e9,
        "Kprime_0": 5.30,
        "theta_B1": 95.2,
        "d_B1": 8.199,
        "m_B1": 0.329,
        "theta_B2": 148.4,
        "d_B2": 4.005,
        "m_B2": 0.383,
        "theta_E1": 214.6,
        "m_E1": 1.211,
        "theta_E2": 140.8,
        "m_E2": 1.077,
        "gamma_0": 2.802,
        "gamma_inf": 1.538,
        "beta": 5.550,
        "a": 160.9e-6,
        "m": 4.06,
        "e": 260.0e-6,
        "f": -1.0,
        "g": 2.4,
        "h": -2.0,
        "H": 32572,
        "S": 0.631,
    },
    "Ta": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 10.851e-6,
        "K_0": 191.39e9,
        "Kprime_0": 3.81,
        "theta_B1": 72.6,
        "d_B1": 5.536,
        "m_B1": 0.117,
        "theta_B2": 101.8,
        "d_B2": 24.513,
        "m_B2": 0.396,
        "theta_E1": 144.0,
        "m_E1": 1.118,
        "theta_E2": 214.9,
        "m_E2": 1.369,
        "gamma_0": 1.714,
        "gamma_inf": 1.241,
        "beta": 6.825,
        "a": 61.9e-6,
        "m": 4.00,
        "e": 167.0e-6,
        "f": -1.0,
        "g": 1.3,
        "h": -2.0,
        "H": 36278,
        "S": 4.910,
    },
    "W": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 9.545e-6,
        "K_0": 306.00e9,
        "Kprime_0": 4.17,
        "theta_B1": 182.8,
        "d_B1": 13.270,
        "m_B1": 0.513,
        "theta_B2": 172.5,
        "d_B2": 3.305,
        "m_B2": 0.174,
        "theta_E1": 287.6,
        "m_E1": 1.166,
        "theta_E2": 213.8,
        "m_E2": 1.145,
        "gamma_0": 1.553,
        "gamma_inf": 0.694,
        "beta": 3.698,
        "a": -39.3e-6,
        "m": 2.67,
        "e": 40.4e-6,
        "f": -1.0,
        "g": 0.2,
        "h": -2.0,
        "H": 14714,
        "S": 0.672,
    },
    "MgO": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 2.0,
        "V_0": 11.248e-6,
        "K_0": 160.31e9,
        "Kprime_0": 4.18,
        "theta_B1": 447.3,
        "d_B1": 11.248,
        "m_B1": 1.429,
        "theta_B2": 384.0,
        "d_B2": 3.593,
        "m_B2": 0.276,
        "theta_E1": 703.8,
        "m_E1": 2.570,
        "theta_E2": 466.0,
        "m_E2": 1.725,
        "gamma_0": 1.522,
        "gamma_inf": 1.111,
        "beta": 4.509,
        "a": 13.56e-6,
        "m": 5.23,
        "e": 0,
        "f": 0.0,
        "g": 0,
        "h": 0.0,
        "H": 0,
        "S": 0,
    },
    "Diamond": {
        "E_0": 0.0,
        "P_0": 0.0,
        "T_0": 298.15,
        "n": 1.0,
        "V_0": 3.417e-6,
        "K_0": 443.16e9,
        "Kprime_0": 3.777,
        "theta_B1": 1202.1,
        "d_B1": 9.604,
        "m_B1": 1.163,
        "theta_B2": 1135.1,
        "d_B2": 3.380,
        "m_B2": 0.218,
        "theta_E1": 1687.2,
        "m_E1": 1.396,
        "theta_E2": 1033.7,
        "m_E2": 0.223,
        "gamma_0": 0.820,
        "gamma_inf": 0.615,
        "beta": 10.121,
        "a": -23.85e-6,
        "m": 1.22,
        "e": 0,
        "f": 0.0,
        "g": 0,
        "h": 0.0,
        "H": 0,
        "S": 0,
    },
}


do_cold_eos = Vinet()


def _F_qhB(m, d, theta, temperature):
    # Equation 8a and equations just after Equation 7
    g = d * np.log(1.0 + theta / (temperature * d))
    b = 1.0 / (np.exp(g) - 1.0)
    return m * ((d - 1.0) / (2.0 * d) * theta - temperature * np.log(1.0 + b))


def _F_qhE(m, theta, temperature):
    # Equation 8b
    return m * (theta / 2.0 + temperature * np.log(1.0 - np.exp(-theta / temperature)))


def _F_anh_factor(m, theta, temperature):
    # Equation 12
    zeta = np.exp(theta / temperature)
    f = np.power(0.5 * theta + theta / (zeta - 1.0), 2.0) + 2.0 * theta * theta * (
        zeta / np.power(zeta - 1.0, 2.0)
    )
    return m * f


def _F_th(temperature, x, params):
    # Equation 10
    F_theta = np.power(x, -params["gamma_inf"]) * np.exp(
        (params["gamma_0"] - params["gamma_inf"])
        / params["beta"]
        * (1.0 - np.power(x, params["beta"]))
    )
    theta_B1 = params["theta_B1"] * F_theta
    theta_B2 = params["theta_B2"] * F_theta
    theta_E1 = params["theta_E1"] * F_theta
    theta_E2 = params["theta_E2"] * F_theta

    # Equation 8
    F_qh = (
        _F_qhB(params["m_B1"], params["d_B1"], theta_B1, temperature)
        + _F_qhB(params["m_B2"], params["d_B2"], theta_B2, temperature)
        + _F_qhE(params["m_E1"], theta_E1, temperature)
        + _F_qhE(params["m_E2"], theta_E2, temperature)
    )

    # Equation 12
    F_anh = (
        params["a"]
        * np.power(x, params["m"])
        / 6.0
        * (
            _F_anh_factor(params["m_B1"], theta_B1, temperature)
            + _F_anh_factor(params["m_B2"], theta_B2, temperature)
            + _F_anh_factor(params["m_E1"], theta_E1, temperature)
            + _F_anh_factor(params["m_E2"], theta_E2, temperature)
        )
    )

    # Equation 13
    F_el = (
        -1.5
        * params["n"]
        * params["e"]
        * np.power(x, params["g"])
        * temperature
        * temperature
    )

    # Equation 14
    F_def = (
        -1.5
        * params["n"]
        * temperature
        * np.exp(
            params["S"] * np.power(x, params["f"])
            - params["H"] * np.power(x, params["h"]) / temperature
        )
    )
    return F_qh + F_anh + F_el + F_def


def _helmholtz_energy(volume, temperature, params):
    # Equation 6
    E = do_cold_eos.molar_internal_energy(0.0, 0.0, volume, params)

    # Just after Equation 6b
    Vrel = volume / params["V_0"]
    Fth = _F_th(temperature, Vrel, params)
    Fth0 = _F_th(params["T_0"], Vrel, params)

    # Equation 5
    return E + gas_constant * (Fth - Fth0)


def _pressure(volume, temperature, params):
    dV = params["V_0"] * 1.0e-7

    F1 = _helmholtz_energy(volume + dV / 2.0, temperature, params)
    F0 = _helmholtz_energy(volume - dV / 2.0, temperature, params)
    return -(F1 - F0) / dV


class Ag(Calibrant):
    """
    The Ag pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Ag"])


class Al(Calibrant):
    """
    The Al pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Al"])


class Au(Calibrant):
    """
    The Au pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Au"])


class Cu(Calibrant):
    """
    The Cu pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Cu"])


class diamond(Calibrant):
    """
    The diamond pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Diamond"])


class MgO(Calibrant):
    """
    The MgO pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["MgO"])


class Pt(Calibrant):
    """
    The Pt pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Pt"])


class Ta(Calibrant):
    """
    The Ta pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["Ta"])


class W(Calibrant):
    """
    The W pressure standard reported by
    Dorogokupets and Oganov (2007; https://doi.org/10.1103/PhysRevB.75.024115).
    """

    def __init__(self):
        Calibrant.__init__(self, _pressure, "pressure", _materials_data["W"])
