# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import scipy.optimize as opt
from . import equation_of_state as eos
from ..utils.math import bracket
import warnings


def bulk_modulus_third_order(volume, params):
    """
    Bulk modulus for the third order Birch-Murnaghan equation of state.

    :param volume: Volume of the material in the same units as
    the reference volume.
    :type volume: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Bulk modulus in the same units as the reference bulk modulus.
    :rtype: float
    """

    x = params["V_0"] / volume
    f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0)

    K = pow(1.0 + 2.0 * f, 5.0 / 2.0) * (
        params["K_0"]
        + (3.0 * params["K_0"] * params["Kprime_0"] - 5 * params["K_0"]) * f
        + 27.0
        / 2.0
        * (params["K_0"] * params["Kprime_0"] - 4.0 * params["K_0"])
        * f
        * f
    )
    return K


def pressure_third_order(invVrel, params):
    """
    Pressure for the third order Birch-Murnaghan equation of state.

    :param invVrel: Reference volume divided by the volume
    :type invVrel: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Pressure in the same units that are supplied
    for the reference bulk modulus (params['K_0']).
    :rtype: float
    """

    return (
        3.0
        * params["K_0"]
        / 2.0
        * (pow(invVrel, 7.0 / 3.0) - pow(invVrel, 5.0 / 3.0))
        * (1.0 - 0.75 * (4.0 - params["Kprime_0"]) * (pow(invVrel, 2.0 / 3.0) - 1.0))
        + params["P_0"]
    )


def volume_third_order(pressure, params):
    """
    Volume for the third order Birch-Murnaghan equation of state.

    :param pressure: Pressure in the same units that are supplied
    for the reference bulk modulus (params['K_0']).
    :type pressure: float
    :param params: Parameter dictionary
    :type params: dict

    :return: Molar volume in the same units as the reference volume.
    :rtype: float
    """

    def delta_pressure(volume):
        return pressure_third_order(params["V_0"] / volume, params) - pressure

    try:
        sol = bracket(delta_pressure, params["V_0"], 1.0e-2 * params["V_0"])
    except ValueError:
        raise ValueError(
            "Cannot find a volume, perhaps you are outside of the "
            "range of validity for the equation of state?"
        )
    return opt.brentq(delta_pressure, sol[0], sol[1])


def bulk_modulus_fourth_order(volume, params):
    """
    Bulk modulus for the fourth order Birch-Murnaghan equation of state.

    :param volume: Volume of the material in the same units as
    the reference volume.
    :type volume: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Bulk modulus in the same units as the reference bulk modulus.
    :rtype: float
    """

    invVrel = params["V_0"] / volume
    f = 0.5 * (pow(invVrel, 2.0 / 3.0) - 1.0)

    Xi = (3.0 / 4.0) * (4.0 - params["Kprime_0"])
    Zeta = (3.0 / 8.0) * (
        (params["K_0"] * params["Kprime_prime_0"])
        + params["Kprime_0"] * (params["Kprime_0"] - 7.0)
        + 143.0 / 9.0
    )

    K = (
        5.0
        * f
        * pow((1.0 + 2.0 * f), 5.0 / 2.0)
        * params["K_0"]
        * (1.0 - (2.0 * Xi * f) + (4.0 * Zeta * pow(f, 2.0)))
    ) + (
        pow(1.0 + (2.0 * f), 7.0 / 2.0)
        * params["K_0"]
        * (1.0 - (4.0 * Xi * f) + (12.0 * Zeta * pow(f, 2.0)))
    )

    return K


def pressure_fourth_order(invVrel, params):
    """
    Pressure for the fourth order Birch-Murnaghan equation of state.

    :param invVrel: Reference volume divided by the volume
    :type invVrel: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Pressure in the same units that are supplied
    for the reference bulk modulus (params['K_0']).
    :rtype: float
    """

    f = 0.5 * (pow(invVrel, 2.0 / 3.0) - 1.0)
    Xi = (3.0 / 4.0) * (4.0 - params["Kprime_0"])
    Zeta = (3.0 / 8.0) * (
        (params["K_0"] * params["Kprime_prime_0"])
        + params["Kprime_0"] * (params["Kprime_0"] - 7.0)
        + 143.0 / 9.0
    )

    return (
        3.0
        * f
        * pow(1.0 + 2.0 * f, 5.0 / 2.0)
        * params["K_0"]
        * (1.0 - (2.0 * Xi * f) + (4.0 * Zeta * pow(f, 2.0)))
        + params["P_0"]
    )


def volume_fourth_order(pressure, params):
    """
    Volume for the fourth order Birch-Murnaghan equation of state.

    :param pressure: Pressure in the same units that are supplied
    for the reference bulk modulus (params['K_0']).
    :type pressure: float
    :param params: Parameter dictionary
    :type params: dict

    :return: Molar volume in the same units as the reference volume.
    :rtype: float
    """

    def delta_pressure(x):
        return pressure_fourth_order(params["V_0"] / x, params) - pressure

    try:
        sol = bracket(delta_pressure, params["V_0"], 1.0e-2 * params["V_0"])
    except ValueError:
        raise ValueError(
            "Cannot find a volume, perhaps you are outside of "
            "the range of validity for the equation of state?"
        )
    return opt.brentq(delta_pressure, sol[0], sol[1])


def shear_modulus_second_order(volume, params):
    """
    Shear modulus for the second order Birch Murnaghan equation of state
    (i.e. expanded to 2nd order in strain).

    :param volume: Molar volume in the same units as the reference volume.
    :type volume: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Shear modulus in the same units as the reference shear modulus.
    :rtype: float
    """

    x = params["V_0"] / volume
    G = (
        params["G_0"]
        * pow(x, 5.0 / 3.0)
        * (
            1.0
            - 0.5
            * (pow(x, 2.0 / 3.0) - 1.0)
            * (5.0 - 3.0 * params["Gprime_0"] * params["K_0"] / params["G_0"])
        )
    )
    return G


def shear_modulus_third_order(volume, params):
    """
    Shear modulus for the third order Birch Murnaghan equation of state.

    :param volume: Molar volume in the same units as the reference volume.
    :type volume: float
    :param params: Parameter dictionary
    :type params: dict
    :return: Shear modulus in the same units as the reference shear modulus.
    :rtype: float
    """

    x = params["V_0"] / volume
    f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0)
    G = pow((1.0 + 2.0 * f), 5.0 / 2.0) * (
        params["G_0"]
        + (3.0 * params["K_0"] * params["Gprime_0"] - 5.0 * params["G_0"]) * f
        + (
            6.0 * params["K_0"] * params["Gprime_0"]
            - 24.0 * params["K_0"]
            - 14.0 * params["G_0"]
            + 9.0 / 2.0 * params["K_0"] * params["Kprime_0"]
        )
        * f
        * f
    )
    return G


class BirchMurnaghanBase(eos.IsothermalEquationOfState):
    """
    Base class for the isothermal Birch Murnaghan equation of state.
    This is third order in strain, and has no temperature dependence.
    However, the shear modulus is sometimes fit to a second order
    function, so if this is the case, you should use that.
    For more see :class:`burnman.birch_murnaghan.BM3Shear2`
    and :class:`burnman.birch_murnaghan.BM3`.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume_third_order(pressure, params)

    def pressure(self, temperature, volume, params):
        return pressure_third_order(params["V_0"] / volume, params)

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]`
        as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus_third_order(volume, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        if self.order == 2:
            return shear_modulus_second_order(volume, params)
        elif self.order == 3:
            return shear_modulus_third_order(volume, params)

    def _molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\\mathcal{E}`
        of the mineral. :math:`[J/mol]`
        """
        x = np.power(volume / params["V_0"], -1.0 / 3.0)
        x2 = x * x
        x4 = x2 * x2
        x6 = x4 * x2

        xi1 = 3.0 * (4.0 - params["Kprime_0"]) / 4.0

        intPdV = (
            -9.0
            / 2.0
            * params["V_0"]
            * params["K_0"]
            * (
                (xi1 + 1.0) * (x4 / 4.0 - x2 / 2.0 + 1.0 / 4.0)
                - xi1 * (x6 / 6.0 - x4 / 4.0 + 1.0 / 12.0)
            )
        )

        return -intPdV + params["E_0"]

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\\mathcal{G}`
        of the mineral. :math:`[J/mol]`
        """
        # G = int VdP = [PV] - int PdV = E + PV

        return (
            self._molar_internal_energy(pressure, temperature, volume, params)
            + volume * pressure
        )

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        if "E_0" not in params:
            params["E_0"] = 0.0
        if "P_0" not in params:
            params["P_0"] = 0.0

        # If G and Gprime are not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if "G_0" not in params:
            params["G_0"] = float("nan")
        if "Gprime_0" not in params:
            params["Gprime_0"] = float("nan")

        # Check that all the required keys are in the dictionary
        expected_keys = ["V_0", "K_0", "Kprime_0", "G_0", "Gprime_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["P_0"] < 0.0:
            warnings.warn("Unusual value for P_0", stacklevel=2)
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-3:
            warnings.warn("Unusual value for V_0", stacklevel=2)
        if params["K_0"] < 1.0e9 or params["K_0"] > 1.0e13:
            warnings.warn("Unusual value for K_0", stacklevel=2)
        if params["Kprime_0"] < 0.0 or params["Kprime_0"] > 20.0:
            warnings.warn("Unusual value for Kprime_0", stacklevel=2)
        if params["G_0"] < 0.0 or params["G_0"] > 1.0e13:
            warnings.warn("Unusual value for G_0", stacklevel=2)
        if params["Gprime_0"] < -5.0 or params["Gprime_0"] > 10.0:
            warnings.warn("Unusual value for Gprime_0", stacklevel=2)


class BM3(BirchMurnaghanBase):
    """
    Third order Birch Murnaghan isothermal equation of state.
    This uses the third order expansion for shear modulus.
    """

    def __init__(self):
        self.order = 3


class BM3Shear2(BirchMurnaghanBase):
    """
    Third order Birch Murnaghan isothermal equation of state.
    This uses the second order expansion for shear modulus.
    """

    def __init__(self):
        self.order = 2


class BM4(BirchMurnaghanBase):
    """
    Base class for the isothermal Birch Murnaghan equation of state.
    This is fourth order in strain, and has no temperature dependence.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume_fourth_order(pressure, params)

    def pressure(self, temperature, volume, params):
        return pressure_fourth_order(params["V_0"] / volume, params)

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus_fourth_order(volume, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return 0.0

    def _molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        x = np.power(volume / params["V_0"], -1.0 / 3.0)
        x2 = x * x
        x4 = x2 * x2
        x6 = x4 * x2
        x8 = x4 * x4

        xi1 = 3.0 * (4.0 - params["Kprime_0"]) / 4.0
        xi2 = (
            3.0
            / 8.0
            * (
                params["K_0"] * params["Kprime_prime_0"]
                + params["Kprime_0"] * (params["Kprime_0"] - 7.0)
            )
            + 143.0 / 24.0
        )

        intPdV = (
            -9.0
            / 2.0
            * params["V_0"]
            * params["K_0"]
            * (
                (xi1 + 1.0) * (x4 / 4.0 - x2 / 2.0 + 1.0 / 4.0)
                - xi1 * (x6 / 6.0 - x4 / 4.0 + 1.0 / 12.0)
                + xi2 * (x8 / 8 - x6 / 2 + 3.0 * x4 / 4.0 - x2 / 2.0 + 1.0 / 8.0)
            )
        )

        return -intPdV + params["E_0"]

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        super().validate_parameters(params)

        if "Kprime_prime_0" not in params:
            raise KeyError("params object missing parameter : Kprime_prime_0")

        if params["Kprime_prime_0"] > 0.0 or params["Kprime_prime_0"] < -10.0:
            warnings.warn("Unusual value for Kprime_prime_0", stacklevel=2)
