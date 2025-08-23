# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from abc import ABC, abstractmethod


def singleton(cls):
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]

    return get_instance


class DebyeTemperatureModelBase(ABC):
    """Abstract base class for Debye temperature models."""

    @abstractmethod
    def value(self, Vrel: float, params: dict) -> float:
        """Return the Debye temperature for the given relative volume."""

    @abstractmethod
    def _gamma(self, Vrel: float, params: dict) -> float:
        """Return the Grueneisen parameter at given relative volume."""

    @abstractmethod
    def _gamma_prime(self, Vrel: float, params: dict) -> float:
        """Return the derivative of Grueneisen parameter wrt relative volume."""

    def dVrel(self, Vrel: float, params: dict) -> float:
        """First derivative of Debye temperature wrt relative volume."""
        theta = self.value(Vrel, params)
        gamma = self._gamma(Vrel, params)
        return -theta * gamma / Vrel

    def dVrel2(self, Vrel: float, params: dict) -> float:
        """Second derivative of Debye temperature wrt relative volume."""
        theta = self.value(Vrel, params)
        gamma = self._gamma(Vrel, params)
        gamma_prime = self._gamma_prime(Vrel, params)
        return theta / Vrel**2 * (gamma * (gamma + 1.0) - Vrel * gamma_prime)


@singleton
class SLB(DebyeTemperatureModelBase):
    """
    A class that provides the Debye temperature as a function of relative volume,
    and its first and second derivatives with respect to relative volume.
    This class is used in the Extended Mie-Grueneisen-Debye equation of state.

    This class is based on a finite strain expansion for the
    squared frequencies of the phonon modes (Stixrude and Lithgow-Bertelloni, 2005).
    The Debye temperature is:
    :math:`\\Theta(V) = \\Theta_0 \\sqrt{1 + a_1 f + 0.5 a_2 f^2}`, where
    :math:`f = 0.5 (f)^{-2/3} - 1`, and :math:`a_1` and :math:`a_2` are
    parameters that define the Grueneisen parameter's dependence on volume and equal
    :math:`a_1 = 6 \\gamma_0`, :math:`a_2 = -12 \\gamma_0 + 36 \\gamma_0^2 - 18 q_0 \\gamma_0`.

    :param Vrel: Relative volume, defined as :math:`V/V_0`, where V_0 is the reference volume.
    :type Vrel: float
    :param params: A dictionary containing the parameters for the model.
    :type params: dict
    """

    def value(self, Vrel, params):
        # This method computes the Debye temperature.
        a1_ii, a2_iikk, f = self.a1_a2_and_f(Vrel, params)
        nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * f * f
        debye_temperature = params["Debye_0"] * np.sqrt(nu_o_nu0_sq)  # EQ 41
        return debye_temperature

    def _gamma(self, Vrel, params):
        a1_ii, a2_iikk, f = self.a1_a2_and_f(Vrel, params)
        nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * f * f
        gr = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)
        return gr

    def _gamma_prime(self, Vrel, params):
        a1_ii, a2_iikk, f = self.a1_a2_and_f(Vrel, params)
        fprime = -(1.0 / 3.0) * Vrel ** (-5.0 / 3.0)

        N = (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)
        D = 6.0 * (1.0 + a1_ii * f + 0.5 * a2_iikk * f * f)

        Nprime = 2.0 * (a1_ii + a2_iikk * f) + (2.0 * f + 1.0) * a2_iikk
        Dprime = 6.0 * (a1_ii + a2_iikk * f)
        gprime = ((Nprime * D - N * Dprime) / (D * D)) * fprime

        return gprime

    def a1_a2_and_f(self, Vrel, params):
        """
        Returns a1, a2 and f for the given relative volume.
        This is useful for debugging and testing purposes.
        """
        a1_ii = 6.0 * params["grueneisen_0"]
        a2_iikk = (
            -12.0 * params["grueneisen_0"]
            + 36.0 * pow(params["grueneisen_0"], 2.0)
            - 18.0 * params["q_0"] * params["grueneisen_0"]
        )
        f = 0.5 * (pow(Vrel, -2.0 / 3.0) - 1.0)
        return a1_ii, a2_iikk, f

    def validate_parameters(self, params):
        # Check for all required keys
        expected_keys = ["grueneisen_0", "q_0", "Debye_0"]
        for key in expected_keys:
            if key not in params:
                raise AttributeError(f"params dictionary must contain a '{key}' key")


@singleton
class PowerLawGammaSimple(DebyeTemperatureModelBase):
    """
    A class that provides the Debye temperature as a function of relative volume,
    and its first and second derivatives with respect to relative volume.
    This class is used in the Extended Mie-Grueneisen-Debye equation of state.

    This class is based on a power-law model for the Grueneisen parameter,
    which is defined as (Matas et al., 2007):
    :math:`\\gamma(V) = \\gamma_0 (V/V_0)^{q_0}`
    where :math:`\\gamma_0` and :math:`q_0`
    are parameters that define the Grueneisen parameter's dependence on volume.

    The Debye temperature is computed as:
    :math:`\\Theta(V) = \\Theta_0 \\exp(-\\int_1^{V/V_0} \\gamma(x)/x dx)`.

    :param Vrel: Relative volume, defined as :math:`V/V_0`, where V_0 is the reference volume.
    :type Vrel: float
    :param params: A dictionary containing the parameters for the model.
    :type params: dict
    """

    def value(self, Vrel, params):
        # This method computes the Debye temperature.
        return params["Debye_0"] * np.exp(
            -params["grueneisen_0"]
            / params["q_0"]
            * (np.power(Vrel, params["q_0"]) - 1.0)
        )

    def _gamma(self, Vrel, params):
        return params["grueneisen_0"] * np.power(Vrel, params["q_0"])

    def _gamma_prime(self, Vrel, params):
        return (
            params["grueneisen_0"] * params["q_0"] * np.power(Vrel, params["q_0"] - 1.0)
        )

    def validate_parameters(self, params):
        # Check for all required keys
        expected_keys = ["grueneisen_0", "q_0", "Debye_0"]
        for key in expected_keys:
            if key not in params:
                raise AttributeError(f"params dictionary must contain a '{key}' key")


@singleton
class PowerLawGamma(DebyeTemperatureModelBase):
    """
    A class that provides the Debye temperature as a function of relative volume,
    and its first and second derivatives with respect to relative volume.
    This class is used in the Extended Mie-Grueneisen-Debye equation of state.

    This class is based on a power-law model for the Grueneisen parameter,
    which is defined as:
    :math:`\\gamma(V) = grueneisen_0 + c_1 ((V/V_0)^{q_1} - 1) + c_2 ((V/V_0)^{q_2} - 1)`
    where :math:`grueneisen_0`, :math:`c_1`, :math:`c_2`, :math:`q_1`, and :math:`q_2`
    are parameters that define the Grueneisen parameter's dependence on volume.

    The Debye temperature is computed as:
    :math:`\\Theta(V) = \\Theta_0 \\exp(-\\int_1^{V/V_0} \\gamma(x)/x dx)`.

    :param Vrel: Relative volume, defined as :math:`V/V_0`, where V_0 is the reference volume.
    :type Vrel: float
    :param params: A dictionary containing the parameters for the model.
    :type params: dict
    """

    def value(self, Vrel, params):
        # This method computes the Debye temperature.
        return (
            params["Debye_0"]
            * np.power(Vrel, params["c_2"] + params["c_1"] - params["grueneisen_0"])
            * np.exp(
                -params["c_1"] / params["q_1"] * (np.power(Vrel, params["q_1"]) - 1.0)
                - params["c_2"] / params["q_2"] * (np.power(Vrel, params["q_2"]) - 1.0)
            )
        )

    def _gamma(self, Vrel, params):
        return (
            params["grueneisen_0"]
            + params["c_1"] * (np.power(Vrel, params["q_1"]) - 1.0)
            + params["c_2"] * (np.power(Vrel, params["q_2"]) - 1.0)
        )

    def _gamma_prime(self, Vrel, params):
        return params["c_1"] * params["q_1"] * np.power(
            Vrel, params["q_1"] - 1.0
        ) + params["c_2"] * params["q_2"] * np.power(Vrel, params["q_2"] - 1.0)

    def validate_parameters(self, params):
        # Check for all required keys
        expected_keys = ["Debye_0", "grueneisen_0", "c_1", "c_2", "q_1", "q_2"]
        for key in expected_keys:
            if key not in params:
                raise AttributeError(f"params dictionary must contain a '{key}' key")
