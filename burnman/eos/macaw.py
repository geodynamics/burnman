from __future__ import absolute_import

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences.
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
import warnings
import numpy as np


# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:

    def jit(nopython=True):
        def decorator(fn):
            return fn

        return decorator


@jit(nopython=True)
def make_params(K0, K0_prime, K_infinity_prime):
    a = (
        16.0 * np.power(K0_prime, 3.0)
        + 84.0 * np.power(K0_prime, 2.0)
        + 192.0 * K0_prime
        - 972.0 * K_infinity_prime
        + 1177.0
    )
    b = 2.0 * np.power(K0_prime, 2.0) + 7.0 * K0_prime - 27.0 * K_infinity_prime + 38.0
    omega = np.power((a + np.sqrt(a * a - 32.0 * b * b * b)), 1.0 / 3.0)
    C = (
        (11.0 / 6.0)
        + (1.0 / 3.0) * K0_prime
        - K_infinity_prime
        + (np.power(2, -1.0 / 3.0) / 6) * omega
        + (np.power(2, 1.0 / 3.0) / 3) * (b / omega)
    )
    B = K_infinity_prime - 1.0
    A = K0 / (B - 0.5 * C + np.power(B + C, 2.0))
    return A, B, C


class MACAW(eos.EquationOfState):
    """
    Class for the MACAW equation of state
    detailed in Lozano and Aslam (2022; https://doi.org/10.1063/5.0076897).

    This equation of state has no temperature dependence.
    """

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        A, B, C = make_params(params["K_0"], params["Kprime_0"], params["Kprime_inf"])
        Vrel = volume / params["V_0"]
        term1 = A * np.power(Vrel, -(B + 1))
        term2 = np.exp((2.0 / 3.0) * C * (1 - np.power(Vrel, 1.5)))
        term3 = np.power(C * np.power(Vrel, 1.5) + B, 2.0) - (
            0.5 * C * np.power(Vrel, 1.5) - B
        )
        return term1 * term2 * term3

    def volume(self, pressure, temperature, params):
        """
        Get the Vinet volume at a reference temperature for a given
        pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
        """

        def delta_pressure(x):
            return self.pressure(0.0, x, params) - pressure

        V = opt.brentq(delta_pressure, 0.1 * params["V_0"], 1.5 * params["V_0"])
        return V

    def pressure(self, temperature, volume, params):
        """
        Returns pressure :math:`[Pa]` as a function of volume :math:`[m^3]`.
        """
        A, B, C = make_params(params["K_0"], params["Kprime_0"], params["Kprime_inf"])
        Vrel = volume / params["V_0"]
        term1 = A * np.power(Vrel, -(B + 1.0))
        term2 = np.exp((2.0 / 3.0) * C * (1.0 - np.power(Vrel, 1.5)))
        term3 = C * np.power(Vrel, 1.5) + B
        return term1 * term2 * term3 - A * (B + C) + params["P_0"]

    def _molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        A, B, C = make_params(params["K_0"], params["Kprime_0"], params["Kprime_inf"])
        Vrel = volume / params["V_0"]
        I1 = -params["V_0"] * (
            np.power(Vrel, -B) * np.exp((2.0 / 3.0) * C * (1.0 - np.power(Vrel, 1.5)))
            - 1.0
        )
        I0 = (-A * (B + C) + params["P_0"]) * params["V_0"] * (Vrel - 1.0)
        return -A * I1 - I0

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        return (
            self._molar_internal_energy(pressure, temperature, volume, params)
            + pressure * volume
        )

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return 1.0e99

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        """
        return 0.0

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, return a very small number. :math:`[J/K/mol]`
        """
        return 1.0e-99

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, return zero. :math:`[1/K]`
        """
        return 0.0

    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, return zero. :math:`[unitless]`
        """
        return 0.0

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters.
        The value for :math:`K'_{\\infty}` is thermodynamically bounded
        between 5/3 and :math:`K'_0` :cite:`StaceyDavis2004`.
        """

        if "E_0" not in params:
            params["E_0"] = 0.0
        if "P_0" not in params:
            params["P_0"] = 1.0e5

        # Check that all the required keys are in the dictionary
        expected_keys = ["V_0", "K_0", "Kprime_0", "Kprime_inf"]
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
        if params["Kprime_0"] < 0.0 or params["Kprime_0"] > 10.0:
            warnings.warn("Unusual value for Kprime_0", stacklevel=2)
        if params["Kprime_inf"] < 1 + 45.0 / 29.0:
            warnings.warn("Value for Kprime_inf below recommended value", stacklevel=2)
        if params["Kprime_inf"] > params["Kprime_0"]:
            warnings.warn("Kprime_inf should be less than Kprime_0", stacklevel=2)
