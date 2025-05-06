from __future__ import absolute_import

# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences.
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from scipy.special import gamma, gammainc, exp1
from . import equation_of_state as eos
import warnings
import numpy as np


# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:

    def jit(fn):
        return fn


def generalised_gammainc(a, x1, x2):
    """
    An implementation of the generalised incomplete gamma
    function. Computed using the relationship with the regularised
    lower incomplete gamma function (scipy.special.gammainc).
    Uses the recurrence relation wherever a<0.

    We could have used mpmath.gammainc(a, x1, x2) directly,
    but it is significantly slower than this implementation.
    """
    n = int(-np.floor(a))
    if n > 0:
        a = a + n
        u_gamma = (
            exp1(x1) - exp1(x2)
            if np.abs(a) < 1.0e-12
            else (gammainc(a, x2) - gammainc(a, x1)) * gamma(a)
        )

        esubx1 = np.exp(-x1)
        esubx2 = np.exp(-x2)
        for _ in range(n):
            a = a - 1.0
            u_gamma = (
                u_gamma + np.power(x2, a) * esubx2 - np.power(x1, a) * esubx1
            ) / a
        return u_gamma
    else:
        return (gammainc(a, x2) - gammainc(a, x1)) * gamma(a)


@jit(nopython=True)
def make_params(K_0, Kp_0, Kp_inf, Kdp_0):
    dKpdlnV_zero = -Kdp_0 * K_0
    c = Kp_inf
    a = dKpdlnV_zero / (Kp_0 - Kp_inf)
    b = (Kp_0 - Kp_inf) / a
    return a, b, c


class SPOCK(eos.EquationOfState):
    """
    Class for the Scaled Power Of Compression K-prime equation of state.
    This equation is derived from the assumption that K' = b*(V/V_0)^a.

    This equation of state has no temperature dependence.
    """

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        ai, bi, ci = make_params(
            params["K_0"], params["Kprime_0"], params["Kprime_inf"], params["Kdprime_0"]
        )

        lnVrel = np.log(volume / params["V_0"])
        return params["K_0"] * np.exp(-bi * (np.exp(ai * lnVrel) - 1.0) - ci * lnVrel)

    def volume(self, pressure, temperature, params):
        """
        Returns volume at a given pressure :math:`[Pa]` in :math:`[m^3]`
        """

        def delta_pressure(x):
            return self.pressure(0.0, x, params) - pressure

        V = opt.brentq(delta_pressure, 0.1 * params["V_0"], 1.5 * params["V_0"])
        return V

    def pressure(self, temperature, volume, params):
        """
        Returns pressure :math:`[Pa]` as a function of volume :math:`[m^3]`.
        """
        ai, bi, ci = make_params(
            params["K_0"], params["Kprime_0"], params["Kprime_inf"], params["Kdprime_0"]
        )
        lnVrel = np.log(volume / params["V_0"])
        return params["P_0"] + (
            params["K_0"]
            * np.exp(bi)
            / ai
            * np.power(bi, ci / ai)
            * (generalised_gammainc(-ci / ai, bi * np.exp(ai * lnVrel), bi))
        )

    def _molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Internal function returning the internal energy
        :math:`\\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        ai, bi, ci = make_params(
            params["K_0"], params["Kprime_0"], params["Kprime_inf"], params["Kdprime_0"]
        )
        lnVrel = np.log(volume / params["V_0"])
        f = (
            -params["V_0"]
            * params["K_0"]
            * np.exp(bi)
            / ai
            * np.power(bi, (ci - 1.0) / ai)
        )

        Vrel = np.exp(lnVrel)
        I1 = (
            np.power(bi, 1.0 / ai)
            * Vrel
            * (generalised_gammainc(-ci / ai, bi * np.exp(ai * lnVrel), bi))
        )
        I2 = generalised_gammainc((1.0 - ci) / ai, bi * np.exp(ai * lnVrel), bi)

        return params["E_0"] + params["P_0"] * (volume - params["V_0"]) + f * (I1 - I2)

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
        This equation of state is athermal, so the function returns a very large number.
        """
        return 1.0e99

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        This equation of state is athermal, so the function returns zero.
        """
        return 0.0

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the molar heat capacity at constant pressure. :math:`[J/K/mol]`
        This equation of state is athermal, so the function returns zero.
        """
        return 1.0e-99

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns the thermal expansivity at constant pressure. :math:`[1/K]`
        This equation of state is athermal, so the function returns zero.
        """
        return 0.0

    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns the grueneisen parameter. Unitless
        This equation of state is athermal, so the function returns zero.
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
        expected_keys = ["V_0", "K_0", "Kprime_0", "Kdprime_0", "Kprime_inf"]
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
        if params["Kdprime_0"] > 0.0:
            warnings.warn("Kdprime_0 should be negative", stacklevel=2)
        if (
            -params["Kdprime_0"] * params["K_0"]
            < params["Kprime_0"] - params["Kprime_inf"]
        ):
            warnings.warn(
                "Kdprime_0*K_0 is expected to be more "
                "negative than (Kprime_0 - Kprime_inf)",
                stacklevel=2,
            )
        if (
            params["Kprime_inf"] < 5.0 / 3.0
            or params["Kprime_inf"] > params["Kprime_0"]
        ):
            warnings.warn(
                "Kprime_inf is expected to be greater "
                "than the Thomas-Fermi limit (5/3)",
                stacklevel=2,
            )
