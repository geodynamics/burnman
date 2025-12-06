# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences.
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
from ..utils.math import generalised_gammainc
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
def _make_params(K_0, Kp_0, Kp_inf, Kdp_0):
    dKpdlnV_zero = -Kdp_0 * K_0
    c = Kp_inf
    a = dKpdlnV_zero / (Kp_0 - Kp_inf)
    b = (Kp_0 - Kp_inf) / a
    return a, b, c


class SPOCK(eos.IsothermalEquationOfState):
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
        ai, bi, ci = _make_params(
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
        ai, bi, ci = _make_params(
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

    def _molar_helmholtz_energy(self, pressure, temperature, volume, params):
        """
        Internal function returning the Helmholtz energy
        :math:`\\mathcal{F}` of the mineral. :math:`[J/mol]`
        """
        ai, bi, ci = _make_params(
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

        return params["F_0"] + params["P_0"] * (volume - params["V_0"]) + f * (I1 - I2)

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        return (
            self._molar_helmholtz_energy(pressure, temperature, volume, params)
            + pressure * volume
        )

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        This equation of state is athermal, so the function returns a very large number.
        """
        return 1.0e99

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters.
        The value for :math:`K'_{\\infty}` is thermodynamically bounded
        between 5/3 and :math:`K'_0` :cite:`StaceyDavis2004`.
        """

        if "F_0" not in params:
            params["F_0"] = 0.0
        if "P_0" not in params:
            params["P_0"] = 1.0e5

        if "E_0" in params:
            raise KeyError(
                "Isothermal equations of state should be "
                "defined in terms of Helmholtz free energy "
                "F_0, not internal energy E_0."
            )

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
