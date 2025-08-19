# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


import warnings
import numpy as np

from . import equation_of_state as eos


def tait_constants(params):
    """
    returns parameters for the modified Tait equation of state
    derived from K_T and its two first pressure derivatives
    EQ 4 from Holland and Powell, 2011
    """
    a = (1.0 + params["Kprime_0"]) / (
        1.0 + params["Kprime_0"] + params["K_0"] * params["Kdprime_0"]
    )
    b = params["Kprime_0"] / params["K_0"] - params["Kdprime_0"] / (
        1.0 + params["Kprime_0"]
    )
    c = (1.0 + params["Kprime_0"] + params["K_0"] * params["Kdprime_0"]) / (
        params["Kprime_0"] * params["Kprime_0"]
        + params["Kprime_0"]
        - params["K_0"] * params["Kdprime_0"]
    )
    return a, b, c


def pressure_modified_tait(Vrel, params):
    """
    Pressure according to the modified Tait equation of state.
    Equation 2 in :cite:`HP2011`.

    :param Vrel: Volume divided by the reference volume.
    :type Vrel: float or numpy array
    :param params: Parameter dictionary
    :type params: dictionary
    :return: pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    :rtype: float or numpy array
    """
    a, b, c = tait_constants(params)
    return (np.power((Vrel + a - 1.0) / a, -1.0 / c) - 1.0) / b + params["P_0"]


def volume(pressure, params):
    """
    Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
    Equation 12 in :cite:`HP2011`.
    """
    a, b, c = tait_constants(params)
    Vrel = 1.0 - a * (1.0 - np.power((1.0 + b * (pressure - params["P_0"])), -1.0 * c))
    return Vrel * params["V_0"]


def bulk_modulus(pressure, params):
    """
    Returns isothermal bulk modulus :math:`K_T` of the mineral. :math:`[Pa]`.
    EQ 13+2
    """
    a, b, c = tait_constants(params)
    return (
        params["K_0"]
        * (1.0 + b * (pressure - params["P_0"]))
        * (a + (1.0 - a) * np.power((1.0 + b * (pressure - params["P_0"])), c))
    )


def intVdP(pressure, params):
    """
    Returns the integral of VdP for the mineral. :math:`[J]`.
    EQ 13
    """
    a, b, c = tait_constants(params)
    psubpth = pressure - params["P_0"]

    if pressure != params["P_0"]:
        intVdP = (
            (pressure - params["P_0"])
            * params["V_0"]
            * (
                1.0
                - a
                + (
                    a
                    * (1.0 - np.power((1.0 + b * (psubpth)), 1.0 - c))
                    / (b * (c - 1.0) * (pressure - params["P_0"]))
                )
            )
        )
    else:
        intVdP = 0.0
    return intVdP


class MT(eos.IsothermalEquationOfState):
    """
    Base class for the generic modified Tait equation of state.
    References for this can be found in :cite:`HC1974`
    and :cite:`HP2011` (followed here).

    An instance "m" of a Mineral can be assigned this
    equation of state with the command m.set_method('mt')
    (or by initialising the class with the param
    equation_of_state = 'mt').
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume(pressure, params)

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        return pressure_modified_tait(volume / params["V_0"], params)

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(pressure, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented in the Modified Tait EoS. :math:`[Pa]`
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.0

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        # G = int VdP = [PV] - int PdV = E + PV
        a, b, c = tait_constants(params)

        intVdP = params["V_0"] * (
            a
            / (b * (1.0 - c))
            * (np.power(b * (pressure - params["P_0"]) + 1.0, 1.0 - c) - 1.0)
            + (1.0 - a) * (pressure - params["P_0"])
        )

        return intVdP + params["F_0"] + params["V_0"] * params["P_0"]

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
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

        # G and Gprime are not defined in this equation of state,
        # We can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if "G_0" not in params:
            params["G_0"] = float("nan")
        if "Gprime_0" not in params:
            params["Gprime_0"] = float("nan")

        # Check that all the required keys are in the dictionary
        expected_keys = ["V_0", "K_0", "Kprime_0", "Kdprime_0", "G_0", "Gprime_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["P_0"] < 0.0:
            warnings.warn("Unusual value for P_0", stacklevel=2)
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-2:
            warnings.warn("Unusual value for V_0", stacklevel=2)
        if params["K_0"] < 1.0e9 or params["K_0"] > 1.0e13:
            warnings.warn("Unusual value for K_0", stacklevel=2)
        if params["Kprime_0"] < 0.0 or params["Kprime_0"] > 40.0:
            warnings.warn("Unusual value for Kprime_0", stacklevel=2)
        if params["G_0"] < 0.0 or params["G_0"] > 1.0e13:
            warnings.warn("Unusual value for G_0", stacklevel=2)
        if params["Gprime_0"] < -5.0 or params["Gprime_0"] > 10.0:
            warnings.warn("Unusual value for Gprime_0", stacklevel=2)
