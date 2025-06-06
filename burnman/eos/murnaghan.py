# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from . import equation_of_state as eos
import warnings
import numpy as np


def volume(pressure, V_0, K_0, Kprime_0):
    return V_0 * np.power(1.0 + (pressure * Kprime_0 / K_0), -1.0 / Kprime_0)


def pressure(volume, V_0, K_0, Kprime_0):
    return K_0 / Kprime_0 * (np.power(volume / V_0, -Kprime_0) - 1.0)


def bulk_modulus(pressure, K_0, Kprime_0):
    return K_0 + pressure * Kprime_0


def energy(volume, E_0, V_0, K_0, Kprime_0):
    Vrel = volume / V_0
    return E_0 + K_0 * V_0 * (
        np.power(Vrel, 1.0 - Kprime_0) / (Kprime_0 * (Kprime_0 - 1))
        + Vrel / Kprime_0
        - 1.0 / (Kprime_0 - 1.0)
    )


def intVdP(pressure, V_0, K_0, Kprime_0):
    return (
        V_0
        * K_0
        * ((np.power(1.0 + (pressure * Kprime_0 / K_0), 1.0 - (1.0 / Kprime_0))) - 1.0)
        / (Kprime_0 - 1.0)
    )


class Murnaghan(eos.IsothermalEquationOfState):
    """
    Base class for the isothermal Murnaghan equation of state,
    as described in :cite:`Murnaghan1944`.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume(pressure, params["V_0"], params["K_0"], params["Kprime_0"])

    def pressure(self, temperature, volume, params):
        return pressure(volume, params["V_0"], params["K_0"], params["Kprime_0"])

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]`
        as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(pressure, params["K_0"], params["Kprime_0"])

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        Currently not included in the Murnghan EOS, so omitted.
        """
        return 0.0

    def _molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\\mathcal{E}` of the mineral.
        :math:`[J/mol]`
        """
        return energy(
            volume, params["E_0"], params["V_0"], params["K_0"], params["Kprime_0"]
        )

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\\mathcal{G}` of the mineral.
        :math:`[J/mol]`
        """
        # G = E + PV
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

        # G is not included in the Murnaghan EOS so we shall set them to NaN's
        if "G_0" not in params:
            params["G_0"] = float("nan")
        if "Gprime_0" not in params:
            params["Gprime_0"] = float("nan")

        # check that all the required keys are in the dictionary
        expected_keys = ["V_0", "K_0", "Kprime_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # now check that the values are reasonable.  I mostly just
        # made up these values from experience, and we are only
        # raising a warning.  Better way to do this? [IR]
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-3:
            warnings.warn("Unusual value for V_0", stacklevel=2)
        if params["K_0"] < 1.0e9 or params["K_0"] > 1.0e13:
            warnings.warn("Unusual value for K_0", stacklevel=2)
        if params["Kprime_0"] < -5.0 or params["Kprime_0"] > 30.0:
            warnings.warn("Unusual value for Kprime_0", stacklevel=2)
