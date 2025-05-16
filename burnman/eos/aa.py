# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.


import numpy as np
from scipy.optimize import brentq
import warnings

from . import equation_of_state as eos
from ..constants import gas_constant


class AA(eos.EquationOfState):
    """
    Class for the liquid metal EOS detailed in :cite:`AA1994`.

    This equation of state is complicated because there is not
    a single set of independent variables.

    The equation of state is based on a reference isentrope
    passing through a defined volume and entropy point.
    Internal energy (:math:`E`) at a given volume is calculated
    along this isentrope using a fourth order BM EoS
    (:math:`V_0`, :math:`KS`, :math:`KS'`, :math:`KS''`).

    The temperature along the isentrope is calculated via integration
    of the Grueneisen parameter:

    :math:`\\gamma = \\partial (\\ln T)/\\partial (\\ln \\rho) |_S`

    which gives:

    :math:`T_S/T_0 = \\exp(\\int( \\gamma/\\rho ) d \\rho)`

    Finally, the internal energy away from the reference isentrope
    is calculated as a function of temperature, using an
    expression for the isochoric heat capacity as a function of
    volume and temperature.

    In this BurnMan implementation, the Helmholtz energy is used
    as the natural potential, with volume and temperature as the
    natural variables.

    We note that :cite:`AA1994` also include a detailed description
    of the Gruneisen parameter as a function of volume and energy
    (Equation 15), and use this to determine the temperature along
    the principal isentrope (Equations B1-B10) and the
    thermal pressure away from that isentrope (Equation 23).
    However, the thermal pressure expression is inconsistent with
    the equation of state away from the principal isentrope.

    Note: the expression for :math:`\\Lambda` (Appendix C) does not
    reproduce Figure 5. We assume that the figure is correct, and
    that the correct expression has the form:
    :math:`F(-325.23 + 302.07 (\\rho/\\rho_0) + 30.45 (\\rho/\\rho_0)^{0.4})`.
    """

    def _ABTheta(self, volume, params):
        """
        Electronic heat capacity functions
        """
        Vfrac = volume / params["V_0"]

        A = params["a"][0] + params["a"][1] * Vfrac  # A2
        B = params["b"][0] + params["b"][1] * Vfrac * Vfrac  # A3
        Theta = params["Theta"][0] * np.power(Vfrac, -params["Theta"][1])  # A4

        return A, B, Theta

    def _lambdaxi(self, volume, params):
        """
        Potential heat capacity functions
        """
        rhofrac = params["V_0"] / volume
        xi = params["xi_0"] * np.power(rhofrac, -0.6)  # A16
        F = 1.0 / (1.0 + np.exp((rhofrac - params["F"][0]) / params["F"][1]))  # A18
        # lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac) + params['lmda'][2])*np.power(rhofrac, 0.4) # A17

        # The following provides a good fit to Figure 5
        lmda = F * (
            params["lmda"][0]
            + params["lmda"][1] * rhofrac
            + params["lmda"][2] * np.power(rhofrac, 0.4)
        )

        return lmda, xi

    def _rhofracxksis(self, volume, params):
        """
        Functions for the fourth order Birch-Murnaghan equation of state
        """
        rhofrac = params["V_0"] / volume
        x = np.power(rhofrac, 1.0 / 3.0)  # Equation 18
        ksi1 = 0.75 * (4.0 - params["Kprime_S"])  # Equation 19
        ksi2 = (
            0.375
            * (
                params["K_S"] * params["Kprime_prime_S"]
                + params["Kprime_S"] * (params["Kprime_S"] - 7.0)
            )
            + 143.0 / 24.0
        )  # Equation 20
        return rhofrac, x, ksi1, ksi2

    def _reference_temperature(self, volume, params):
        """
        Temperature along the reference isentrope
        """

        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(volume, params)

        # Equation B6 -- B10
        a1 = ksi2 / 8.0
        a2 = (ksi1 + 3.0 * ksi2) / 6.0
        a3 = (1.0 + 2.0 * ksi1 + 3.0 * ksi2) / 4.0
        a4 = (1.0 + ksi1 + ksi2) / 2.0
        a5 = (6.0 + 4.0 * ksi1 + 3.0 * ksi2) / 24.0
        n = params["grueneisen_n"]

        # Equation B5
        Ts = params["T_0"] * np.exp(
            params["grueneisen_0"] * np.log(rhofrac)
            + 13.5
            * params["grueneisen_prime"]
            * params["V_0"]
            * params["K_S"]
            * (
                (a1 / (3 * n + 8.0)) * (np.power(x, (3 * n + 8.0)) - 1.0)
                - (a2 / (3 * n + 6.0)) * (np.power(x, (3 * n + 6.0)) - 1.0)
                + (a3 / (3 * n + 4.0)) * (np.power(x, (3 * n + 4.0)) - 1.0)
                - (a4 / (3 * n + 2.0)) * (np.power(x, (3 * n + 2.0)) - 1.0)
                + (a5 / (3 * n + 0.0)) * (np.power(x, (3 * n + 0.0)) - 1.0)
            )
        )

        return Ts

    '''
    def _reference_pressure(self, volume, params):
        """
        Pressure along the reference isentrope (Eq. 17).
        Currently unused as the other pressure terms have not yet been derived.
        """
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(volume, params)
        x2 = x * x
        x3 = x * x * x
        x5 = x3 * x2
        x7 = x5 * x2

        Ps = (
            1.5
            * params["K_S"]
            * (x7 - x5)
            * (1.0 + ksi1 - ksi1 * x2 + ksi2 * (x2 - 1.0) * (x2 - 1.0))
        )  # Eq. 17

        return Ps
    '''

    def _isentropic_energy_change(self, volume, params):
        """
        Birch Murnaghan equation of state expression for the energy change along an isentrope
        """
        _, x, ksi1, ksi2 = self._rhofracxksis(volume, params)
        x2 = x * x
        x4 = x2 * x2
        x6 = x4 * x2
        x8 = x4 * x4

        E_S = (
            4.5
            * params["V_0"]
            * params["K_S"]
            * (
                (ksi1 + 1.0) * (x4 / 4.0 - x2 / 2.0 + 1.0 / 4.0)
                - ksi1 * (x6 / 6.0 - x4 / 4.0 + 1.0 / 12.0)
                + ksi2 * (x8 / 8.0 - x6 / 2.0 + 3.0 * x4 / 4.0 - x2 / 2.0 + 1.0 / 8.0)
            )
        )  # Eq. 21
        return E_S

    def _molar_heat_capacity_v(self, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """

        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)

        # HT limit of kinetic contribution (just after Equation 29.)
        C_kin = 1.5 * params["n"] * gas_constant
        C_e = A * (
            1.0 - (Theta * Theta) / (Theta * Theta + temperature * temperature)
        ) + B * np.power(
            temperature, 0.6
        )  # Equation A1
        C_pot = (lmda * temperature + xi * params["theta"]) / (
            params["theta"] + temperature
        )  # Equation A15

        return C_kin + C_e + C_pot

    def entropy(self, pressure, temperature, volume, params):
        Trel = self._reference_temperature(volume, params)
        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)

        C_kin = 1.5 * params["n"] * gas_constant
        S_0 = params["S_0"]
        S_k = C_kin * np.log(temperature / Trel)

        S_pot = (lmda - xi) * (
            np.log(params["theta"] + temperature) - np.log(params["theta"] + Trel)
        ) + xi * np.log(temperature / Trel)

        sqTh = Theta * Theta
        S_e = 0.5 * A * (
            np.log(sqTh + temperature * temperature) - np.log(sqTh + Trel * Trel)
        ) + 5.0 / 3.0 * B * (np.power(temperature, 0.6) - np.power(Trel, 0.6))

        return S_0 + S_k + S_pot + S_e

    def _helmholtz_thermal(self, temperature, volume, params):
        Trel = self._reference_temperature(volume, params)
        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)

        C_kin = 1.5 * params["n"] * gas_constant
        F_0 = -params["S_0"] * (temperature - Trel)
        F_k = -C_kin * (temperature * (np.log(temperature / Trel) - 1.0) + Trel)

        F_T = (lmda - xi) * (temperature + params["theta"]) * np.log(
            temperature + params["theta"]
        ) + temperature * (
            -lmda
            + xi * np.log(temperature / Trel)
            - (lmda - xi) * np.log(params["theta"] + Trel)
        )
        F_Trel = (lmda - xi) * params["theta"] * np.log(
            Trel + params["theta"]
        ) - Trel * lmda
        F_pot = F_Trel - F_T

        sqTh = Theta * Theta

        F_e = (
            A * (temperature - Trel)
            - 25.0 / 24.0 * B * (np.power(temperature, 1.6) - np.power(Trel, 1.6))
            + 5.0 / 3.0 * B * (temperature - Trel) * np.power(Trel, 0.6)
            - A * Theta * (np.arctan(temperature / Theta) - np.arctan(Trel / Theta))
            - 0.5
            * A
            * temperature
            * (np.log(sqTh + temperature * temperature) - np.log(sqTh + Trel * Trel))
        )

        return F_0 + F_k + F_pot + F_e

    def _helmholtz_energy(self, temperature, volume, params):
        """
        Returns the Helmholtz energy at the temperature and volume of the mineral [J/mol]
        F = E - TS
        """
        Eref = params["E_0"] + self._isentropic_energy_change(volume, params)
        Tref = self._reference_temperature(volume, params)
        Fref = Eref - Tref * self.entropy(0.0, Tref, volume, params)
        Fth = self._helmholtz_thermal(temperature, volume, params)
        return Fref + Fth

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume [Pa]
        """
        dV = volume * 1.0e-4
        V1 = volume - dV / 2.0
        V2 = volume + dV / 2.0
        F1 = self._helmholtz_energy(temperature, V1, params)
        F2 = self._helmholtz_energy(temperature, V2, params)
        return -(F2 - F1) / dV

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """

        def delta_pressure(volume):
            return pressure - self.pressure(temperature, volume, params)

        V_0 = params["V_0"]
        return brentq(delta_pressure, 0.1 * V_0, 2.0 * V_0)

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        # K_T = -V * dP/dV
        dV = volume * 1.0e-4
        P0 = self.pressure(temperature, volume - 0.5 * dV, params)
        P1 = self.pressure(temperature, volume + 0.5 * dV, params)

        K_T = -volume * (P1 - P0) / dV
        return K_T

    def _aKT(self, temperature, volume, params):
        """
        Returns the product of the thermal expansivity and isothermal bulk modulus. :math:`[Pa/K]`
        """

        # aKT = dP/dT at constant volume
        dT = 0.1
        P0 = self.pressure(temperature - 0.5 * dT, volume, params)
        P1 = self.pressure(temperature + 0.5 * dT, volume, params)

        return (P1 - P0) / dT

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns the thermal expansivity. :math:`[1/K]`
        """

        aK_T = self._aKT(temperature, volume, params)
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)

        return aK_T / K_T

    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns the Grueneisen parameter. [unitless]
        """
        aK_T = self._aKT(temperature, volume, params)
        C_V = self._molar_heat_capacity_v(temperature, volume, params)

        return aK_T * volume / C_V

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """

        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self._molar_heat_capacity_v(temperature, volume, params)
        C_p = C_v * (1.0 + gr * alpha * temperature)

        return C_p

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        F + PV
        """
        F = self._helmholtz_energy(temperature, volume, params)
        return F + pressure * self.volume(pressure, temperature, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus. :math:`[Pa]`
        Zero for a liquid
        """
        return 0.0

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["P_0", "T_0", "S_0", "molar_mass", "grueneisen_0"]

        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["T_0"] < 0.0:
            warnings.warn("Unusual value for T_0", stacklevel=2)
        if params["molar_mass"] < 0.001 or params["molar_mass"] > 10.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)
        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
