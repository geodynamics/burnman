# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import warnings
import pkgutil
from scipy.optimize import brentq
from scipy.special import binom

from . import equation_of_state as eos
from ..constants import gas_constant
from ..utils.math import bracket


class BroshCalphad(eos.EquationOfState):
    """
    Class for the high pressure CALPHAD equation of state by
    :cite:`Brosh2007`.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        X = [
            1.0
            / (
                1.0
                - params["a"][i - 2]
                + params["a"][i - 2]
                * np.power(
                    1.0 + i / (3.0 * params["a"][i - 2]) * pressure / params["K_0"],
                    1.0 / float(i),
                )
            )
            for i in range(2, 6)
        ]
        V_c = params["V_0"] * np.sum(
            [params["c"][i - 2] * np.power(X[i - 2], 3.0) for i in range(2, 6)]
        )

        nu = self._theta(pressure, params) / temperature
        dthetadP = self._dtheta_dP(pressure, params)
        V_qh = (
            3.0
            * params["n"]
            * gas_constant
            * np.exp(-nu)
            / (1.0 - np.exp(-nu))
            * dthetadP
        )  # eq. 6

        f = np.sqrt(
            1.0
            + 2.0
            * params["b"][1]
            * (1.0 + params["delta"][1])
            * pressure
            / params["K_0"]
        )
        dIdP = (
            (1.0 + params["delta"][1])
            / (params["K_0"] * (1.0 + params["b"][1]))
            * np.exp((1.0 - f) / params["b"][1])
        )
        V_th = self._C_T(temperature, params) * dIdP

        # V = dG_c/dP + dG_qh/dP - C_T*(dI_P/dP)
        return V_c + V_qh + V_th

    def dvolume_dP(self, pressure, temperature, params):
        K_0 = params["K_0"]
        n = params["n"]
        b = params["b"][1]
        delta = params["delta"][1]

        dVc_dP = 0.0

        for i in range(2, 6):
            a_i = params["a"][i - 2]
            c_i = params["c"][i - 2]
            power_term = 1.0 + i * pressure / (3.0 * a_i * K_0)
            z_i = power_term ** (1.0 / i)
            x_i = 1.0 / (1.0 - a_i + a_i * z_i)
            dx_i_dP = -(x_i**2) / (3.0 * K_0) * power_term ** ((1.0 / i) - 1.0)
            dVc_dP_i = c_i * 3.0 * x_i**2 * dx_i_dP
            dVc_dP += dVc_dP_i

        dVc_dP *= params["V_0"]

        theta = self._theta(pressure, params)
        dthetadP = self._dtheta_dP(pressure, params)
        d2thetadP2 = self._d2theta_dP2(pressure, params)

        nu = theta / temperature
        exp_neg_nu = np.exp(-nu)
        denom = 1.0 - exp_neg_nu

        factor1 = -(exp_neg_nu) / (denom**2) * (dthetadP**2) / temperature
        factor2 = (exp_neg_nu) / denom * d2thetadP2

        dV_qh_dP = 3.0 * n * gas_constant * (factor1 + factor2)

        alpha = 2.0 * b * (1.0 + delta) / K_0
        f = np.sqrt(1.0 + alpha * pressure)
        exp_term = np.exp((1.0 - f) / b)

        numerator = self._C_T(temperature, params) * (1.0 + delta) ** 2 * exp_term
        denominator = K_0**2 * (1.0 + b) * f

        dV_th_dP = -numerator / denominator

        return dVc_dP + dV_qh_dP + dV_th_dP

    def pressure(self, temperature, volume, params):
        def _delta_volume(pressure):
            return self.volume(pressure, temperature, params) - volume

        try:
            sol = bracket(_delta_volume, 300.0e9, 1.0e5, ())
        except ValueError:
            raise Exception(
                "Cannot find a pressure, perhaps you are outside "
                "the range of validity for the equation of state?"
            )

        return brentq(_delta_volume, sol[0], sol[1])

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns the isothermal bulk modulus :math:`K_T` :math:`[Pa]`
        as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return -volume / self.dvolume_dP(pressure, temperature, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns the shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return 0.0

    def _Cp_1bar(self, temperature, params):
        # first, identify which of the piecewise segments we're in
        i = np.argmax(
            [T > temperature for T in list(zip(*params["gibbs_coefficients"]))[0]]
        )

        # select the appropriate coefficients
        coeffs = params["gibbs_coefficients"][i][1]
        Cp = -(
            coeffs[2]
            + 2.0 * coeffs[3] / temperature / temperature
            + 6.0 * coeffs[4] / (temperature * temperature * temperature)
            + 12.0 * coeffs[5] * np.power(temperature, -4.0)
            + 90.0 * coeffs[6] * np.power(temperature, -10.0)
            + 2.0 * coeffs[7] * temperature
            + 6.0 * coeffs[8] * temperature * temperature
            + 12.0 * coeffs[9] * temperature * temperature * temperature
            + 42.0 * coeffs[10] * np.power(temperature, 6.0)
            - 0.25 * coeffs[11] / np.sqrt(temperature)
            - coeffs[12] / temperature
        )
        return Cp

    def _S_1bar(self, temperature, params):
        # first, identify which of the piecewise segments we're in
        i = np.argmax(
            [T > temperature for T in list(zip(*params["gibbs_coefficients"]))[0]]
        )

        # select the appropriate coefficients
        coeffs = params["gibbs_coefficients"][i][1]
        S = -(
            coeffs[1]
            + coeffs[2] * (1.0 + np.log(temperature))
            - coeffs[3] / temperature / temperature
            - 2.0 * coeffs[4] / (temperature * temperature * temperature)
            - 3.0 * coeffs[5] * np.power(temperature, -4.0)
            - 9.0 * coeffs[6] * np.power(temperature, -10.0)
            + 2.0 * coeffs[7] * temperature
            + 3.0 * coeffs[8] * temperature * temperature
            + 4.0 * coeffs[9] * temperature * temperature * temperature
            + 7.0 * coeffs[10] * np.power(temperature, 6.0)
            + 0.5 * coeffs[11] / np.sqrt(temperature)
            + coeffs[12] / temperature
        )
        return S

    def _gibbs_1bar(self, temperature, params):
        # first, identify which of the piecewise segments we're in
        i = np.argmax(
            [T > temperature for T in list(zip(*params["gibbs_coefficients"]))[0]]
        )

        # select the appropriate coefficients
        coeffs = params["gibbs_coefficients"][i][1]
        gibbs = (
            coeffs[0]
            + coeffs[1] * temperature
            + coeffs[2] * temperature * np.log(temperature)
            + coeffs[3] / temperature
            + coeffs[4] / (temperature * temperature)
            + coeffs[5] / (temperature * temperature * temperature)
            + coeffs[6] * np.power(temperature, -9.0)
            + coeffs[7] * temperature * temperature
            + coeffs[8] * temperature * temperature * temperature
            + coeffs[9] * np.power(temperature, 4.0)
            + coeffs[10] * np.power(temperature, 7.0)
            + coeffs[11] * np.sqrt(temperature)
            + coeffs[12] * np.log(temperature)
        )
        return gibbs

    def _X(self, pressure, params):
        return [
            1.0
            / (
                1.0
                - params["a"][n - 2]
                + params["a"][n - 2]
                * np.power(
                    1.0
                    + float(n) / (3.0 * params["a"][n - 2]) * pressure / params["K_0"],
                    1.0 / float(n),
                )
            )
            for n in range(2, 6)
        ]  # eq. A2

    def _Gamma(self, n, an, Xn):
        def d(k, Xn):
            return (
                np.power(Xn, 3.0 - float(k)) * float(k) / (float(k) - 3.0)
                if k != 3
                else -3.0 * np.log(Xn)
            )  # eq. A9

        return (
            3.0
            * np.power(an, 1.0 - float(n))
            / float(n)
            * np.sum(
                [
                    binom(n, k) * np.power(an - 1.0, float(n - k)) * d(k, Xn)
                    for k in range(0, n + 1)
                ]
            )
        )  # eq. A9, CHECKED

    def _dGamma_dXT2(self, n, an, Xn):
        total = 0.0
        for k in range(0, n + 1):
            coeff = binom(n, k) * (an - 1.0) ** (n - k)
            if k != 3:
                term = (3.0 - k) * Xn ** (2.0 - k) * k / (k - 3.0)
            else:
                term = -3.0 / Xn
            total += coeff * term
        return 3.0 * an ** (1.0 - n) / float(n) * total

    def _d2Gamma_dXT2_2(self, n, an, Xn):
        total = 0.0
        for k in range(0, n + 1):
            coeff = binom(n, k) * (an - 1.0) ** (n - k)
            if k != 3:
                term = (2.0 - k) * (3.0 - k) * Xn ** (1.0 - k) * k / (k - 3.0)
            else:
                term = 3.0 / (Xn**2)
            total += coeff * term
        return 3.0 * an ** (1.0 - n) / float(n) * total

    def _theta(self, pressure, params):
        # Theta (for quasiharmonic term)
        ab2 = 1.0 / (3.0 * params["b"][0] - 1.0)
        K0b = params["K_0"] / (1.0 + params["delta"][0])  # eq. B1b
        XT2 = 1.0 / (
            1.0 - ab2 + ab2 * np.power(1.0 + 2.0 / (3.0 * ab2) * pressure / K0b, 0.5)
        )  # eq. 6 b of SE2015

        # eq. B1 (6 of SE2015)
        return params["theta_0"] * np.exp(
            params["grueneisen_0"]
            / (1.0 + params["delta"][0])
            * (self._Gamma(2, ab2, XT2) - self._Gamma(2, ab2, 1.0))
        )

    def _dtheta_dP(self, pressure, params):
        ab2 = 1.0 / (3.0 * params["b"][0] - 1.0)
        K0b = params["K_0"] / (1.0 + params["delta"][0])
        A = 2.0 / (3.0 * ab2 * K0b)
        y = np.sqrt(1.0 + A * pressure)
        XT2 = 1.0 / (1.0 - ab2 + ab2 * y)

        dG_dXT2 = self._dGamma_dXT2(2, ab2, XT2)

        # dXT2/dy
        dXT2_dy = -ab2 / (1.0 - ab2 + ab2 * y) ** 2

        # dy/dP
        dy_dP = 0.5 * A / np.sqrt(1.0 + A * pressure)

        # Prefactor
        B = params["grueneisen_0"] / (1.0 + params["delta"][0])

        # Theta
        theta_val = self._theta(pressure, params)

        return theta_val * B * dG_dXT2 * dXT2_dy * dy_dP

    def _d2theta_dP2(self, pressure, params):
        ab2 = 1.0 / (3.0 * params["b"][0] - 1.0)
        K0b = params["K_0"] / (1.0 + params["delta"][0])
        A = 2.0 / (3.0 * ab2 * K0b)
        y = np.sqrt(1.0 + A * pressure)
        dy_dP = 0.5 * A / np.sqrt(1.0 + A * pressure)
        d2y_dP2 = -0.25 * A * A / (1.0 + A * pressure) ** (3.0 / 2.0)

        XT2 = 1.0 / (1.0 - ab2 + ab2 * y)
        D = 1.0 - ab2 + ab2 * y
        dXT2_dy = -ab2 / (D**2)
        d2XT2_dy2 = 2.0 * ab2 * ab2 / (D**3)

        XT2_prime = dXT2_dy * dy_dP
        XT2_double_prime = d2XT2_dy2 * (dy_dP**2) + dXT2_dy * d2y_dP2

        # Gamma derivatives
        dGamma_dXT2 = self._dGamma_dXT2(2, ab2, XT2)
        d2Gamma_dXT2_2 = self._d2Gamma_dXT2_2(2, ab2, XT2)

        B = params["grueneisen_0"] / (1.0 + params["delta"][0])
        theta_val = self._theta(pressure, params)
        dtheta_dP = self._dtheta_dP(pressure, params)

        return (
            theta_val
            * B
            * (d2Gamma_dXT2_2 * XT2_prime**2 + dGamma_dXT2 * XT2_double_prime)
            + dtheta_dP * B * dGamma_dXT2 * XT2_prime
        )

    def _interpolating_function(self, pressure, params):
        f = np.sqrt(
            1.0
            + 2.0
            * params["b"][1]
            * (1.0 + params["delta"][1])
            * pressure
            / params["K_0"]
        )

        # eq. D2 (9 of SE2015)
        return (
            1.0
            / (1.0 + params["b"][1])
            * (params["b"][1] + f)
            * np.exp((1.0 - f) / params["b"][1])
        )

    def _gibbs_qh(self, temperature, theta, n):
        return (
            3.0
            * n
            * gas_constant
            * temperature
            * np.log(1.0 - np.exp(-theta / temperature))
        )  # eq. 5

    def _S_qh(self, temperature, theta, n):
        nu = theta / temperature
        return (
            3.0
            * n
            * gas_constant
            * (nu / (np.exp(nu) - 1.0) - np.log(1.0 - np.exp(-nu)))
        )

    def _C_T(self, temperature, params):
        # C, which is the (G_qh(t,p0) - G_sgte(t,p0)) term
        G_SGTE = self._gibbs_1bar(temperature, params)
        G_qh0 = self._gibbs_qh(temperature, params["theta_0"], params["n"])
        if temperature < params["T_0"]:
            C_T = (
                temperature * temperature / (2.0 * params["T_0"]) * params["delta_Cpr"]
            )

        else:
            C_T = (
                (G_qh0 - G_SGTE)
                + params["delta_Gr"]
                - (temperature - params["T_0"]) * params["delta_Sr"]
                + (temperature - params["T_0"] / 2.0) * params["delta_Cpr"]
            )

        return C_T

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy of the mineral. :math:`[J/mol]`
        """
        # Cold compression term, eq. A8
        X = self._X(pressure, params)
        G_c = (
            params["K_0"]
            * params["V_0"]
            * np.sum(
                [
                    params["c"][n - 2]
                    * (
                        self._Gamma(n, params["a"][n - 2], X[n - 2])
                        - self._Gamma(n, params["a"][n - 2], 1.0)
                    )
                    for n in range(2, 6)
                ]
            )
        )

        # G_SGTE
        G_SGTE = self._gibbs_1bar(temperature, params)

        # G_qh
        theta = self._theta(pressure, params)
        G_qh = self._gibbs_qh(temperature, theta, params["n"])
        G_qh0 = self._gibbs_qh(temperature, params["theta_0"], params["n"])

        C_T = self._C_T(temperature, params)
        I_P = self._interpolating_function(pressure, params)
        return G_SGTE + G_c + G_qh - G_qh0 + C_T * (1.0 - I_P)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy of the mineral. :math:`[J/K/mol]`
        """

        S_SGTE = self._S_1bar(temperature, params)

        # S_qh
        theta = self._theta(pressure, params)
        S_qh = self._S_qh(temperature, theta, params["n"])
        S_qh0 = self._S_qh(temperature, params["theta_0"], params["n"])

        # dCdT, which is the (S_qh(t,p0) - S_sgte(t,p0)) term
        if temperature < params["T_0"]:
            dC_TdT = temperature / params["T_0"] * params["delta_Cpr"]

        else:
            dC_TdT = -(S_qh0 - S_SGTE) - params["delta_Sr"] + params["delta_Cpr"]

        I_P = self._interpolating_function(pressure, params)
        return S_SGTE + S_qh - S_qh0 - dC_TdT * (1.0 - I_P)

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the molar isobaric heat capacity :math:`[J/K/mol]`.
        For now, this is calculated by numerical differentiation.
        """
        dT = 0.1
        if temperature < dT / 2.0:
            return 0.0
        else:
            dS = self.entropy(
                pressure, temperature + dT / 2.0, volume, params
            ) - self.entropy(pressure, temperature - dT / 2.0, volume, params)
            return temperature * dS / dT

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns the volumetric thermal expansivity :math:`[1/K]`.
        For now, this is calculated by numerical differentiation.
        """
        dT = 0.1
        if temperature < dT / 2.0:
            return 0.0
        else:
            dV = self.volume(pressure, temperature + dT / 2.0, params) - self.volume(
                pressure, temperature - dT / 2.0, params
            )
            return dV / dT / volume

    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns the grueneisen parameter.
        This is a dependent thermodynamic variable in this equation of state.
        """
        Cv = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        if Cv == 0.0:
            return 0.0
        else:
            return (
                self.thermal_expansivity(pressure, temperature, volume, params)
                * self.isothermal_bulk_modulus_reuss(
                    pressure, temperature, volume, params
                )
                * self.volume(pressure, temperature, params)
            ) / Cv

    def calculate_transformed_parameters(self, params):
        """
        This function calculates the "c" parameters of the :cite:`Brosh2007`
        equation of state.
        """
        Zs = pkgutil.get_data("burnman", "data/input_masses/atomic_numbers.dat")
        Zs = Zs.decode("ascii").split("\n")
        Z = {
            str(sl[0]): int(sl[1])
            for sl in [line.split() for line in Zs if len(line) > 0 and line[0] != "#"]
        }

        nZs = [(n_at, float(Z[el])) for (el, n_at) in params["formula"].items()]

        # eq. A2 at 300 TPa
        X3_300TPa = [
            np.power(
                1.0
                - params["a"][i - 2]
                + params["a"][i - 2]
                * np.power(
                    (
                        1.0
                        + float(i)
                        / (3.0 * params["a"][i - 2])
                        * 300.0e12
                        / params["K_0"]
                    ),
                    1.0 / float(i),
                ),
                -3.0,
            )
            for i in range(2, 6)
        ]

        # eq. A2 at 330 TPa
        X3_330TPa = [
            np.power(
                1.0
                - params["a"][i - 2]
                + params["a"][i - 2]
                * np.power(
                    (
                        1.0
                        + float(i)
                        / (3.0 * params["a"][i - 2])
                        * 330.0e12
                        / params["K_0"]
                    ),
                    1.0 / float(i),
                ),
                -3.0,
            )
            for i in range(2, 6)
        ]

        # eq. A6a, m^3/mol
        V_QSM_300TPa = (
            np.sum(
                [
                    n_at
                    * (
                        0.02713
                        * np.exp(
                            0.97626 * np.log(Zi) - 0.057848 * np.log(Zi) * np.log(Zi)
                        )
                    )
                    for (n_at, Zi) in nZs
                ]
            )
            * 1.0e-6
        )

        # eq. A6b, m^3/mol
        V_QSM_330TPa = (
            np.sum(
                [
                    n_at
                    * (
                        0.025692
                        * np.exp(
                            0.97914 * np.log(Zi) - 0.057741 * np.log(Zi) * np.log(Zi)
                        )
                    )
                    for (n_at, Zi) in nZs
                ]
            )
            * 1.0e-6
        )

        A = np.array(
            [
                [1.0, 1.0, 1.0, 1.0],  # eq A3
                [0.0, 6.0, 8.0, 9.0],  # eq A4
                X3_300TPa,  # eq A5a
                X3_330TPa,
            ]
        )  # eq A5b

        b = np.array(
            [1.0, 8.0, V_QSM_300TPa / params["V_0"], V_QSM_330TPa / params["V_0"]]
        )

        # does not quite reproduce the published values of c
        # A.c consistently gives b[2], b[3] ~1% larger than Brosh
        return np.linalg.solve(A, b)

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        params["T_0"] = 298.15
        if "P_0" not in params:
            params["P_0"] = 1.0e5

        if "a" not in params:
            params["a"] = [
                (float(i) - 1.0) / (3.0 * params["Kprime_0"] - 1.0) for i in range(2, 6)
            ]  # eq. A2

        if "c" not in params:
            params["c"] = self.calculate_transformed_parameters(params)

        # Calculate reference values for gibbs free energy and heat capacity
        nur = params["theta_0"] / params["T_0"]
        G_qhr = self._gibbs_qh(params["T_0"], params["theta_0"], params["n"])
        S_qhr = self._S_qh(params["T_0"], params["theta_0"], params["n"])
        Cp_qhr = (
            3.0
            * params["n"]
            * gas_constant
            * nur
            * nur
            * np.exp(nur)
            / np.power(np.exp(nur) - 1.0, 2.0)
        )

        G_SGTEr = self._gibbs_1bar(params["T_0"], params)
        S_SGTEr = self._S_1bar(params["T_0"], params)
        Cp_SGTEr = self._Cp_1bar(params["T_0"], params)

        params["delta_Cpr"] = Cp_SGTEr - Cp_qhr
        params["delta_Gr"] = G_SGTEr - G_qhr
        params["delta_Sr"] = S_SGTEr - S_qhr

        # Check that all the required keys are in the dictionary
        expected_keys = [
            "gibbs_coefficients",
            "V_0",
            "K_0",
            "Kprime_0",
            "theta_0",
            "grueneisen_0",
            "delta",
            "b",
        ]
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
