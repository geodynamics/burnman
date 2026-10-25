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

    This equation of state includes contributions from
    cold compression, quasiharmonic lattice vibrations,
    and thermal excitation of electrons.

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``gibbs_coefficients``
          - Coefficients for the piecewise polynomial fit to the Gibbs free energy at 1 bar.
            Format is a list of lists, where each sublist contains the maximum temperature for
            that segment and the coefficients in the following order:
            [const, T, T*ln(T), T^(-1), T^(-2), T^(-3), T^(-9), T^2, T^3, T^4, T^7, T^(1/2), ln(T)]
            The final segment is cautiously extrapolated above its maximum temperature.
          - varies
        * - ``V_0``
          - Reference volume.
          - :math:`\\text{m}^3`
        * - ``K_0``
          - Reference bulk modulus.
          - :math:`\\text{Pa}`
        * - ``n``
          - Number of atoms per formula unit.
          - Dimensionless
        * - ``a``
          - Parameter array for cold compression (length 4 array).
          - Dimensionless
        * - ``b``
          - Parameter array for thermal excitation (length 2 array).
          - Dimensionless
        * - ``c``
          - [OPTIONAL] Parameter array for cold compression (length 4 array).
            If you do not know how to calculate it, do not add it to the
            parameters dictionary.
            It will be calculated from the other input parameters.
          - Dimensionless
        * - ``delta``
          - Parameter array for thermal excitation (length 2 array).
          - Dimensionless
        * - ``theta_0``
          - Reference Einstein temperature.
          - :math:`\\text{K}`
        * - ``grueneisen_0``
          - Reference Grüneisen parameter.
          - Dimensionless

    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        X = self._X(pressure, params)
        V_c = params["V_0"] * np.dot(params["c"], np.power(X, 3.0))

        if temperature == 0.0:
            V_qh = 0.0
        else:
            nu = self._theta(pressure, params) / temperature
            V_qh = (
                3.0
                * params["n"]
                * gas_constant
                * np.exp(-nu)
                / -np.expm1(-nu)
                * self._dtheta_dP(pressure, params)
            )  # eq. 6

        V_th = -self._C_T(temperature, params) * self._dI_dP(pressure, params)

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

        if temperature == 0.0:
            dV_qh_dP = 0.0
        else:
            nu = theta / temperature
            exp_neg_nu = np.exp(-nu)
            denom = -np.expm1(-nu)
            factor1 = -exp_neg_nu / denom**2 * dthetadP**2 / temperature
            factor2 = exp_neg_nu / denom * d2thetadP2
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

    def _sgte_coefficients(self, temperature, params):
        # Interior boundaries select the next segment. Extrapolate the last
        # segment above its upper limit, rather than wrapping to the first.
        for upper_temperature, coefficients in params["gibbs_coefficients"]:
            if temperature < upper_temperature:
                return coefficients
        return params["gibbs_coefficients"][-1][1]

    def _Cp_1bar(self, temperature, params):
        coeffs = self._sgte_coefficients(temperature, params)
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
        coeffs = self._sgte_coefficients(temperature, params)
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
        coeffs = self._sgte_coefficients(temperature, params)
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

    def _dlogtheta_dP(self, pressure, params):
        # Appendix B: d(log(theta))/dP = gamma_0 * (X_2^T)^3 / K_0.
        a = 1.0 / (3.0 * params["b"][0] - 1.0)
        K = params["K_0"] / (1.0 + params["delta"][0])
        X = 1.0 / (1.0 - a + a * np.sqrt(1.0 + 2.0 * pressure / (3.0 * a * K)))
        return params["grueneisen_0"] / params["K_0"] * X**3

    def _dtheta_dP(self, pressure, params):
        return self._theta(pressure, params) * self._dlogtheta_dP(pressure, params)

    def _d2theta_dP2(self, pressure, params):
        a = 1.0 / (3.0 * params["b"][0] - 1.0)
        K = params["K_0"] / (1.0 + params["delta"][0])
        y = np.sqrt(1.0 + 2.0 * pressure / (3.0 * a * K))
        X = 1.0 / (1.0 - a + a * y)
        dlogtheta_dP = params["grueneisen_0"] / params["K_0"] * X**3
        return self._theta(pressure, params) * (
            dlogtheta_dP**2 - dlogtheta_dP * X / (K * y)
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

    def _dI_dP(self, pressure, params):
        b = params["b"][1]
        delta = params["delta"][1]
        f = np.sqrt(1.0 + 2.0 * b * (1.0 + delta) * pressure / params["K_0"])
        return -(1.0 + delta) / (params["K_0"] * (1.0 + b)) * np.exp((1.0 - f) / b)

    def _gibbs_qh(self, temperature, theta, n):
        if temperature == 0.0:
            return 0.0
        return (
            3.0
            * n
            * gas_constant
            * temperature
            * np.log(-np.expm1(-theta / temperature))
        )  # eq. 5

    def _S_qh(self, temperature, theta, n):
        if temperature == 0.0:
            return 0.0
        nu = theta / temperature
        denominator = -np.expm1(-nu)
        return (
            3.0
            * n
            * gas_constant
            * (nu * np.exp(-nu) / denominator - np.log(denominator))
        )

    def _Cp_qh(self, temperature, theta, n):
        if temperature == 0.0:
            return 0.0
        nu = theta / temperature
        return 3.0 * n * gas_constant * nu**2 * np.exp(-nu) / np.expm1(-nu) ** 2

    def _C_T(self, temperature, params):
        # C, which is the (G_qh(t,p0) - G_sgte(t,p0)) term
        if temperature < params["T_0"]:
            C_T = (
                temperature * temperature / (2.0 * params["T_0"]) * params["delta_Cpr"]
            )

        else:
            G_SGTE = self._gibbs_1bar(temperature, params)
            G_qh0 = self._gibbs_qh(temperature, params["theta_0"], params["n"])
            C_T = (
                (G_qh0 - G_SGTE)
                + params["delta_Gr"]
                - (temperature - params["T_0"]) * params["delta_Sr"]
                + (temperature - params["T_0"] / 2.0) * params["delta_Cpr"]
            )

        return C_T

    def _dC_dT(self, temperature, params):
        if temperature < params["T_0"]:
            return temperature / params["T_0"] * params["delta_Cpr"]
        return (
            self._S_1bar(temperature, params)
            - self._S_qh(temperature, params["theta_0"], params["n"])
            - params["delta_Sr"]
            + params["delta_Cpr"]
        )

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

        theta = self._theta(pressure, params)
        G_qh = self._gibbs_qh(temperature, theta, params["n"])
        C_T = self._C_T(temperature, params)
        I_P = self._interpolating_function(pressure, params)

        # Eq. 4 with Appendix C's continuation.
        reference = (
            params["delta_Gr"]
            + params["T_0"] * params["delta_Sr"]
            - params["T_0"] / 2.0 * params["delta_Cpr"]
            + temperature * (params["delta_Cpr"] - params["delta_Sr"])
        )
        return G_c + G_qh - C_T * I_P + reference

    def entropy(self, pressure, temperature, volume, params):
        """Returns the molar entropy :math:`[J/K/mol]`."""
        theta = self._theta(pressure, params)
        return (
            self._S_qh(temperature, theta, params["n"])
            + self._dC_dT(temperature, params)
            * self._interpolating_function(pressure, params)
            + params["delta_Sr"]
            - params["delta_Cpr"]
        )

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """Returns the molar isobaric heat capacity :math:`[J/K/mol]`."""
        Cp_qh = self._Cp_qh(temperature, self._theta(pressure, params), params["n"])
        if temperature < params["T_0"]:
            correction = temperature / params["T_0"] * params["delta_Cpr"]
        else:
            correction = self._Cp_1bar(temperature, params) - self._Cp_qh(
                temperature, params["theta_0"], params["n"]
            )
        return Cp_qh + correction * self._interpolating_function(pressure, params)

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """Returns the volumetric thermal expansivity :math:`[1/K]`."""
        Cp_qh = self._Cp_qh(temperature, self._theta(pressure, params), params["n"])
        return (
            Cp_qh * self._dlogtheta_dP(pressure, params)
            - self._dC_dT(temperature, params) * self._dI_dP(pressure, params)
        ) / volume

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
        G_qhr = self._gibbs_qh(params["T_0"], params["theta_0"], params["n"])
        S_qhr = self._S_qh(params["T_0"], params["theta_0"], params["n"])
        Cp_qhr = self._Cp_qh(params["T_0"], params["theta_0"], params["n"])

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
