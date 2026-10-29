# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU GPL v2 or later.


import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
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

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``E_0``
          - Reference internal energy.
          - :math:`\\text{J/mol}`
        * - ``V_0``
          - Reference volume.
          - :math:`\\text{m}^3/\\text{mol}`
        * - ``P_0``
          - Reference pressure on the principal isentrope.
          - :math:`\\text{Pa}`
        * - ``K_S``
          - Reference isentropic bulk modulus.
          - :math:`\\text{Pa}`
        * - ``Kprime_S``
          - First pressure derivative of the isentropic bulk modulus.
          - :math:`\\text{unitless}`
        * - ``Kprime_prime_S``
          - Second pressure derivative of the isentropic bulk modulus.
          - :math:`\\text{Pa}^{-1}`
        * - ``T_0``
          - Reference temperature on the principal isentrope.
          - :math:`\\text{K}`
        * - ``S_0``
          - Reference entropy on the principal isentrope.
          - :math:`\\text{J/K/mol}`
        * - ``n``
          - Number of atoms per formula unit.
          - :math:`\\text{unitless}`
        * - ``grueneisen_0``
          - Reference Gruneisen parameter on the principal isentrope.
          - :math:`\\text{unitless}`
        * - ``grueneisen_n``
          - Volume dependence of the Gruneisen parameter on the principal isentrope.
          - :math:`\\text{unitless}`
        * - ``grueneisen_prime``
          - Internal-energy derivative of the Gruneisen parameter at reference density.
          - :math:`\\text{mol/J}`
        * - ``theta``
          - Temperature scale of the potential heat-capacity contribution.
          - :math:`\\text{K}`
        * - ``a``
          - Electronic heat capacity parameter A (length 2 array).
          - :math:`\\text{J/(mol K)}`
        * - ``b``
          - Electronic heat capacity parameter B (length 2 array).
          - :math:`\\text{J/(mol K}^{8/5}\\text{)}`
        * - ``Theta``
          - Volume dependence parameters for electronic heat capacity (length 2 array).
          - K for the first entry; dimensionless for the second.
        * - ``lmda``
          - Potential heat capacity parameter lambda (length 3 array).
          - :math:`\\text{J/(mol K)}`
        * - ``xi_0``
          - Reference xi parameter.
          - :math:`\\text{J/(mol K)}`
        * - ``F``
          - Reference F parameter (length 2 array).
          - :math:`\\text{unitless}`

    """

    @staticmethod
    def _check_state(temperature, volume):
        if np.isscalar(temperature) and np.isscalar(volume):
            if (
                temperature > 0.0
                and volume > 0.0
                and np.isfinite(temperature)
                and np.isfinite(volume)
            ):
                return
        if not np.all(np.isfinite(temperature) & (np.asarray(temperature) > 0.0)):
            raise ValueError("AA requires positive, finite temperature.")
        if not np.all(np.isfinite(volume) & (np.asarray(volume) > 0.0)):
            raise ValueError("AA requires positive, finite volume.")

    @staticmethod
    def _integrated_power(exponent, log_rhofrac):
        """Integral of exp(exponent*l), including its removable singularity."""
        if exponent == 0.0:
            return log_rhofrac
        return np.expm1(exponent * log_rhofrac) / exponent

    @staticmethod
    def _thermal_functions(
        temperature, reference_temperature, Theta, theta, helmholtz=True
    ):
        """Caloric basis functions at fixed volume and reference temperature.

        These are the same Cv and Cv/T integrals as in Appendix A. log1p,
        expm1 and the arctangent difference retain precision close to T_s.
        """
        T, Ts = temperature, reference_temperature
        dT = T - Ts
        lnT = np.log1p(dT / Ts)
        ln_pot = np.log1p(dT / (theta + Ts))
        S_e = 0.5 * np.log1p(dT * (T + Ts) / (Theta * Theta + Ts * Ts))
        S_b = (5.0 / 3.0) * np.power(Ts, 0.6) * np.expm1(0.6 * lnT)
        functions = {
            "lnT": lnT,
            "ln_pot": ln_pot,
            "S_e": S_e,
            "S_b": S_b,
        }
        if not helmholtz:
            return functions
        atan_difference = np.arctan2(Theta * dT, Theta * Theta + T * Ts)
        functions.update(
            {
                "F_kin": dT - T * lnT,
                "F_e": dT - Theta * atan_difference - T * S_e,
                "F_b": (5.0 / 8.0) * np.power(Ts, 1.6) * np.expm1(1.6 * lnT) - T * S_b,
                "F_lmda": dT - (theta + T) * ln_pot,
                "F_xi": (theta + T) * ln_pot - T * lnT,
                "atan_difference": atan_difference,
            }
        )
        return functions

    @staticmethod
    def _heat_capacity(temperature, A, B, Theta, lmda, xi, theta, C_kin):
        """Kinetic, electronic (A1) and potential (A15) contributions."""
        return (
            C_kin
            + A * temperature**2 / (Theta**2 + temperature**2)
            + B * np.power(temperature, 0.6)
            + (lmda * temperature + xi * theta) / (theta + temperature)
        )

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
        F = expit((params["F"][0] - rhofrac) / params["F"][1])  # A18
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

        rhofrac, _, ksi1, ksi2 = self._rhofracxksis(volume, params)
        log_rhofrac = np.log(rhofrac)

        # Equation B6 -- B10
        a1 = ksi2 / 8.0
        a2 = (ksi1 + 3.0 * ksi2) / 6.0
        a3 = (1.0 + 2.0 * ksi1 + 3.0 * ksi2) / 4.0
        a4 = (1.0 + ksi1 + ksi2) / 2.0
        a5 = (6.0 + 4.0 * ksi1 + 3.0 * ksi2) / 24.0
        n = params["grueneisen_n"]

        # B5 written as integrals of powers in log density, without division
        # by zero when 3*n+j vanishes. P_0 work also enters E_s in gamma_s.
        integral = (
            4.5
            * params["V_0"]
            * params["K_S"]
            * sum(
                coefficient * self._integrated_power(n + power / 3.0, log_rhofrac)
                for coefficient, power in [
                    (a1, 8.0),
                    (-a2, 6.0),
                    (a3, 4.0),
                    (-a4, 2.0),
                    (a5, 0.0),
                ]
            )
        )
        integral += (
            params["P_0"]
            * params["V_0"]
            * (
                self._integrated_power(n, log_rhofrac)
                - self._integrated_power(n - 1.0, log_rhofrac)
            )
        )
        return params["T_0"] * np.exp(
            params["grueneisen_0"] * log_rhofrac + params["grueneisen_prime"] * integral
        )

    def _reference_pressure(self, volume, params):
        """
        Pressure along the reference isentrope (Eq. 17), including P_0.
        """
        rhofrac, _, ksi1, ksi2 = self._rhofracxksis(volume, params)
        q = np.expm1((2.0 / 3.0) * np.log(rhofrac))
        return params["P_0"] + 1.5 * params["K_S"] * np.power(
            rhofrac, 5.0 / 3.0
        ) * q * (1.0 - ksi1 * q + ksi2 * q * q)

    def _reference_bulk_modulus(self, volume, params):
        rhofrac, _, ksi1, ksi2 = self._rhofracxksis(volume, params)
        q = np.expm1((2.0 / 3.0) * np.log(rhofrac))
        return (
            1.5
            * params["K_S"]
            * np.power(rhofrac, 5.0 / 3.0)
            * (
                (5.0 / 3.0) * q * (1.0 - ksi1 * q + ksi2 * q * q)
                + (2.0 / 3.0) * (q + 1.0) * (1.0 - 2.0 * ksi1 * q + 3.0 * ksi2 * q * q)
            )
        )

    def _isentropic_energy_change(self, volume, params):
        """
        Birch Murnaghan equation of state expression for the energy change along an isentrope
        """
        rhofrac, _, ksi1, ksi2 = self._rhofracxksis(volume, params)
        log_rhofrac = np.log(rhofrac)
        q = np.expm1((2.0 / 3.0) * log_rhofrac)
        # Eq. 21 factored in strain; avoid cancellation at V=V_0.
        return 4.5 * params["V_0"] * params["K_S"] * q * q * (
            0.25 - ksi1 * q / 6.0 + ksi2 * q * q / 8.0
        ) - params["P_0"] * params["V_0"] * np.expm1(-log_rhofrac)

    def _molar_heat_capacity_v(self, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        self._check_state(temperature, volume)
        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)

        return self._heat_capacity(
            temperature,
            A,
            B,
            Theta,
            lmda,
            xi,
            params["theta"],
            1.5 * params["n"] * gas_constant,
        )

    def entropy(self, pressure, temperature, volume, params):
        self._check_state(temperature, volume)
        Trel = self._reference_temperature(volume, params)
        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)
        f = self._thermal_functions(
            temperature, Trel, Theta, params["theta"], helmholtz=False
        )
        return (
            params["S_0"]
            + (1.5 * params["n"] * gas_constant + xi) * f["lnT"]
            + (lmda - xi) * f["ln_pot"]
            + A * f["S_e"]
            + B * f["S_b"]
        )

    def _helmholtz_thermal(self, temperature, volume, params, Trel=None):
        if Trel is None:
            self._check_state(temperature, volume)
            Trel = self._reference_temperature(volume, params)
        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)
        f = self._thermal_functions(temperature, Trel, Theta, params["theta"])
        return (
            -params["S_0"] * (temperature - Trel)
            + 1.5 * params["n"] * gas_constant * f["F_kin"]
            + (A * f["F_e"] + B * f["F_b"] + lmda * f["F_lmda"] + xi * f["F_xi"])
        )

    def _helmholtz_energy(self, temperature, volume, params):
        """
        Returns the Helmholtz energy at the temperature and volume of the mineral [J/mol]
        F = E - TS
        """
        self._check_state(temperature, volume)
        Eref = params["E_0"] + self._isentropic_energy_change(volume, params)
        Tref = self._reference_temperature(volume, params)
        Fref = Eref - Tref * params["S_0"]
        Fth = self._helmholtz_thermal(temperature, volume, params, Tref)
        return Fref + Fth

    def _thermal_derivatives(
        self, temperature, volume, params, bulk_modulus=True, thermal_pressure=True
    ):
        """Return P, K_T and P_T from the existing Helmholtz potential.

        Subscripts l denote derivatives with respect to l=ln(V_0/V).
        Coefficient derivatives hold T and T_s fixed. Moving the lower
        integration limit T_s supplies the explicit Cv_s boundary terms.
        No off-isentrope use of the paper's inconsistent Eq. 23 is made.
        Unrequested derivatives are returned as None and are not calculated.
        """
        self._check_state(temperature, volume)
        r = params["V_0"] / volume
        T, Ts = temperature, self._reference_temperature(volume, params)
        dT = T - Ts
        A, B, Theta = self._ABTheta(volume, params)
        xi = params["xi_0"] * np.power(r, -0.6)
        theta = params["theta"]
        A_l = -params["a"][1] / r
        B_l = -2.0 * params["b"][1] / (r * r)
        Theta_l = params["Theta"][1] * Theta
        xi_l = -0.6 * xi

        logistic = expit((params["F"][0] - r) / params["F"][1])
        logistic_l = -logistic * (1.0 - logistic) * r / params["F"][1]
        c = params["lmda"][2] * np.power(r, 0.4)
        h = params["lmda"][0] + params["lmda"][1] * r + c
        h_l = params["lmda"][1] * r + 0.4 * c
        lmda_l = logistic_l * h + logistic * h_l
        # Reuse these coefficients rather than evaluating _lambdaxi again.
        lmda = logistic * h

        f = self._thermal_functions(T, Ts, Theta, theta)
        denominator = Theta * Theta + Ts * Ts
        # First Theta derivative of the electronic F/A basis.
        Fe_Theta = -f["atan_difference"] + Theta * dT / denominator
        F_l = (
            A_l * f["F_e"]
            + A * Theta_l * Fe_Theta
            + B_l * f["F_b"]
            + (lmda_l * f["F_lmda"] + xi_l * f["F_xi"])
        )
        Cv_s = self._heat_capacity(
            Ts, A, B, Theta, lmda, xi, theta, 1.5 * params["n"] * gas_constant
        )
        Ps = self._reference_pressure(volume, params)
        Es = self._isentropic_energy_change(volume, params)
        g = params["grueneisen_prime"] * np.power(r, params["grueneisen_n"])
        gamma_s = params["grueneisen_0"] + g * Es
        F_l += gamma_s * Cv_s * dT
        P = Ps + F_l / volume
        aKT = None
        if thermal_pressure:
            S_l = (
                A_l * f["S_e"]
                + A
                * Theta_l
                * Theta
                * (1.0 / (Theta * Theta + T * T) - 1.0 / denominator)
                + B_l * f["S_b"]
                + (lmda_l - xi_l) * f["ln_pot"]
                + xi_l * f["lnT"]
            )
            aKT = (gamma_s * Cv_s - S_l) / volume
        if not bulk_modulus:
            return P, None, aKT

        # Second volume derivatives are only needed for K_T, not for the
        # pressure evaluations used by the volume solver.
        logistic_ll = (
            logistic
            * (1.0 - logistic)
            * (1.0 - 2.0 * logistic)
            * (r / params["F"][1]) ** 2
            + logistic_l
        )
        h_ll = params["lmda"][1] * r + 0.16 * c
        lmda_ll = logistic_ll * h + 2.0 * logistic_l * h_l + logistic * h_ll
        Fe_ThetaTheta = (
            T / (Theta * Theta + T * T)
            - Ts / denominator
            + dT * (Ts * Ts - Theta * Theta) / denominator**2
        )
        F_ll = (
            -A_l * f["F_e"]
            + (2.0 * A_l * Theta_l + A * params["Theta"][1] * Theta_l) * Fe_Theta
            + A * Theta_l**2 * Fe_ThetaTheta
            - 2.0 * B_l * f["F_b"]
            + lmda_ll * f["F_lmda"]
            + 0.36 * xi * f["F_xi"]
        )
        Cv_l = (
            A_l * Ts * Ts / denominator
            - 2.0 * A * Ts * Ts * Theta * Theta_l / denominator**2
            + (B_l * np.power(Ts, 0.6) + (lmda_l * Ts + xi_l * theta) / (theta + Ts))
        )
        Cv_T = (
            2.0 * A * Ts * Theta * Theta / denominator**2
            + 0.6 * B * np.power(Ts, -0.4)
            + (theta * (lmda - xi) / (theta + Ts) ** 2)
        )
        Ks = self._reference_bulk_modulus(volume, params)
        gamma_l = g * (params["grueneisen_n"] * Es + volume * Ps)

        F_ll += (
            2.0 * gamma_s * dT * Cv_l
            + gamma_s**2 * Ts * dT * Cv_T
            + Cv_s * (dT * gamma_l - Ts * gamma_s**2)
        )
        return P, Ks + (F_l + F_ll) / volume, aKT

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume [Pa]
        """
        return self._thermal_derivatives(
            temperature, volume, params, bulk_modulus=False, thermal_pressure=False
        )[0]

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """

        self._check_state(temperature, params["V_0"])
        if not np.isfinite(pressure):
            raise ValueError("AA requires finite pressure.")
        V_0 = params["V_0"]

        def delta_pressure(relative_volume):
            return pressure - self.pressure(temperature, relative_volume * V_0, params)

        # Keep the existing bracket; solve in V/V_0 so accuracy is independent
        # of formula-unit size and the small SI molar-volume scale.
        return V_0 * brentq(
            delta_pressure, 0.1, 2.0, xtol=1.0e-14, rtol=4.0 * np.finfo(float).eps
        )

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        return self._thermal_derivatives(temperature, volume, params)[1]

    def _aKT(self, temperature, volume, params):
        """
        Returns the product of the thermal expansivity and isothermal bulk modulus. :math:`[Pa/K]`
        """

        return self._thermal_derivatives(
            temperature, volume, params, bulk_modulus=False
        )[2]

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns the thermal expansivity. :math:`[1/K]`
        """

        _, K_T, aK_T = self._thermal_derivatives(temperature, volume, params)
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

        _, K_T, aK_T = self._thermal_derivatives(temperature, volume, params)
        C_v = self._molar_heat_capacity_v(temperature, volume, params)
        return C_v + temperature * volume * aK_T * aK_T / K_T

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        return self._molar_heat_capacity_v(temperature, volume, params)

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        return self._grueneisen_parameter(pressure, temperature, volume, params)

    def isentropic_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        _, K_T, aK_T = self._thermal_derivatives(temperature, volume, params)
        C_v = self._molar_heat_capacity_v(temperature, volume, params)
        return K_T + temperature * volume * aK_T * aK_T / C_v

    def helmholtz_energy(self, pressure, temperature, volume, params):
        return self._helmholtz_energy(temperature, volume, params)

    def molar_internal_energy(self, pressure, temperature, volume, params):
        return self._helmholtz_energy(
            temperature, volume, params
        ) + temperature * self.entropy(pressure, temperature, volume, params)

    def enthalpy(self, pressure, temperature, volume, params):
        return (
            self.molar_internal_energy(pressure, temperature, volume, params)
            + pressure * volume
        )

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        F + PV
        """
        F = self._helmholtz_energy(temperature, volume, params)
        return F + pressure * volume

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

        scalar_keys = [
            "P_0",
            "T_0",
            "S_0",
            "E_0",
            "V_0",
            "K_S",
            "Kprime_S",
            "Kprime_prime_S",
            "molar_mass",
            "n",
            "grueneisen_0",
            "grueneisen_prime",
            "grueneisen_n",
            "theta",
            "xi_0",
        ]
        array_sizes = {"a": 2, "b": 2, "Theta": 2, "lmda": 3, "F": 2}
        for key in scalar_keys + list(array_sizes):
            if key not in params:
                raise KeyError("params object missing parameter : " + key)
            shape = (array_sizes[key],) if key in array_sizes else ()
            try:
                valid = (
                    np.shape(params[key]) == shape
                    and np.isrealobj(params[key])
                    and np.all(np.isfinite(params[key]))
                )
            except TypeError:
                valid = False
            if not valid:
                raise ValueError(
                    "AA parameter "
                    + key
                    + " must have shape "
                    + str(shape)
                    + " and contain finite real values."
                )
        for key in ["T_0", "V_0", "K_S", "molar_mass", "n", "theta"]:
            if params[key] <= 0.0:
                raise ValueError("AA parameter " + key + " must be positive.")
        if params["Theta"][0] <= 0.0 or params["F"][1] <= 0.0:
            raise ValueError("AA requires positive Theta[0] and F[1].")

        if params["molar_mass"] < 0.001 or params["molar_mass"] > 10.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)
        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
