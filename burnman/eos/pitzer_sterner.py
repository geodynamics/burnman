# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import scipy.optimize as opt

from . import debye
from . import equation_of_state as eos
from .. import constants


class PitzerSterner(eos.EquationOfState):
    """
    Pitzer and Sterner fluid equation of state (1994, doi:10.1063/1.467624).

    Ported from BurnMan's water_filter branch (commit ffcb3688). The stable
    density-root selection follows burnman_cpp's water implementation.
    Parameters may specify either the original Debye thermal reference or
    an ideal-gas reference with H_0, S_0 and Holland-Powell Cp coefficients.
    The latter is required when combining water with HP mineral datasets.
    """

    def _compute_cs(self, c_coeffs, temperature):
        cs = np.empty(10)
        for i, ci in enumerate(c_coeffs):
            cs[i] = (
                ci[0] * np.power(temperature, -4.0)
                + ci[1] * np.power(temperature, -2.0)
                + ci[2] * np.power(temperature, -1.0)
                + ci[3]
                + ci[4] * temperature
                + ci[5] * np.power(temperature, 2)
            )
        return cs

    def _compute_dcdTs(self, c_coeffs, temperature):
        cs = np.empty(10)
        for i, ci in enumerate(c_coeffs):
            cs[i] = (
                -4.0 * ci[0] * np.power(temperature, -5.0)
                - 2.0 * ci[1] * np.power(temperature, -3.0)
                - ci[2] * np.power(temperature, -2.0)
                + ci[4]
                + 2.0 * ci[5] * temperature
            )
        return cs

    def _delta_pressure(self, volume, pressure, temperature, params):
        return self.pressure(temperature, volume, params) - pressure

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """
        if not (
            np.isfinite(pressure)
            and 0.0 < pressure <= params.get("P_max", np.inf)
            and np.isfinite(temperature)
            and temperature > 0.0
            and params.get("T_min", 0.0) <= temperature
            and temperature <= params.get("T_max", np.inf)
        ):
            raise ValueError("Pressure or temperature outside the PS1994 domain.")
        # Scan density, rather than bracketing the entire isotherm once.
        # Below the critical temperature both liquid and vapour roots may
        # be stable mechanically; the equilibrium root has the lowest G.
        ideal_density = pressure / (constants.gas_constant * temperature * 1.0e6)
        low = max(np.finfo(float).tiny, min(ideal_density * 0.01, 1.0e-6))
        densities = np.geomspace(low, 0.16, 321)
        residuals = self.pressure(temperature, 1.0e-6 / densities, params) - pressure
        volumes = []
        for i in np.flatnonzero((residuals[:-1] <= 0.0) & (residuals[1:] >= 0.0)):
            log_density = opt.brentq(
                lambda log_rho: self._delta_pressure(
                    1.0e-6 / np.exp(log_rho), pressure, temperature, params
                ),
                np.log(densities[i]),
                np.log(densities[i + 1]),
                xtol=1.0e-14,
            )
            volume = 1.0e-6 / np.exp(log_density)
            if (
                self.isothermal_bulk_modulus_reuss(
                    pressure, temperature, volume, params
                )
                > 0.0
            ):
                volumes.append(volume)
        if not volumes:
            raise ValueError("No stable PS1994 density root found.")
        return min(
            volumes,
            key=lambda v: self.gibbs_energy(pressure, temperature, v, params),
        )

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given
        temperature and volume [Pa]
        """
        # convert volume in m^3/mol to specific density in mol/cm^3
        rho_sp = 1.0e-6 / volume
        c = self._compute_cs(params["c_coeffs"], temperature)
        reduced_P = (
            rho_sp
            + c[0] * np.power(rho_sp, 2.0)
            - np.power(rho_sp, 2.0)
            * (
                (
                    c[2]
                    + 2.0 * c[3] * rho_sp
                    + 3.0 * c[4] * np.power(rho_sp, 2.0)
                    + 4.0 * c[5] * np.power(rho_sp, 3.0)
                )
                / np.power(
                    c[1]
                    + c[2] * rho_sp
                    + c[3] * np.power(rho_sp, 2.0)
                    + c[4] * np.power(rho_sp, 3.0)
                    + c[5] * np.power(rho_sp, 4.0),
                    2.0,
                )
            )
            + c[6] * np.power(rho_sp, 2.0) * np.exp(-c[7] * rho_sp)
            + c[8] * np.power(rho_sp, 2.0) * np.exp(-c[9] * rho_sp)
        )

        # Expression returns reduced pressure in MPa, so multiply by 1.e6
        return reduced_P * constants.gas_constant * temperature * 1.0e6

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]`
        """
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        return alpha * K_T * volume / C_v

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        dV = volume * 1.0e-5
        dPdV = (
            self.pressure(temperature, volume + dV / 2.0, params)
            - self.pressure(temperature, volume - dV / 2.0, params)
        ) / dV
        return -volume * dPdV

    def isentropic_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus. :math:`[Pa]`
        """
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        return K_T * C_p / C_v

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus. :math:`[Pa]`
        Gas phases have a shear modulus of 0.
        """
        return 0.0

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        dT = 0.01
        dSdT = (
            self.entropy(pressure, temperature + dT / 2.0, volume, params)
            - self.entropy(pressure, temperature - dT / 2.0, volume, params)
        ) / dT
        return temperature * dSdT

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        return C_v + volume * temperature * alpha * alpha * K_T

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        dT = 0.01
        dPdT = (
            self.pressure(temperature + dT / 2.0, volume, params)
            - self.pressure(temperature - dT / 2.0, volume, params)
        ) / dT
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        return dPdT / K_T

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the volume and temperature
        of the mineral [J/mol]
        """
        G = (
            self.helmholtz_energy(pressure, temperature, volume, params)
            + pressure * volume
        )
        return G

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy at the volume and temperature
        of the mineral [J/mol]
        """
        return self.helmholtz_energy(pressure, temperature, volume, params) + (
            temperature * self.entropy(pressure, temperature, volume, params)
        )

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the volume and temperature
        of the mineral [J/mol]
        """
        # convert volume in m^3/mol to specific density in mol/cm^3
        rho_sp = 1.0e-6 / volume
        c = self._compute_cs(params["c_coeffs"], temperature)

        g = (
            c[1]
            + c[2] * rho_sp
            + c[3] * np.power(rho_sp, 2.0)
            + c[4] * np.power(rho_sp, 3.0)
            + c[5] * np.power(rho_sp, 4.0)
        )

        reduced_F = (
            np.log(rho_sp)
            + c[0] * rho_sp
            + (1.0 / g - 1.0 / c[1])
            - (c[6] / c[7]) * np.expm1(-c[7] * rho_sp)
            - (c[8] / c[9]) * np.expm1(-c[9] * rho_sp)
        )

        dcdT = self._compute_dcdTs(params["c_coeffs"], temperature)

        dgdT = (
            dcdT[1]
            + dcdT[2] * rho_sp
            + dcdT[3] * np.power(rho_sp, 2.0)
            + dcdT[4] * np.power(rho_sp, 3.0)
            + dcdT[5] * np.power(rho_sp, 4.0)
        )
        h1 = (dcdT[6] / c[7] - c[6] * dcdT[7] / c[7] ** 2) * np.expm1(
            -c[7] * rho_sp
        ) - c[6] / c[7] * rho_sp * dcdT[7] * np.exp(-c[7] * rho_sp)
        h2 = (dcdT[8] / c[9] - c[8] * dcdT[9] / c[9] ** 2) * np.expm1(
            -c[9] * rho_sp
        ) - c[8] / c[9] * rho_sp * dcdT[9] * np.exp(-c[9] * rho_sp)

        reduced_dFdT = (
            +dcdT[0] * rho_sp - dgdT / (g * g) + dcdT[1] / (c[1] * c[1]) - h1 - h2
        )

        if "Cp" in params:
            _, ideal_entropy = self._ideal_gas_reference(temperature, params)
            S_thermal = ideal_entropy - constants.gas_constant * np.log(
                constants.gas_constant * temperature * 1.0e6 / params["P_0"]
            )
        else:
            S_thermal = params["Cv_0"] * np.log(temperature) + debye.entropy(
                temperature, params["Debye_0"], n=params["Debye_n"]
            )

        S = (
            -constants.gas_constant * (reduced_F + reduced_dFdT * temperature)
            + S_thermal
        )
        return S

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the volume and temperature
        of the mineral [J/mol]
        """

        return self.helmholtz_energy(pressure, temperature, volume, params) + (
            temperature * self.entropy(pressure, temperature, volume, params)
            + pressure * volume
        )

    def helmholtz_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the volume and temperature
        of the mineral [J/mol]
        """
        # convert volume in m^3/mol to specific density in mol/cm^3
        rho_sp = 1.0e-6 / volume
        c = self._compute_cs(params["c_coeffs"], temperature)
        reduced_F = (
            np.log(rho_sp)
            + c[0] * rho_sp
            + (
                1.0
                / (
                    c[1]
                    + c[2] * rho_sp
                    + c[3] * np.power(rho_sp, 2.0)
                    + c[4] * np.power(rho_sp, 3.0)
                    + c[5] * np.power(rho_sp, 4.0)
                )
                - 1.0 / c[1]
            )
            - (c[6] / c[7]) * np.expm1(-c[7] * rho_sp)
            - (c[8] / c[9]) * np.expm1(-c[9] * rho_sp)
        )

        if "Cp" in params:
            enthalpy, entropy = self._ideal_gas_reference(temperature, params)
            F_thermal = (
                enthalpy
                - temperature * entropy
                + constants.gas_constant
                * temperature
                * (
                    np.log(constants.gas_constant * temperature * 1.0e6 / params["P_0"])
                    - 1.0
                )
            )
        else:
            F_thermal = (
                params["Cv_0"] * (temperature - temperature * np.log(temperature))
                + debye.helmholtz_energy(
                    temperature, params["Debye_0"], n=params["Debye_n"]
                )
                + params["F_0"]
            )
        return reduced_F * constants.gas_constant * temperature + F_thermal

    def _ideal_gas_reference(self, temperature, params):
        """Ideal-gas enthalpy and entropy using the HP2011 Cp polynomial."""
        t, t0 = temperature, params["T_0"]
        a, b, c, d = params["Cp"]
        enthalpy = (
            params["H_0"]
            + a * (t - t0)
            + b * (t * t - t0 * t0) / 2.0
            - c * (1.0 / t - 1.0 / t0)
            + 2.0 * d * (np.sqrt(t) - np.sqrt(t0))
        )
        entropy = (
            params["S_0"]
            + a * np.log(t / t0)
            + b * (t - t0)
            - c * (1.0 / t**2 - 1.0 / t0**2) / 2.0
            - 2.0 * d * (1.0 / np.sqrt(t) - 1.0 / np.sqrt(t0))
        )
        return enthalpy, entropy

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = [
            "molar_mass",
            "n",
            "formula",
            "c_coeffs",
        ]
        if "Cp" in params:
            params.setdefault("T_0", 298.15)
            params.setdefault("P_0", 1.0e5)
            expected_keys += ["H_0", "S_0"]
            if len(params["Cp"]) != 4:
                raise ValueError("PS1994 Cp must have four coefficients.")
            if params["T_0"] <= 0.0 or params["P_0"] <= 0.0:
                raise ValueError(
                    "PS1994 reference temperature and pressure must be positive."
                )
        else:
            expected_keys += ["Debye_n", "Debye_0", "Cv_0", "F_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)
        if np.asarray(params["c_coeffs"]).shape != (10, 6):
            raise ValueError("PS1994 requires a 10 by 6 coefficient array.")
