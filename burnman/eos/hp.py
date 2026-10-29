# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import warnings

from . import modified_tait as mt
from . import murnaghan as murn
from . import equation_of_state as eos

from . import einstein


class HP_TMT(eos.EquationOfState):
    """
    Base class for the thermal equation of state based on
    the isothermal modified Tait equation of state (class MT),
    as described in :cite:`HP2011`.

    .. math::
        \\mathcal{G}(P,T) &= \\mathcal{H}_{\\textrm{1 bar, T}} - T\\mathcal{S}_{\\textrm{1 bar, T}} + \\int_{\\textrm{1 bar}}^P V(P,T)~dP, \\\\
        \\mathcal{H}_{\\textrm{1 bar, T}} &= \\Delta_f\\mathcal{H}_{\\textrm{1 bar, 298 K}} + \\int_{298}^T C_P~dT, \\\\
        \\mathcal{S}_{\\textrm{1 bar, T}} &= \\mathcal{S}_{\\textrm{1 bar, 298 K}} + \\int_{298}^T \\frac{C_P}{T}~dT, \\\\
        \\int_{\\textrm{1 bar}}^P V(P,T)~dP &= P V_0 \\left( 1 - a + \\left( a \\frac{(1-b P_{th})^{1-c} - (1 + b(P-P_{th}))^{1-c}}{b (c-1) P} \\right) \\right)

    The heat capacity at one bar is given by an empirical polynomial fit to
    experimental data:

    .. math::
        C_p = a + bT + cT^{-2} + dT^{-0.5}

    The thermal pressure is calculated using a Mie-Grüneisen-like function:

    .. math::
        P_{\\textrm{th}} &= \\frac{\\alpha_0 K_0 E_{\\textrm{th}} }{C_{V0}}, \\\\
        E_{\\textrm{th}} &= 3 n R \\Theta \\left(0.5 + \\frac{1}{ \\exp(\\frac{\\Theta}{T}) - 1 }\\right), \\\\
        C_{V} &= 3 n R \\frac{(\\frac{\\Theta}{T})^2\\exp(\\frac{\\Theta}{T})}{(\\exp(\\frac{\\Theta}{T})-1)^2}


    :math:`\\Theta` is the Einstein temperature of the crystal in Kelvin,
    approximated for a substance :math:`i` with :math:`n_i` atoms in the
    unit formula and a molar entropy :math:`S_i` using the empirical formula:

    .. math::
        \\Theta_i=\\frac{10636}{S_i/n_i + 6.44}

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``P_0``
          - Reference pressure.
          - Pa
        * - ``T_0``
          - Reference temperature.
          - K
        * - ``H_0``
          - Enthalpy at the reference state.
          - J/mol
        * - ``S_0``
          - Entropy at the reference state.
          - J/(mol K)
        * - ``Cp``
          - Heat capacity coefficients at the reference pressure (length 4 list).
          - Various
        * - ``V_0``
          - Volume at the reference state.
          - m^3/mol
        * - ``K_0``
          - Isothermal bulk modulus at the reference state.
          - Pa
        * - ``Kprime_0``
          - Pressure derivative of bulk modulus at the reference state.
          - Dimensionless
        * - ``Kdprime_0``
          - Second pressure derivative of bulk modulus at the reference state.
          - Dimensionless
        * - ``a_0``
          - Thermal expansivity at reference state.
          - 1/K
        * - ``T_einstein``
          - Einstein temperature.
          - K
        * - ``n``
          - Number of atoms in the unit formula.
          - Dimensionless

    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        # Equation 12 of HP2011
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.volume(pressure - Pth, params)

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        # Between Equations 11 and 12 of HP2011
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.pressure_modified_tait(volume / params["V_0"], params) + Pth

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        # Two equations after Equations 13 of HP2011
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.bulk_modulus(pressure - Pth, params)

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.0

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]. This function replaces -Pth in Equation 13+1
        in :cite:`HP2011` with P-Pth for non-ambient temperature
        """
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)
        psubpth = pressure - params["P_0"] - Pth

        C_V0 = einstein.molar_heat_capacity_v(
            params["T_0"], params["T_einstein"], params["n"]
        )
        C_V = einstein.molar_heat_capacity_v(
            temperature, params["T_einstein"], params["n"]
        )
        alpha = (
            params["a_0"]
            * (C_V / C_V0)
            * 1.0
            / ((1.0 + b * psubpth) * (a + (1.0 - a) * np.power((1 + b * psubpth), c)))
        )
        return alpha

    def _molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]. Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP2011`.
        """
        Cp = (
            params["Cp"][0]
            + params["Cp"][1] * temperature
            + params["Cp"][2] * np.power(temperature, -2.0)
            + params["Cp"][3] * np.power(temperature, -0.5)
        )
        return Cp

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # Calculate temperature and pressure integrals
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)

        psubpth = pressure - params["P_0"] - Pth

        # Equation 13 in HP2011
        if pressure != params["P_0"]:
            intVdP = (
                (pressure - params["P_0"])
                * params["V_0"]
                * (
                    1.0
                    - a
                    + (
                        a
                        * (
                            np.power((1.0 - b * Pth), 1.0 - c)
                            - np.power((1.0 + b * (psubpth)), 1.0 - c)
                        )
                        / (b * (c - 1.0) * (pressure - params["P_0"]))
                    )
                )
            )
        else:
            intVdP = 0.0
        return (
            params["H_0"]
            + self.__intCpdT(temperature, params)
            - temperature * (params["S_0"] + self.__intCpoverTdT(temperature, params))
            + intVdP
        )

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)

        ksi_over_ksi_0 = einstein.molar_heat_capacity_v(
            temperature, params["T_einstein"], params["n"]
        ) / einstein.molar_heat_capacity_v(
            params["T_0"], params["T_einstein"], params["n"]
        )

        dintVdpdT = (
            params["V_0"] * params["a_0"] * params["K_0"] * a * ksi_over_ksi_0
        ) * (
            np.power((1.0 + b * (pressure - params["P_0"] - Pth)), 0.0 - c)
            - np.power((1.0 - b * Pth), 0.0 - c)
        )
        return params["S_0"] + self.__intCpoverTdT(temperature, params) + dintVdpdT

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = mt.tait_constants(params)
        T = temperature
        T_e = params["T_einstein"]
        n = params["n"]
        Pth = self.__relative_thermal_pressure(T, params)

        ksi_over_ksi_0 = einstein.molar_heat_capacity_v(
            T, T_e, n
        ) / einstein.molar_heat_capacity_v(params["T_0"], T_e, n)

        dintVdpdT = (
            params["V_0"] * params["a_0"] * params["K_0"] * a * ksi_over_ksi_0
        ) * (
            np.power((1.0 + b * (pressure - params["P_0"] - Pth)), 0.0 - c)
            - np.power((1.0 - b * Pth), 0.0 - c)
        )

        dSdT0 = (
            params["V_0"]
            * params["K_0"]
            * np.power((ksi_over_ksi_0 * params["a_0"]), 2.0)
            * (
                np.power((1.0 + b * (pressure - params["P_0"] - Pth)), -1.0 - c)
                - np.power((1.0 + b * (-Pth)), -1.0 - c)
            )
        )

        x = T_e / T
        dCv_einstdT = -(
            einstein.molar_heat_capacity_v(T, T_e, n)
            * (1 - 2.0 / x + 2.0 / (np.exp(x) - 1.0))
            * x
            / T
        )

        dSdT1 = -dintVdpdT * dCv_einstdT / einstein.molar_heat_capacity_v(T, T_e, n)

        dSdT = dSdT0 + dSdT1
        return self._molar_heat_capacity_p0(temperature, params) + temperature * dSdT

    def __thermal_pressure(self, T, params):
        """
        Returns thermal pressure [Pa] as a function of T [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR. This function is *not* used to calculate the
        # heat capacity - Holland and Powell (2011) prefer the additional
        # freedom provided by their polynomial expression.

        E_th = einstein.thermal_energy(T, params["T_einstein"], params["n"])
        C_V0 = einstein.molar_heat_capacity_v(
            params["T_0"], params["T_einstein"], params["n"]
        )
        P_th = params["a_0"] * params["K_0"] / C_V0 * E_th
        return P_th

    def __relative_thermal_pressure(self, T, params):
        """
        Returns relative thermal pressure [Pa] as a function of T-params['T_0'] [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """
        return self.__thermal_pressure(T, params) - self.__thermal_pressure(
            params["T_0"], params
        )

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * temperature
            + 0.5 * params["Cp"][1] * np.power(temperature, 2.0)
            - params["Cp"][2] / temperature
            + 2.0 * params["Cp"][3] * np.sqrt(temperature)
        ) - (
            params["Cp"][0] * params["T_0"]
            + 0.5 * params["Cp"][1] * params["T_0"] * params["T_0"]
            - params["Cp"][2] / params["T_0"]
            + 2.0 * params["Cp"][3] * np.sqrt(params["T_0"])
        )

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * np.log(temperature)
            + params["Cp"][1] * temperature
            - 0.5 * params["Cp"][2] / np.power(temperature, 2.0)
            - 2.0 * params["Cp"][3] / np.sqrt(temperature)
        ) - (
            params["Cp"][0] * np.log(params["T_0"])
            + params["Cp"][1] * params["T_0"]
            - 0.5 * params["Cp"][2] / (params["T_0"] * params["T_0"])
            - 2.0 * params["Cp"][3] / np.sqrt(params["T_0"])
        )

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if "T_0" not in params:
            params["T_0"] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if "H_0" not in params:
            params["H_0"] = float("nan")
        if "S_0" not in params:
            params["S_0"] = float("nan")

        # First, let's check the EoS parameters for Tref
        mt.MT.validate_parameters(mt.MT(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["H_0", "S_0", "V_0", "Cp", "a_0", "n", "molar_mass"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # The following line estimates the Einstein temperature
        # according to the empirical equation of
        # Holland and Powell, 2011; base of p.346, para.1
        if "T_einstein" not in params:
            params["T_einstein"] = 10636.0 / (params["S_0"] / params["n"] + 6.44)

        # Finally, check that the values are reasonable.
        if params["T_0"] < 0.0:
            warnings.warn("Unusual value for T_0", stacklevel=2)
        if params["G_0"] is not float("nan") and (
            params["G_0"] < 0.0 or params["G_0"] > 1.0e13
        ):
            warnings.warn("Unusual value for G_0", stacklevel=2)
        if params["Gprime_0"] is not float("nan") and (
            params["Gprime_0"] < -5.0 or params["Gprime_0"] > 10.0
        ):
            warnings.warn("Unusual value for Gprime_0", stacklevel=2)

        # no test for H_0
        if params["S_0"] is not float("nan") and params["S_0"] < 0.0:
            warnings.warn("Unusual value for S_0", stacklevel=2)
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-2:
            warnings.warn("Unusual value for V_0", stacklevel=2)

        if self._molar_heat_capacity_p0(params["T_0"], params) < 0.0:
            warnings.warn("Negative heat capacity at T_0", stacklevel=2)
        if self._molar_heat_capacity_p0(2000.0, params) < 0.0:
            warnings.warn("Negative heat capacity at 2000K", stacklevel=2)

        if params["a_0"] < 0.0 or params["a_0"] > 1.0e-3:
            warnings.warn("Unusual value for a_0", stacklevel=2)

        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
        if params["molar_mass"] < 0.001 or params["molar_mass"] > 10.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)


class HP_TMTL(eos.EquationOfState):
    """
    Base class for the thermal equation of state
    described in :cite:`HP1998`, but with the Modified Tait as the static part,
    as described in :cite:`HP2011`.

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``P_0``
          - Reference pressure.
          - Pa
        * - ``T_0``
          - Reference temperature.
          - K
        * - ``H_0``
          - Enthalpy at the reference state.
          - J/mol
        * - ``S_0``
          - Entropy at the reference state.
          - J/(mol K)
        * - ``Cp``
          - Heat capacity coefficients at the reference pressure (length 4 list).
          - Various
        * - ``V_0``
          - Volume at the reference state.
          - m^3/mol
        * - ``K_0``
          - Isothermal bulk modulus at the reference state.
          - Pa
        * - ``Kprime_0``
          - Pressure derivative of bulk modulus at the reference state.
          - Dimensionless
        * - ``Kdprime_0``
          - Second pressure derivative of bulk modulus at the reference state.
          - Dimensionless
        * - ``a_0``
          - Thermal expansivity at reference state.
          - 1/K
        * - ``dKdT_0``
          - Temperature derivative of bulk modulus at reference state.
          - Pa/K
    """

    def _V_T_1bar(self, temperature, params):
        # Constant thermal expansivity at standard state pressure
        # (p.348 of HP2011)
        return params["V_0"] * np.exp(params["a_0"] * (temperature - params["T_0"]))

    def _K_T_1bar(self, temperature, params):
        # Linear bulk modulus dependence as in HP1998 (p.348 of HP2011)
        return params["K_0"] + params["dKdT_0"] * (temperature - params["T_0"])

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        self.static_params["V_0"] = self._V_T_1bar(temperature, params)
        self.static_params["K_0"] = self._K_T_1bar(temperature, params)
        return mt.volume(pressure, self.static_params)

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        self.static_params["V_0"] = self._V_T_1bar(temperature, params)
        self.static_params["K_0"] = self._K_T_1bar(temperature, params)
        return mt.bulk_modulus(pressure, self.static_params)

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.0

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]
        """
        # The derivation of the high pressure thermal expansivity is tedious,
        # so here we take a numerical derivative.
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        self.static_params["V_0"] = self._V_T_1bar(temperature + dT / 2.0, params)
        self.static_params["K_0"] = self._K_T_1bar(temperature + dT / 2.0, params)
        volume1 = mt.volume(pressure, self.static_params)
        self.static_params["V_0"] = self._V_T_1bar(temperature - dT / 2.0, params)
        self.static_params["K_0"] = self._K_T_1bar(temperature - dT / 2.0, params)
        volume0 = mt.volume(pressure, self.static_params)

        return 2.0 * (volume1 - volume0) / (volume1 + volume0) / dT

    def _molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP1998`.
        """
        Cp = (
            params["Cp"][0]
            + params["Cp"][1] * temperature
            + params["Cp"][2] * np.power(temperature, -2.0)
            + params["Cp"][3] * np.power(temperature, -0.5)
        )
        return Cp

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        self.static_params["V_0"] = self._V_T_1bar(temperature, params)
        self.static_params["K_0"] = self._K_T_1bar(temperature, params)
        return (
            params["H_0"]
            + self.__intCpdT(temperature, params)
            - temperature * (params["S_0"] + self.__intCpoverTdT(temperature, params))
            + mt.intVdP(pressure, self.static_params)
        )

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The derivation of the entropy is tedious,
        # so here we take a numerical derivative.
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        G1 = self.gibbs_energy(pressure, temperature + dT / 2.0, volume, params)
        G0 = self.gibbs_energy(pressure, temperature - dT / 2.0, volume, params)

        return (G0 - G1) / dT

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The differentiation is tedious, so for now we just take the
        # numerical derivative of S
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        S1 = self.entropy(pressure, temperature + dT / 2.0, volume, params)
        S0 = self.entropy(pressure, temperature - dT / 2.0, volume, params)
        return temperature * (S1 - S0) / dT

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * temperature
            + 0.5 * params["Cp"][1] * np.power(temperature, 2.0)
            - params["Cp"][2] / temperature
            + 2.0 * params["Cp"][3] * np.sqrt(temperature)
        ) - (
            params["Cp"][0] * params["T_0"]
            + 0.5 * params["Cp"][1] * params["T_0"] * params["T_0"]
            - params["Cp"][2] / params["T_0"]
            + 2.0 * params["Cp"][3] * np.sqrt(params["T_0"])
        )

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * np.log(temperature)
            + params["Cp"][1] * temperature
            - 0.5 * params["Cp"][2] / np.power(temperature, 2.0)
            - 2.0 * params["Cp"][3] / np.sqrt(temperature)
        ) - (
            params["Cp"][0] * np.log(params["T_0"])
            + params["Cp"][1] * params["T_0"]
            - 0.5 * params["Cp"][2] / (params["T_0"] * params["T_0"])
            - 2.0 * params["Cp"][3] / np.sqrt(params["T_0"])
        )

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if "T_0" not in params:
            params["T_0"] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if "H_0" not in params:
            params["H_0"] = float("nan")
        if "S_0" not in params:
            params["S_0"] = float("nan")

        # First, let's check the EoS parameters for Tref
        mt.MT.validate_parameters(mt.MT(), params)
        self.static_params = {
            "V_0": params["V_0"],
            "K_0": params["K_0"],
            "Kprime_0": params["Kprime_0"],
            "Kdprime_0": params["Kdprime_0"],
            "P_0": params["P_0"],
        }

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["H_0", "S_0", "V_0", "Cp", "a_0", "dKdT_0", "n", "molar_mass"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["T_0"] < 0.0:
            warnings.warn("Unusual value for T_0", stacklevel=2)
        if params["G_0"] is not float("nan") and (
            params["G_0"] < 0.0 or params["G_0"] > 1.0e13
        ):
            warnings.warn("Unusual value for G_0", stacklevel=2)
        if params["Gprime_0"] is not float("nan") and (
            params["Gprime_0"] < -5.0 or params["Gprime_0"] > 10.0
        ):
            warnings.warn("Unusual value for Gprime_0", stacklevel=2)

        # no test for H_0 or S_0 (several HP endmembers have S_0 < 0)
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-2:
            warnings.warn("Unusual value for V_0", stacklevel=2)

        if self._molar_heat_capacity_p0(params["T_0"], params) < 0.0:
            warnings.warn("Negative heat capacity at T_0", stacklevel=2)
        if self._molar_heat_capacity_p0(2000.0, params) < 0.0:
            warnings.warn("Negative heat capacity at 2000K", stacklevel=2)

        if params["a_0"] < 0.0 or params["a_0"] > 1.0e-3:
            warnings.warn("Unusual value for a_0", stacklevel=2)

        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
        if params["molar_mass"] < 0.001 or params["molar_mass"] > 10.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)


class HP98(eos.EquationOfState):
    """
    Base class for the thermal equation of state
    described in :cite:`HP1998`.

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``P_0``
          - Reference pressure.
          - Pa
        * - ``T_0``
          - Reference temperature.
          - K
        * - ``H_0``
          - Enthalpy at the reference state.
          - J/mol
        * - ``S_0``
          - Entropy at the reference state.
          - J/(mol K)
        * - ``Cp``
          - Heat capacity coefficients at the reference pressure (length 4 list).
          - Various
        * - ``V_0``
          - Volume at the reference state.
          - m^3/mol
        * - ``K_0``
          - Isothermal bulk modulus at the reference state.
          - Pa
        * - ``Kprime_0``
          - Pressure derivative of bulk modulus at the reference state.
          - Dimensionless
        * - ``a_0``
          - Thermal expansivity at reference state.
          - 1/K
        * - ``dKdT_0``
          - Temperature derivative of bulk modulus at reference state.
          - Pa/K
    """

    def _V_T_1bar(self, temperature, params):
        return params["V_0"] * (
            1.0
            + params["a_0"] * (temperature - params["T_0"])
            - 20.0 * params["a_0"] * (np.sqrt(temperature) - np.sqrt(params["T_0"]))
        )

    def _K_T_1bar(self, temperature, params):
        return params["K_0"] + params["dKdT_0"] * (temperature - params["T_0"])

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        return murn.volume(
            pressure,
            self._V_T_1bar(temperature, params),
            self._K_T_1bar(temperature, params),
            params["Kprime_0"],
        )

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        return murn.bulk_modulus(
            pressure, self._K_T_1bar(temperature, params), params["Kprime_0"]
        )

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.0

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]
        """
        VT = self._V_T_1bar(temperature, params)
        KT = self._K_T_1bar(temperature, params)
        volume = murn.volume(pressure, VT, KT, params["Kprime_0"])
        dVTdT = params["V_0"] * params["a_0"] * (1.0 - 10.0 / np.sqrt(temperature))
        g = volume / VT
        dgdKT = (
            pressure
            * np.power(
                1.0 + pressure * params["Kprime_0"] / KT, -1 - 1.0 / params["Kprime_0"]
            )
            / (KT * KT)
        )
        dVdT = dVTdT * g + VT * dgdKT * params["dKdT_0"]
        return dVdT / volume

    def _molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP1998`.
        """
        Cp = (
            params["Cp"][0]
            + params["Cp"][1] * temperature
            + params["Cp"][2] * np.power(temperature, -2.0)
            + params["Cp"][3] * np.power(temperature, -0.5)
        )
        return Cp

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        return (
            params["H_0"]
            + self.__intCpdT(temperature, params)
            - temperature * (params["S_0"] + self.__intCpoverTdT(temperature, params))
            + murn.intVdP(
                pressure,
                self._V_T_1bar(temperature, params),
                self._K_T_1bar(temperature, params),
                params["Kprime_0"],
            )
        )

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The entropy involves differentiating intdVdp
        # with respect to temperature.
        # We do this using the product and chain rules

        VT = self._V_T_1bar(temperature, params)
        KT = self._K_T_1bar(temperature, params)
        dVTdT = params["V_0"] * params["a_0"] * (1.0 - 10.0 / np.sqrt(temperature))
        g = murn.intVdP(pressure, VT, KT, params["Kprime_0"]) / VT
        dgdKT = (
            (pressure / KT + 1.0)
            * np.power(
                1.0 + pressure * params["Kprime_0"] / KT, -1.0 / params["Kprime_0"]
            )
            - 1.0
        ) / (params["Kprime_0"] - 1.0)
        dintVdpdT = dVTdT * g + VT * dgdKT * params["dKdT_0"]
        return params["S_0"] + self.__intCpoverTdT(temperature, params) - dintVdpdT

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The differentiation is tedious, so for now we just take the
        # numerical derivative of S
        # TODO calculate the analytical derivative
        dT = 0.1
        S1 = self.entropy(pressure, temperature + dT / 2.0, volume, params)
        S0 = self.entropy(pressure, temperature - dT / 2.0, volume, params)
        return temperature * (S1 - S0) / dT

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * temperature
            + 0.5 * params["Cp"][1] * np.power(temperature, 2.0)
            - params["Cp"][2] / temperature
            + 2.0 * params["Cp"][3] * np.sqrt(temperature)
        ) - (
            params["Cp"][0] * params["T_0"]
            + 0.5 * params["Cp"][1] * params["T_0"] * params["T_0"]
            - params["Cp"][2] / params["T_0"]
            + 2.0 * params["Cp"][3] * np.sqrt(params["T_0"])
        )

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return (
            params["Cp"][0] * np.log(temperature)
            + params["Cp"][1] * temperature
            - 0.5 * params["Cp"][2] / np.power(temperature, 2.0)
            - 2.0 * params["Cp"][3] / np.sqrt(temperature)
        ) - (
            params["Cp"][0] * np.log(params["T_0"])
            + params["Cp"][1] * params["T_0"]
            - 0.5 * params["Cp"][2] / (params["T_0"] * params["T_0"])
            - 2.0 * params["Cp"][3] / np.sqrt(params["T_0"])
        )

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if "T_0" not in params:
            params["T_0"] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if "H_0" not in params:
            params["H_0"] = float("nan")
        if "S_0" not in params:
            params["S_0"] = float("nan")

        # First, let's check the EoS parameters for Tref
        murn.Murnaghan.validate_parameters(murn.Murnaghan(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["H_0", "S_0", "V_0", "Cp", "a_0", "dKdT_0", "n", "molar_mass"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["T_0"] < 0.0:
            warnings.warn("Unusual value for T_0", stacklevel=2)
        if params["G_0"] is not float("nan") and (
            params["G_0"] < 0.0 or params["G_0"] > 1.0e13
        ):
            warnings.warn("Unusual value for G_0", stacklevel=2)
        if params["Gprime_0"] is not float("nan") and (
            params["Gprime_0"] < -5.0 or params["Gprime_0"] > 10.0
        ):
            warnings.warn("Unusual value for Gprime_0", stacklevel=2)

        # no test for H_0 or S_0 (several HP endmembers have S_0 < 0)
        if params["V_0"] < 1.0e-7 or params["V_0"] > 1.0e-2:
            warnings.warn("Unusual value for V_0", stacklevel=2)

        if self._molar_heat_capacity_p0(params["T_0"], params) < 0.0:
            warnings.warn("Negative heat capacity at T_0", stacklevel=2)
        if self._molar_heat_capacity_p0(2000.0, params) < 0.0:
            warnings.warn("Negative heat capacity at 2000K", stacklevel=2)

        if params["a_0"] < 0.0 or params["a_0"] > 1.0e-3:
            warnings.warn("Unusual value for a_0", stacklevel=2)

        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
        if params["molar_mass"] < 0.001 or params["molar_mass"] > 10.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)
