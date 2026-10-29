# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2023 by the BurnMan team, released under the GNU
# GPL v2 or later.

# TO DO: Correct heat capacity, volume where internal order-disorder is
# implemented (Landau and Bragg-Williams models)


import numpy as np

from . import equation_of_state as eos
from .. import constants

import warnings


def _cork_variables(cork, cork_P, cork_T, temperature):
    a = (
        cork[0][0] * cork_T ** (2.5) / cork_P
        + cork[0][1] * cork_T ** (1.5) / cork_P * temperature
    )
    b = cork[1][0] * cork_T / cork_P
    c = (
        cork[2][0] * cork_T / cork_P ** (1.5)
        + cork[2][1] / cork_P ** (1.5) * temperature
    )
    d = (
        cork[3][0] * cork_T / cork_P ** (2.0)
        + cork[3][1] / cork_P ** (2.0) * temperature
    )
    return [a, b, c, d]


class CORK(eos.EquationOfState):
    """
    Class for the CoRK equation of state detailed in :cite:`HP1991`. The
    CoRK EoS is a simple virial-type extension to the modified Redlich-Kwong
    (MRK) equation of state. It was designed to compensate for the tendency of
    the MRK equation of state to overestimate volumes at high pressures and
    accommodate the volume behaviour of coexisting gas and liquid phases along
    the saturation curve.

    .. math::
        V &= \\frac{RT}{P} + c_1 - \\frac{c_0 R T^{0.5}}{(RT + c_1 P)(RT + 2 c_1 P)} + c_2 P^{0.5} + c_3 P, \\\\
        c_0 &= c_{0,0} T_c^{2.5}/P_c + c_{0,1} T_c^{1.5}/P_c T, \\\\
        c_1 &= c_{1,0} T_c/P_c, \\\\
        c_2 &= c_{2,0} T_c/P_c^{1.5} + c_{2,1}/P_c^{1.5} T, \\\\
        c_3 &= c_{3,0} T_c/P_c^2 + c_{3,1}/P_c^2 T

    where :math:`c_{i,j}` are the CoRK parameters,
    :math:`T_c` is the critical temperature,
    and :math:`P_c` is the critical pressure.

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``H_0``
          - Reference enthalpy.
          - :math:`\\text{J/mol}`
        * - ``S_0``
          - Reference entropy.
          - :math:`\\text{J/(mol K)}`
        * - ``Cp``
          - Heat capacity parameters (length 4 array).
          - :math:`\\text{J/(mol K)}`
        * - ``T_0``
          - Reference temperature.
          - :math:`\\text{K}`
        * - ``P_0``
          - Reference pressure.
          - :math:`\\text{Pa}`
        * - ``cork_params``
          - CoRK parameters (length 4x2 array).
          - varies
        * - ``cork_T``
          - Critical temperature.
          - :math:`\\text{K}`
        * - ``cork_P``
          - Critical pressure.
          - :math:`\\text{Pa}`

    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        Eq. 7 in Holland and Powell, 1991
        """
        cork = _cork_variables(
            params["cork_params"], params["cork_P"], params["cork_T"], temperature
        )
        V = constants.gas_constant * temperature / pressure + (
            cork[1]
            - cork[0]
            * constants.gas_constant
            * np.sqrt(temperature)
            / (
                (constants.gas_constant * temperature + cork[1] * pressure)
                * (constants.gas_constant * temperature + 2.0 * cork[1] * pressure)
            )
            + cork[2] * np.sqrt(pressure)
            + cork[3] * pressure
        )
        return V

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        cork = _cork_variables(
            params["cork_params"], params["cork_P"], params["cork_T"], temperature
        )
        dVdP = -constants.gas_constant * temperature / (pressure * pressure) + (
            -cork[0]
            * cork[1]
            / np.sqrt(temperature)
            * (
                1.0
                / np.power(
                    constants.gas_constant * temperature + cork[1] * pressure, 2.0
                )
                - 4.0
                / np.power(
                    constants.gas_constant * temperature + 2.0 * cork[1] * pressure, 2.0
                )
            )
            + 0.5 * cork[2] / np.sqrt(pressure)
            + cork[3]
        )
        return -volume / dVdP

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
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        Replace -Pth in EQ 13+1 with P-Pth for non-ambient temperature
        """
        cork = _cork_variables(
            params["cork_params"], params["cork_P"], params["cork_T"], temperature
        )
        dcorkdT = params["dcorkdT"]
        RT = constants.gas_constant * temperature
        c1P = cork[1] * pressure
        dVdT = constants.gas_constant / pressure + (
            -dcorkdT[0]
            * constants.gas_constant
            * np.sqrt(temperature)
            / ((RT + c1P) * (RT + 2.0 * c1P))
            + 0.5
            * cork[0]
            * constants.gas_constant
            * (-2.0 * np.power(c1P, 2.0) + 3.0 * RT * (c1P + RT))
            / np.sqrt(temperature)
            / (np.power(RT + c1P, 2.0) * np.power(RT + 2.0 * c1P, 2.0))
            + dcorkdT[2] * np.sqrt(pressure)
            + dcorkdT[3] * pressure
        )
        return dVdT / volume

    # Heat capacity at ambient pressure
    def _molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in Holland and Powell, 2011
        """
        Cp = (
            params["Cp"][0]
            + params["Cp"][1] * temperature
            + params["Cp"][2] * np.power(temperature, -2.0)
            + params["Cp"][3] * np.power(temperature, -0.5)
        )
        return Cp

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        P_relative = pressure - params["P_0"]

        Cp0 = self._molar_heat_capacity_p0(temperature, params)

        if params["cork_T"] == 0:
            d2RTlnfdTdT = 0.0
        else:
            cork = _cork_variables(
                params["cork_params"], params["cork_P"], params["cork_T"], temperature
            )
            dcorkdT = params["dcorkdT"]

            RT = constants.gas_constant * temperature
            c1P = cork[1] * P_relative
            logterm = np.log(RT + c1P) - np.log(RT + 2.0 * c1P)
            dlogtermdT = constants.gas_constant * (
                (1.0 / (RT + c1P) - 1.0 / (RT + 2.0 * c1P))
            )
            d2logtermdTdT = -(
                constants.gas_constant
                * constants.gas_constant
                * (
                    (
                        1.0 / np.power(RT + c1P, 2.0)
                        - 1.0 / np.power(RT + 2.0 * c1P, 2.0)
                    )
                )
            )

            a = dcorkdT[0] / (cork[1] * np.sqrt(temperature)) - cork[0] / (
                2.0 * cork[1] * temperature * np.sqrt(temperature)
            )
            dadT = (
                -0.5 * dcorkdT[0] / (cork[1] * temperature * np.sqrt(temperature))
                - dcorkdT[0] / (2.0 * cork[1] * temperature * np.sqrt(temperature))
                + 3.0
                * cork[0]
                / (4.0 * cork[1] * temperature * temperature * np.sqrt(temperature))
            )
            b = cork[0] / (cork[1] * np.sqrt(temperature))
            dbdT = dcorkdT[0] / (cork[1] * np.sqrt(temperature)) - 0.5 * b / temperature
            d2RTlnfdTdT = (
                dadT * logterm + a * dlogtermdT + dbdT * dlogtermdT + b * d2logtermdTdT
            )  # Second temperature derivative of Eq. 8 in Holland and Powell, 1991

        return Cp0 - temperature * d2RTlnfdTdT

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        T_0 = params["T_0"]
        P_relative = pressure - params["P_0"]

        # Calculate temperature and pressure integrals
        intCpdT = (
            params["Cp"][0] * temperature
            + 0.5 * params["Cp"][1] * np.power(temperature, 2.0)
            - params["Cp"][2] / temperature
            + 2.0 * params["Cp"][3] * np.sqrt(temperature)
        ) - (
            params["Cp"][0] * T_0
            + 0.5 * params["Cp"][1] * T_0 * T_0
            - params["Cp"][2] / T_0
            + 2.0 * params["Cp"][3] * np.sqrt(T_0)
        )

        intCpoverTdT = (
            params["Cp"][0] * np.log(temperature)
            + params["Cp"][1] * temperature
            - 0.5 * params["Cp"][2] / np.power(temperature, 2.0)
            - 2.0 * params["Cp"][3] / np.sqrt(temperature)
        ) - (
            params["Cp"][0] * np.log(T_0)
            + params["Cp"][1] * T_0
            - 0.5 * params["Cp"][2] / (T_0 * T_0)
            - 2.0 * params["Cp"][3] / np.sqrt(T_0)
        )

        if params["cork_T"] == 0:
            RTlnf = 0.0
        else:
            cork = _cork_variables(
                params["cork_params"], params["cork_P"], params["cork_T"], temperature
            )

            RTlnf = (
                constants.gas_constant * temperature * np.log(1e-5 * P_relative)
                + cork[1] * P_relative
                + cork[0]
                / (cork[1] * np.sqrt(temperature))
                * (
                    np.log(constants.gas_constant * temperature + cork[1] * P_relative)
                    - np.log(
                        constants.gas_constant * temperature
                        + 2.0 * cork[1] * P_relative
                    )
                )
                + 2.0 / 3.0 * cork[2] * P_relative * np.sqrt(P_relative)
                + cork[3] / 2.0 * P_relative * P_relative
            )  # Eq. 8 in Holland and Powell, 1991

        return (
            params["H_0"]
            + intCpdT
            - temperature * (params["S_0"] + intCpoverTdT)
            + RTlnf
        )

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        T_0 = params["T_0"]
        P_relative = pressure - params["P_0"]

        intCpoverTdT = (
            params["Cp"][0] * np.log(temperature)
            + params["Cp"][1] * temperature
            - 0.5 * params["Cp"][2] / np.power(temperature, 2.0)
            - 2.0 * params["Cp"][3] / np.sqrt(temperature)
        ) - (
            params["Cp"][0] * np.log(T_0)
            + params["Cp"][1] * T_0
            - 0.5 * params["Cp"][2] / (T_0 * T_0)
            - 2.0 * params["Cp"][3] / np.sqrt(T_0)
        )

        if params["cork_T"] == 0:
            dRTlnfdT = 0.0
        else:
            cork = _cork_variables(
                params["cork_params"], params["cork_P"], params["cork_T"], temperature
            )
            dcorkdT = params["dcorkdT"]

            logterm = np.log(
                constants.gas_constant * temperature + cork[1] * P_relative
            ) - np.log(
                constants.gas_constant * temperature + 2.0 * cork[1] * P_relative
            )
            dRTlnfdT = (
                constants.gas_constant * np.log(1e-5 * P_relative)
                + (
                    dcorkdT[0] / (cork[1] * np.sqrt(temperature))
                    - cork[0] / (2.0 * cork[1] * temperature * np.sqrt(temperature))
                )
                * logterm
                + (cork[0] / (cork[1] * np.sqrt(temperature)))
                * constants.gas_constant
                * (
                    (
                        1.0
                        / (constants.gas_constant * temperature + cork[1] * P_relative)
                        - 1.0
                        / (
                            constants.gas_constant * temperature
                            + 2.0 * cork[1] * P_relative
                        )
                    )
                )
                + 2.0 / 3.0 * dcorkdT[2] * P_relative * np.sqrt(P_relative)
                + dcorkdT[3] / 2.0 * P_relative * P_relative
            )  # Temperature derivative of Eq. 8 in Holland and Powell, 1991

        return params["S_0"] + intCpoverTdT - dRTlnfdT

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        if "T_0" not in params:
            params["T_0"] = 298.15
        if "P_0" not in params:
            params["P_0"] = 0.0

        # if G and Gprime are not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if "H_0" not in params:
            params["H_0"] = float("nan")
        if "S_0" not in params:
            params["S_0"] = float("nan")

        cork, cork_P, cork_T = params["cork_params"], params["cork_P"], params["cork_T"]
        if "dcorkdT" not in params:
            params["dcorkdT"] = [
                cork[0][1] * cork_T ** (1.5) / cork_P,
                0.0,
                cork[2][1] / cork_P ** (1.5),
                cork[3][1] / cork_P ** (2.0),
            ]

        # check that all the required keys are in the dictionary
        expected_keys = ["cork_params", "cork_T", "cork_P", "Cp"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # now check that the values are reasonable.  I mostly just
        # made up these values from experience, and we are only
        # raising a warning.  Better way to do this? [IR]

        # no test for H_0
        if params["S_0"] is not float("nan") and params["S_0"] < 0.0:
            warnings.warn("Unusual value for S_0", stacklevel=2)

        if params["cork_T"] < -1.0:
            warnings.warn("Unusual value for cork_T", stacklevel=2)
        if params["cork_P"] < 1.0e4 or params["cork_P"] > 1.0e8:
            warnings.warn("Unusual value for cork_P", stacklevel=2)

        if self._molar_heat_capacity_p0(params["T_0"], params) < 0.0:
            warnings.warn("Negative heat capacity at T_0", stacklevel=2)
        if self._molar_heat_capacity_p0(2000.0, params) < 0.0:
            warnings.warn("Negative heat capacity at 2000K", stacklevel=2)
