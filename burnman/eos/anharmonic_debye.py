# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


class AnharmonicDebye:
    """
    Class providing methods to compute the anharmonic contribution to the
    Helmholtz free energy, pressure, entropy, isochoric heat capacity,
    isothermal bulk modulus and thermal expansion coefficient
    multiplied by the isothermal bulk modulus.

    The Helmholtz energy is defined as
    :math:`F(V, T) = A(V) (F_a(T, \\Theta(V)) - F_a(T_0, \\Theta(V)))`
    where :math:`A(V)` is the anharmonic prefactor,
    :math:`F_a(T, \\Theta)` is the anharmonic Helmholtz energy,
    and :math:`\\Theta(V)` is the Debye temperature.

    These three functions and their derivatives with respect to
    their arguments are provided as class instances contained within
    the params dictionary, with keys "anharmonic_prefactor_model",
    "debye_temperature_model", and "anharmonic_thermal_model".

    :return: _description_
    :rtype: _type_
    """

    @staticmethod
    def helmholtz_energy(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)
        F_a = anharmonic_model.nondimensional_helmholtz_energy(
            temperature, debye_T, params
        )
        F_a0 = anharmonic_model.nondimensional_helmholtz_energy(
            params["T_0"], debye_T, params
        )
        return A * (F_a - F_a0)

    @staticmethod
    def entropy(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)
        S_a = anharmonic_model.nondimensional_entropy(temperature, debye_T, params)
        return A * S_a

    @staticmethod
    def heat_capacity_v(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)
        Cv_a = anharmonic_model.nondimensional_heat_capacity(
            temperature, debye_T, params
        )
        return A * Cv_a

    @staticmethod
    def pressure(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        dAdV = params["anharmonic_prefactor_model"].dVrel(x, params) / params["V_0"]
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)
        F_a = anharmonic_model.nondimensional_helmholtz_energy(
            temperature, debye_T, params
        )
        F_a0 = anharmonic_model.nondimensional_helmholtz_energy(
            params["T_0"], debye_T, params
        )
        F_ad = anharmonic_model.nondimensional_dhelmholtz_dTheta(
            temperature, debye_T, params
        )
        F_ad0 = anharmonic_model.nondimensional_dhelmholtz_dTheta(
            params["T_0"], debye_T, params
        )
        return -(
            dAdV * (F_a - F_a0)
            + A * (theta_model.dVrel(x, params) / params["V_0"]) * (F_ad - F_ad0)
        )

    @staticmethod
    def isothermal_bulk_modulus(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        dAdV = params["anharmonic_prefactor_model"].dVrel(x, params) / params["V_0"]
        d2AdV2 = (
            params["anharmonic_prefactor_model"].dVrel2(x, params) / params["V_0"] ** 2
        )
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)

        F_a = anharmonic_model.nondimensional_helmholtz_energy(
            temperature, debye_T, params
        )
        F_a0 = anharmonic_model.nondimensional_helmholtz_energy(
            params["T_0"], debye_T, params
        )
        F_ad = anharmonic_model.nondimensional_dhelmholtz_dTheta(
            temperature, debye_T, params
        )
        F_ad0 = anharmonic_model.nondimensional_dhelmholtz_dTheta(
            params["T_0"], debye_T, params
        )
        F_add = anharmonic_model.nondimensional_d2helmholtz_dTheta2(
            temperature, debye_T, params
        )
        F_add0 = anharmonic_model.nondimensional_d2helmholtz_dTheta2(
            params["T_0"], debye_T, params
        )

        return volume * (
            d2AdV2 * (F_a - F_a0)
            + 2 * dAdV * (F_ad - F_ad0) * theta_model.dVrel(x, params) / params["V_0"]
            + A * (F_add - F_add0) * (theta_model.dVrel(x, params) / params["V_0"]) ** 2
            + A * (F_ad - F_ad0) * theta_model.dVrel2(x, params) / params["V_0"] ** 2
        )

    @staticmethod
    def dSdV(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["anharmonic_prefactor_model"].value(x, params)
        dAdV = params["anharmonic_prefactor_model"].dVrel(x, params) / params["V_0"]
        theta_model = params["debye_temperature_model"]
        anharmonic_model = params["anharmonic_thermal_model"]
        debye_T = theta_model.value(x, params)

        S_a = anharmonic_model.nondimensional_entropy(temperature, debye_T, params)
        S_ad = anharmonic_model.nondimensional_dentropy_dTheta(
            temperature, debye_T, params
        )

        aK_T = dAdV * S_a + A * (theta_model.dVrel(x, params) / params["V_0"]) * S_ad

        return aK_T

    @staticmethod
    def validate_parameters(params):
        # Check for all required keys
        expected_keys = [
            "debye_temperature_model",
            "anharmonic_prefactor_model",
            "anharmonic_thermal_model",
            "V_0",
            "T_0",
        ]
        for key in expected_keys:
            if key not in params:
                raise AttributeError(f"params dictionary must contain an '{key}' key")

        # Validate the three models:
        models = [
            params["debye_temperature_model"],
            params["anharmonic_prefactor_model"],
            params["anharmonic_thermal_model"],
        ]
        for model in models:
            model.validate_parameters(params)

        # Check that the required methods are present in the
        # debye_temperature_model
        expected_methods = ["value", "dVrel", "dVrel2"]
        for method in expected_methods:
            if not hasattr(params["debye_temperature_model"], method) or not callable(
                getattr(params["debye_temperature_model"], method)
            ):
                raise AttributeError(
                    f"params['debye_temperature_model'] must have a {method} method that takes arguments Vrel and params"
                )

        # Check that the required methods are present in the
        # anharmonic_prefactor_model
        expected_methods = ["value", "dVrel", "dVrel2"]
        for method in expected_methods:
            if not hasattr(
                params["anharmonic_prefactor_model"], method
            ) or not callable(getattr(params["anharmonic_prefactor_model"], method)):
                raise AttributeError(
                    f"params['anharmonic_prefactor_model'] must have a {method} method that takes arguments Vrel and params"
                )

        # Check that the required methods are present in the
        # anharmonic_thermal_model
        expected_methods = [
            "nondimensional_helmholtz_energy",
            "nondimensional_dhelmholtz_dTheta",
            "nondimensional_d2helmholtz_dTheta2",
            "nondimensional_entropy",
            "nondimensional_dentropy_dTheta",
            "nondimensional_heat_capacity",
        ]
        for method in expected_methods:
            if not hasattr(params["anharmonic_thermal_model"], method) or not callable(
                getattr(params["anharmonic_thermal_model"], method)
            ):
                raise AttributeError(
                    f"params['anharmonic_thermal_model'] must have a {method} method that takes arguments temperature, debye_T, and params"
                )
