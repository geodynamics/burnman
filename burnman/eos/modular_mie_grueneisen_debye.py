# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:

    def jit(nopython=True):
        def decorator(fn):
            return fn

        return decorator


import numpy as np
from . import debye
from . import equation_of_state as eos
from ..utils.math import bracket
from ..utils.misc import copy_documentation
from . import bukowinski_electronic as el
from .anharmonic_debye import AnharmonicDebye as Anharmonic


class ModularMGD(eos.EquationOfState):
    """
    Base class for a modular Mie-Grueneisen-Debye equation of state.

    Parameters
    ----------
    params : dict
        Dictionary of parameters. Must contain all the parameters needed by
        the reference EoS and Debye temperature model, parameters for
        the electronic and anharmonic contributions (if necessary), plus:

        V_0 : float
            Reference volume [m^3].
        T_0 : float
            Reference temperature [K].
        n : int
            Number of formula units in the mineral.
        molar_mass : float
            Molar mass of the mineral [kg/mol].
        reference_eos : EquationOfState
            Instance of an equation of state that has methods for pressure,
            isothermal bulk modulus, and Gibbs energy.
        debye_temperature_model : callable
            An instance providing the Debye temperature as a function of relative
            volume via ``value``, and its first and second derivatives with
            respect to relative volume via ``dVrel`` and ``dVrel2``.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """

        def func(V):
            P = self.pressure(temperature, V, params)
            return P - pressure

        try:
            sol = bracket(func, params["V_0"], 1.0e-2 * params["V_0"])
        except:
            raise ValueError(
                "Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?"
            )
        volume = opt.brentq(func, sol[0], sol[1])
        return volume

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume
        [Pa]
        """
        Vrel = volume / params["V_0"]
        T_0 = params["T_0"]

        P_ref = params["reference_eos"].pressure(T_0, volume, params)
        Debye_T = params["debye_temperature_model"].value(Vrel, params)
        dThetadV = params["debye_temperature_model"].dVrel(Vrel, params) / params["V_0"]
        P_th = -dThetadV * (
            debye.dhelmholtz_dTheta(temperature, Debye_T, params["n"])
            - debye.dhelmholtz_dTheta(T_0, Debye_T, params["n"])
        )

        P = P_ref + P_th

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            bel_0 = params["bel_0"]
            gel = params["gel"]
            P += (
                0.5
                * gel
                * bel_0
                * pow(Vrel, gel)
                * (temperature * temperature - T_0 * T_0)
                / volume
            )

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            P += Anharmonic.pressure(temperature, volume, params)

        return P

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        KT_ref = params["reference_eos"].isothermal_bulk_modulus_reuss(
            0.0, params["T_0"], volume, params
        )
        Debye_T = params["debye_temperature_model"].value(
            volume / params["V_0"], params
        )
        V = volume
        d2ThetadV2 = (
            params["debye_temperature_model"].dVrel2(V / params["V_0"], params)
            / params["V_0"] ** 2
        )
        dThetadV = (
            params["debye_temperature_model"].dVrel(V / params["V_0"], params)
            / params["V_0"]
        )
        dFthdTheta = debye.dhelmholtz_dTheta(temperature, Debye_T, params["n"])
        dFthdTheta_T0 = debye.dhelmholtz_dTheta(params["T_0"], Debye_T, params["n"])
        d2FthdTheta2 = debye.d2helmholtz_dTheta2(temperature, Debye_T, params["n"])
        d2FthdTheta2_T0 = debye.d2helmholtz_dTheta2(params["T_0"], Debye_T, params["n"])
        d2FthdV2 = d2ThetadV2 * (dFthdTheta - dFthdTheta_T0) + dThetadV**2 * (
            d2FthdTheta2 - d2FthdTheta2_T0
        )

        KT = KT_ref + V * d2FthdV2

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            KT += volume * el.KToverV(
                temperature,
                volume,
                params["T_0"],
                params["V_0"],
                params["bel_0"],
                params["gel"],
            )

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            KT += Anharmonic.isothermal_bulk_modulus(temperature, volume, params)

        return KT

    def _molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        debye_T = params["debye_temperature_model"].value(
            volume / params["V_0"], params
        )
        C_v = debye.molar_heat_capacity_v(temperature, debye_T, params["n"])

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            bel_0 = params["bel_0"]
            gel = params["gel"]
            C_v += temperature * el.CVoverT(volume, params["V_0"], bel_0, gel)

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            C_v += Anharmonic.heat_capacity_v(temperature, volume, params)

        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        debye_T = params["debye_temperature_model"].value(
            volume / params["V_0"], params
        )
        dThetadV = (
            params["debye_temperature_model"].dVrel(volume / params["V_0"], params)
            / params["V_0"]
        )
        dSdTheta = debye.dentropy_dTheta(temperature, debye_T, params["n"])
        aKT = dSdTheta * dThetadV

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            bel_0 = params["bel_0"]
            gel = params["gel"]
            aKT += el.aKT(temperature, volume, params["V_0"], bel_0, gel)

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            aKT += Anharmonic.dSdV(temperature, volume, params)

        KT = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)

        return aKT / KT

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature
        of the mineral [J/K/mol]
        """
        Debye_T = params["debye_temperature_model"].value(
            volume / params["V_0"], params
        )
        S = debye.entropy(temperature, Debye_T, params["n"])

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            S += el.entropy(
                temperature, volume, params["V_0"], params["bel_0"], params["gel"]
            )

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            S += Anharmonic.entropy(temperature, volume, params)

        return S

    def _helmholtz_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature
        of the mineral [J/mol]
        """
        P_ref = params["reference_eos"].pressure(params["T_0"], volume, params)
        G_ref = params["reference_eos"].gibbs_energy(
            P_ref, params["T_0"], volume, params
        )
        F_ref = G_ref - volume * params["reference_eos"].pressure(
            params["T_0"], volume, params
        )
        Debye_T = params["debye_temperature_model"].value(
            volume / params["V_0"], params
        )

        F_th = debye.helmholtz_energy(
            temperature, Debye_T, params["n"]
        ) - debye.helmholtz_energy(params["T_0"], Debye_T, params["n"])

        F = F_ref + F_th

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            F += el.helmholtz(
                temperature,
                volume,
                params["T_0"],
                params["V_0"],
                params["bel_0"],
                params["gel"],
            )

        # If the material has an anharmonic component, add it
        if "anharmonic_thermal_model" in params:
            F += Anharmonic.helmholtz_energy(temperature, volume, params)

        return F

    # Derived properties from here
    # These functions can use pressure as an argument,
    # as pressure should have been determined self-consistently
    # by the point at which these functions are called.
    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter (-dlnT/dlnV at fixed entropy) :math:`[unitless]`
        The grueneisen parameter is a product of partial derivatives of the
        Helmholtz energy, but the product involves a division by the heat capacity
        (which goes to zero at low temperatures).
        Therefore we reset the temperature to a small value
        if it is too low, to avoid division by zero.
        """
        temperature = max(temperature, 1.0e-6)
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        C_v = self._molar_heat_capacity_v(pressure, temperature, volume, params)
        return alpha * K_T * volume / C_v

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        C_v = self._molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v + alpha * alpha * K_T * volume * temperature
        return C_p

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature
        of the mineral [J/mol]
        """
        G = (
            self._helmholtz_energy(pressure, temperature, volume, params)
            + pressure * volume
        )
        return G

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        # First, let's check the EoS parameters for the reference EoS
        params["reference_eos"].validate_parameters(params)

        # And check that the reference EoS has the required methods
        required_methods = ["pressure", "isothermal_bulk_modulus_reuss", "gibbs_energy"]
        for method in required_methods:
            if not hasattr(params["reference_eos"], method):
                raise AttributeError(
                    "Reference EoS does not have the required method: " + method
                )

        # Now check that the params dictionary contains a model for the Debye temperature.
        if "debye_temperature_model" not in params:
            raise KeyError(
                "params object missing 'debye_temperature_model' key, which must be a class instance."
            )

        # Check that the Debye temperature model has the required methods
        expected_methods = ["value", "dVrel", "dVrel2"]
        for method in expected_methods:
            if not hasattr(params["debye_temperature_model"], method):
                raise AttributeError(
                    f"params['debye_temperature_model'] must have a {method} method."
                )

        params["debye_temperature_model"].validate_parameters(params)

        # Make some checks if anharmonicity is included
        if "anharmonic_thermal_model" in params:
            Anharmonic.validate_parameters(params)

        # Now check all the other required keys are in the dictionary
        expected_keys = ["molar_mass", "n", "T_0", "V_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Check if material is conductive
        if "bel_0" not in params:
            params["bel_0"] = None
            params["gel"] = None


class ModularMGDWithAnharmonicity(ModularMGD):
    """
    This class extends the ModularMGD class to include anharmonicity effects
    according to a simplification of the model proposed by
    Wu and Wentzcovitch, 2009.

    The basis of the anharmonic model is the definition of a scaled volume, :math:`V'`:
    :math:`\\ln (V'/V) =  = - c \\ln (V/V_0)`.
    In this expression, :math:`V` is the target volume,
    :math:`V_0` is a first order approximation to the volume at the same pressure
    and reference temperature :math:`T_0`, and :math:`c` is an anharmonicity parameter
    provided in the params dictionary as `c_anh`.

    The anharmonic Helmholtz energy :math:`F` is related to the scaled volume by the equation:

    .. math::

        F(V,T) = F_h(V',T) + F_h(V,T_0) - F_h(V',T_0)

    where :math:`F_h` is the harmonic Helmholtz energy,
    potentially with electronic contributions.

    Note: This model is not the same as that published in
    Wu and Wentzcovitch (2009). The results are expected to be similar,
    but the :math:`c` parameter will in general need to be tweaked.
    This is because only a local approximation to the volume change
    between 0 K and the target temperature is used. This does not mean that
    the model is less able to capture the essential physics of the problem;
    indeed, the model of Wu and Wentzcovitch (2009) is only intended to be
    an effective ansatz.
    """

    def lnVoverV0_approx(self, temperature, volume, params):
        a_h = ModularMGD.thermal_expansivity(
            super(), np.nan, temperature, volume, params
        )
        KT_h = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, temperature, volume, params
        )
        KT0_h = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, params["T_0"], volume, params
        )
        P_th = a_h * KT_h * (temperature - params["T_0"])
        return P_th / KT0_h  # assume a constant bulk modulus

    @copy_documentation(ModularMGD.pressure)
    def pressure(self, temperature, volume, params):
        dV = 1e-5 * params["V_0"]

        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        dlnVoverV0_dV = (
            self.lnVoverV0_approx(temperature, volume + dV, params)
            - self.lnVoverV0_approx(temperature, volume - dV, params)
        ) / (2 * dV)

        lnVprimeoverV = -params["c_anh"] * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)

        dlnVprimeoverV_dV = -params["c_anh"] * dlnVoverV0_dV
        dVprime_dV = Vprime * (1.0 / volume + dlnVprimeoverV_dV)

        T0 = params["T_0"]
        P_V_T0 = ModularMGD.pressure(super(), T0, volume, params)
        P_Vp_T = ModularMGD.pressure(super(), temperature, Vprime, params)
        P_Vp_T0 = ModularMGD.pressure(super(), T0, Vprime, params)

        P = P_V_T0 + dVprime_dV * (P_Vp_T - P_Vp_T0)

        return P

    @copy_documentation(ModularMGD.isothermal_bulk_modulus_reuss)
    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        dV = 1e-5 * params["V_0"]

        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        lnVoverV0_plus_dV = self.lnVoverV0_approx(temperature, volume + dV, params)
        lnVoverV0_minus_dV = self.lnVoverV0_approx(temperature, volume - dV, params)

        dlnVoverV0_dV = (lnVoverV0_plus_dV - lnVoverV0_minus_dV) / (2 * dV)
        d2lnVoverV0_dV2 = (lnVoverV0_plus_dV - 2 * lnVoverV0 + lnVoverV0_minus_dV) / (
            dV**2
        )

        lnVprimeoverV = -params["c_anh"] * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)

        dlnVprimeoverV_dV = -params["c_anh"] * dlnVoverV0_dV
        dVprime_dV = Vprime * (1.0 / volume + dlnVprimeoverV_dV)
        d2Vprime_dV2 = Vprime * (
            (1.0 / volume + dlnVprimeoverV_dV) ** 2
            - 1.0 / volume**2
            - params["c_anh"] * d2lnVoverV0_dV2
        )

        T0 = params["T_0"]
        P_Vp_T = ModularMGD.pressure(super(), temperature, Vprime, params)
        P_Vp_T0 = ModularMGD.pressure(super(), T0, Vprime, params)

        KT_V_T0 = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, T0, volume, params
        )
        KT_Vp_T = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, temperature, Vprime, params
        )
        KT_Vp_T0 = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, T0, Vprime, params
        )

        KT = (
            KT_V_T0
            - volume * d2Vprime_dV2 * (P_Vp_T - P_Vp_T0)
            + volume / Vprime * (dVprime_dV) ** 2 * (KT_Vp_T - KT_Vp_T0)
        )
        return KT

    @copy_documentation(ModularMGD._molar_heat_capacity_v)
    def _molar_heat_capacity_v(self, pressure, temperature, volume, params):
        dT = 1e-2

        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        lnVoverV0_plus_dT = self.lnVoverV0_approx(temperature + dT, volume, params)
        lnVoverV0_minus_dT = self.lnVoverV0_approx(temperature - dT, volume, params)

        dlnVoverV0_dT = (lnVoverV0_plus_dT - lnVoverV0_minus_dT) / (2 * dT)
        d2lnVoverV0_dT2 = (lnVoverV0_plus_dT - 2 * lnVoverV0 + lnVoverV0_minus_dT) / (
            dT**2
        )

        c = params["c_anh"]
        lnVprimeoverV = -c * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)

        dlnVprimeoverV_dT = -c * dlnVoverV0_dT
        d2lnVprimeoverV_dT2 = -c * d2lnVoverV0_dT2

        dVprime_dT = Vprime * dlnVprimeoverV_dT
        d2Vprime_dT2 = Vprime * (dlnVprimeoverV_dT**2 + d2lnVprimeoverV_dT2)

        T0 = params["T_0"]

        Cvh_Vprime_T = ModularMGD._molar_heat_capacity_v(
            super(), np.nan, temperature, Vprime, params
        )
        P_Vprime_T = ModularMGD.pressure(super(), temperature, Vprime, params)
        P_Vprime_T0 = ModularMGD.pressure(super(), T0, Vprime, params)
        KT_Vprime_T = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, temperature, Vprime, params
        )
        KT_Vprime_T0 = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, T0, Vprime, params
        )
        dPdT_Vprime_T = (
            ModularMGD.thermal_expansivity(super(), np.nan, temperature, Vprime, params)
            * KT_Vprime_T
        )

        C_v = Cvh_Vprime_T + temperature * (
            2.0 * dPdT_Vprime_T * dVprime_dT
            - (P_Vprime_T0 - P_Vprime_T) * d2Vprime_dT2
            - (KT_Vprime_T - KT_Vprime_T0) / Vprime * dVprime_dT**2
        )
        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        dV = 1e-5 * params["V_0"]
        dT = 1e-2

        # lnVoverV0 derivatives
        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        lnVoverV0_plusV = self.lnVoverV0_approx(temperature, volume + dV, params)
        lnVoverV0_minusV = self.lnVoverV0_approx(temperature, volume - dV, params)
        lnVoverV0_plusT = self.lnVoverV0_approx(temperature + dT, volume, params)
        lnVoverV0_minusT = self.lnVoverV0_approx(temperature - dT, volume, params)
        lnVoverV0_plusV_plusT = self.lnVoverV0_approx(
            temperature + dT, volume + dV, params
        )
        lnVoverV0_plusV_minusT = self.lnVoverV0_approx(
            temperature - dT, volume + dV, params
        )
        lnVoverV0_minusV_plusT = self.lnVoverV0_approx(
            temperature + dT, volume - dV, params
        )
        lnVoverV0_minusV_minusT = self.lnVoverV0_approx(
            temperature - dT, volume - dV, params
        )

        dlnVoverV0_dV = (lnVoverV0_plusV - lnVoverV0_minusV) / (2 * dV)
        dlnVoverV0_dT = (lnVoverV0_plusT - lnVoverV0_minusT) / (2 * dT)
        d2lnVoverV0_dTdV = (
            lnVoverV0_plusV_plusT
            - lnVoverV0_plusV_minusT
            - lnVoverV0_minusV_plusT
            + lnVoverV0_minusV_minusT
        ) / (4 * dV * dT)

        c = params["c_anh"]
        lnVprimeoverV = -c * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)
        dlnVprime_dV = -c * dlnVoverV0_dV
        dlnVprime_dT = -c * dlnVoverV0_dT
        d2lnVprime_dTdV = -c * d2lnVoverV0_dTdV

        dVprime_dV = Vprime * (1.0 / volume + dlnVprime_dV)
        dVprime_dT = Vprime * dlnVprime_dT
        d2Vprime_dTdV = Vprime * (
            dlnVprime_dT * (1.0 / volume + dlnVprime_dV) + d2lnVprime_dTdV
        )

        T0 = params["T_0"]
        P_Vprime_T = ModularMGD.pressure(super(), temperature, Vprime, params)
        P_Vprime_T0 = ModularMGD.pressure(super(), T0, Vprime, params)
        KT_Vprime_T = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, temperature, Vprime, params
        )
        KT_Vprime_T0 = ModularMGD.isothermal_bulk_modulus_reuss(
            super(), np.nan, T0, Vprime, params
        )
        alphaKT_Vprime_T = (
            ModularMGD.thermal_expansivity(super(), np.nan, temperature, Vprime, params)
            * KT_Vprime_T
        )

        aKT = (
            alphaKT_Vprime_T * dVprime_dV
            - (KT_Vprime_T - KT_Vprime_T0) * dVprime_dV * dVprime_dT / Vprime
            + (P_Vprime_T - P_Vprime_T0) * d2Vprime_dTdV
        )
        KT = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)

        return aKT / KT

    def entropy(self, pressure, temperature, volume, params):
        dT = 1e-2

        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        dlnVoverV0_dT = (
            self.lnVoverV0_approx(temperature + dT, volume, params)
            - self.lnVoverV0_approx(temperature - dT, volume, params)
        ) / (2 * dT)

        c = params["c_anh"]
        lnVprimeoverV = -c * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)

        dlnVprimeoverV_dT = -c * dlnVoverV0_dT
        dVprime_dT = Vprime * dlnVprimeoverV_dT

        Sh_Vprime_T = ModularMGD.entropy(super(), np.nan, temperature, Vprime, params)
        P_Vprime_T = ModularMGD.pressure(super(), temperature, Vprime, params)
        P_Vprime_T0 = ModularMGD.pressure(super(), params["T_0"], Vprime, params)

        S = Sh_Vprime_T + (P_Vprime_T - P_Vprime_T0) * dVprime_dT
        return S

    def _helmholtz_energy(self, pressure, temperature, volume, params):
        lnVoverV0 = self.lnVoverV0_approx(temperature, volume, params)
        lnVprimeoverV = -params["c_anh"] * lnVoverV0
        Vprime = volume * np.exp(lnVprimeoverV)
        Fh_Vprime_T = ModularMGD._helmholtz_energy(
            super(), np.nan, temperature, Vprime, params
        )
        F_V_T0 = ModularMGD._helmholtz_energy(
            super(), np.nan, params["T_0"], volume, params
        )
        F_Vprime_T0 = ModularMGD._helmholtz_energy(
            super(), np.nan, params["T_0"], Vprime, params
        )

        F = Fh_Vprime_T + F_V_T0 - F_Vprime_T0
        return F

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        ModularMGD.validate_parameters(self, params)

        if "anharmonic_thermal_model" in params:
            raise ValueError(
                "The ModularMGDWithAnharmonicity class "
                "does not allow for a separate anharmonic "
                "contribution. Anharmonicity is incorporated "
                "through a single 'c_anh' parameter."
            )

        if "c_anh" not in params:
            raise ValueError(
                "The ModularMGDWithAnharmonicity class "
                "requires a 'c_anh' parameter to be set."
            )
