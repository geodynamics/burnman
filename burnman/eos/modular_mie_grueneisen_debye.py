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


from . import debye
from . import equation_of_state as eos
from ..utils.math import bracket
from . import bukowinski_electronic as el
from .anharmonic_debye_pade import AnharmonicDebyePade as Anharmonic


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
            volume via ``__call__``, and its first and second derivatives with
            respect to relative volume via ``dVrel`` and ``d2dVrel2``.
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
        Debye_T = params["debye_temperature_model"](Vrel, params)
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
        if params["a_anh"] is not None:
            P += Anharmonic.pressure(temperature, volume, params)

        return P

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        KT_ref = params["reference_eos"].isothermal_bulk_modulus_reuss(
            0.0, params["T_0"], volume, params
        )
        Debye_T = params["debye_temperature_model"](volume / params["V_0"], params)
        V = volume
        d2ThetadV2 = (
            params["debye_temperature_model"].d2dVrel2(V / params["V_0"], params)
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
        if params["a_anh"] is not None:
            KT += Anharmonic.isothermal_bulk_modulus(temperature, volume, params)

        return KT

    def _molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        debye_T = params["debye_temperature_model"](volume / params["V_0"], params)
        C_v = debye.molar_heat_capacity_v(temperature, debye_T, params["n"])

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            bel_0 = params["bel_0"]
            gel = params["gel"]
            C_v += temperature * el.CVoverT(volume, params["V_0"], bel_0, gel)

        # If the material has an anharmonic component, add it
        if params["a_anh"] is not None:
            C_v += Anharmonic.heat_capacity_v(temperature, volume, params)

        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        debye_T = params["debye_temperature_model"](volume / params["V_0"], params)
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
        if params["a_anh"] is not None:
            aKT += Anharmonic.dSdV(temperature, volume, params)

        KT = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)

        return aKT / KT

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature
        of the mineral [J/K/mol]
        """
        Debye_T = params["debye_temperature_model"](volume / params["V_0"], params)
        S = debye.entropy(temperature, Debye_T, params["n"])

        # If the material is conductive, add the electronic contribution
        if params["bel_0"] is not None:
            S += el.entropy(
                temperature, volume, params["V_0"], params["bel_0"], params["gel"]
            )

        # If the material has an anharmonic component, add it
        if params["a_anh"] is not None:
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
        Debye_T = params["debye_temperature_model"](volume / params["V_0"], params)

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
        if params["a_anh"] is not None:
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

        if not hasattr(params["debye_temperature_model"], "__call__"):
            raise AttributeError(
                "params['debye_temperature_model'] must have a __call__ method to compute the Debye temperature as a function of Vrel and params."
            )
        if not hasattr(params["debye_temperature_model"], "dVrel"):
            raise AttributeError(
                "params['debye_temperature_model'] must have a dVrel method to compute the first derivative of the Debye temperature with respect to Vrel. Arguments should be Vrel and the params dictionary."
            )
        if not hasattr(params["debye_temperature_model"], "d2dVrel2"):
            raise AttributeError(
                "params['debye_temperature_model'] must have a d2dVrel2 method to compute the second derivative of the Debye temperature with respect to Vrel. Arguments should be Vrel and the params dictionary."
            )

        # Now check all the other required keys are in the dictionary
        expected_keys = ["molar_mass", "n", "T_0", "V_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Check if material is conductive
        if "bel_0" not in params:
            params["bel_0"] = None
            params["gel"] = None

        # Check if material has an anharmonic component
        if "a_anh" not in params:
            params["a_anh"] = None
            params["m_anh"] = None
