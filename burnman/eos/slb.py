# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
import scipy.optimize as opt
import warnings

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:

    def jit(fn):
        return fn


from . import birch_murnaghan as bm
from . import debye
from . import equation_of_state as eos
from . import bukowinski_electronic as el
from ..utils.math import bracket


@jit(nopython=True)
def _grueneisen_parameter_slb(V_0, volume, gruen_0, q_0):
    """global function with plain parameters so jit will work"""
    x = V_0 / volume
    f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
    a1_ii = 6.0 * gruen_0  # EQ 47
    a2_iikk = -12.0 * gruen_0 + 36.0 * gruen_0 * gruen_0 - 18.0 * q_0 * gruen_0  # EQ 47
    nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * f * f  # EQ 41
    return 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)


@jit(nopython=True)
def _delta_pressure(
    x,
    pressure,
    temperature,
    V_0,
    T_0,
    Debye_0,
    n,
    a1_ii,
    a2_iikk,
    b_iikk,
    b_iikkmm,
    bel_0,
    gel,
):
    f = 0.5 * (pow(V_0 / x, 2.0 / 3.0) - 1.0)
    nu_o_nu0_sq = 1.0 + a1_ii * f + 1.0 / 2.0 * a2_iikk * f * f
    debye_temperature = Debye_0 * np.sqrt(nu_o_nu0_sq)
    E_th = debye.thermal_energy(
        temperature, debye_temperature, n
    )  # thermal energy at temperature T
    E_th_ref = debye.thermal_energy(
        T_0, debye_temperature, n
    )  # thermal energy at reference temperature
    nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * f * f  # EQ 41
    gr = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)

    Pel = (
        0.5
        * gel
        * bel_0
        * np.power(x / V_0, gel)
        * (temperature * temperature - T_0 * T_0)
        / x
    )

    return (
        (1.0 / 3.0)
        * (pow(1.0 + 2.0 * f, 5.0 / 2.0))
        * ((b_iikk * f) + (0.5 * b_iikkmm * f * f))
        + gr * (E_th - E_th_ref) / x
        + Pel
        - pressure
    )  # EQ 21


class SLBBase(eos.EquationOfState):
    """
    Base class for the finite strain-Mie-Grueneiesen-Debye equation of state
    detailed in :cite:`Stixrude2005`.  For the most part the equations are
    all third order in strain, but see further the :class:`burnman.slb.SLB2`
    and :class:`burnman.slb.SLB3` classes.
    """

    def _debye_temperature(self, x, params):
        """
        Finite strain approximation for Debye Temperature [K]
        x = ref_vol/vol
        """
        f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
        a1_ii = 6.0 * params["grueneisen_0"]  # EQ 47
        a2_iikk = (
            -12.0 * params["grueneisen_0"]
            + 36.0 * pow(params["grueneisen_0"], 2.0)
            - 18.0 * params["q_0"] * params["grueneisen_0"]
        )  # EQ 47
        nu_o_nu0_sq = 1.0 + a1_ii * f + 1.0 / 2.0 * a2_iikk * f * f
        if nu_o_nu0_sq > 0.0:
            return params["Debye_0"] * np.sqrt(nu_o_nu0_sq)
        else:
            raise Exception(
                f"This volume (V = {1./x:.2f}*V_0) exceeds the "
                "valid range of the thermal "
                "part of the slb equation of state."
            )

    def volume_dependent_q(self, x, params):
        """
        Finite strain approximation for :math:`q`, the isotropic volume strain
        derivative of the grueneisen parameter.
        """
        f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
        a1_ii = 6.0 * params["grueneisen_0"]  # EQ 47
        a2_iikk = (
            -12.0 * params["grueneisen_0"]
            + 36.0 * pow(params["grueneisen_0"], 2.0)
            - 18.0 * params["q_0"] * params["grueneisen_0"]
        )  # EQ 47
        nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * f * f  # EQ 41
        gr = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)
        # avoids divide by zero if grueneisen_0 = 0.
        if np.abs(params["grueneisen_0"]) < 1.0e-10:
            q = 1.0 / 9.0 * (18.0 * gr - 6.0)
        else:
            q = (
                1.0
                / 9.0
                * (
                    18.0 * gr
                    - 6.0
                    - 1.0
                    / 2.0
                    / nu_o_nu0_sq
                    * (2.0 * f + 1.0)
                    * (2.0 * f + 1.0)
                    * a2_iikk
                    / gr
                )
            )
        return q

    def _isotropic_eta_s(self, x, params):
        """
        Finite strain approximation for :math:`eta_{s0}`, the isotropic shear
        strain derivative of the grueneisen parameter.
        """
        f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
        a2_s = -2.0 * params["grueneisen_0"] - 2.0 * params["eta_s_0"]  # EQ 47
        a1_ii = 6.0 * params["grueneisen_0"]  # EQ 47
        a2_iikk = (
            -12.0 * params["grueneisen_0"]
            + 36.0 * pow(params["grueneisen_0"], 2.0)
            - 18.0 * params["q_0"] * params["grueneisen_0"]
        )  # EQ 47
        nu_o_nu0_sq = 1.0 + a1_ii * f + (1.0 / 2.0) * a2_iikk * pow(f, 2.0)  # EQ 41
        gr = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f)
        # EQ 46 NOTE the typo from Stixrude 2005:
        eta_s = -gr - (
            1.0 / 2.0 * pow(nu_o_nu0_sq, -1.0) * pow((2.0 * f) + 1.0, 2.0) * a2_s
        )

        return eta_s

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """
        T_0 = params["T_0"]
        Debye_0 = params["Debye_0"]
        V_0 = params["V_0"]
        dV = 1.0e-2 * params["V_0"]
        n = params["n"]

        a1_ii = 6.0 * params["grueneisen_0"]  # EQ 47
        a2_iikk = (
            -12.0 * params["grueneisen_0"]
            + 36.0 * pow(params["grueneisen_0"], 2.0)
            - 18.0 * params["q_0"] * params["grueneisen_0"]
        )  # EQ 47

        b_iikk = 9.0 * params["K_0"]  # EQ 28
        b_iikkmm = 27.0 * params["K_0"] * (params["Kprime_0"] - 4.0)  # EQ 29z

        bel_0, gel = 0.0, 1.0
        if self.conductive:
            bel_0, gel = params["bel_0"], params["gel"]

        # Finding the volume at a given pressure requires a
        # root-finding scheme. Here we use brentq to find the root.

        # Root-finding using brentq requires bounds to be specified.
        # We do this using a bracketing function.
        args = (
            pressure,
            temperature,
            V_0,
            T_0,
            Debye_0,
            n,
            a1_ii,
            a2_iikk,
            b_iikk,
            b_iikkmm,
            bel_0,
            gel,
        )

        try:
            # The first attempt to find a bracket for
            # root finding uses V_0 as a starting point
            sol = bracket(_delta_pressure, V_0, dV, args)
        except Exception:
            # At high temperature, the naive bracketing above may
            # try a volume guess that exceeds the point at which the
            # bulk modulus goes negative at that temperature.
            # In this case, we try a more nuanced approach by
            # first finding the volume at which the bulk modulus goes
            # negative, and then either (a) raising an exception if the
            # desired pressure is less than the pressure at that volume,
            # or (b) using that pressure to create a better bracket for
            # brentq.
            def _K_T(V, T, params):
                return self.isothermal_bulk_modulus_reuss(0.0, T, V, params)

            sol_K_T = bracket(_K_T, V_0, dV, args=(temperature, params))
            V_crit = opt.brentq(
                _K_T, sol_K_T[0], sol_K_T[1], args=(temperature, params)
            )
            P_min = self.pressure(temperature, V_crit, params)
            if P_min > pressure:
                raise Exception(
                    "The desired pressure is not achievable "
                    "at this temperature. The minimum pressure "
                    f"achievable is {P_min:.2e} Pa."
                )
            else:
                try:
                    sol = bracket(_delta_pressure, V_crit - dV, dV, args)
                except Exception:
                    raise Exception(
                        "Cannot find a volume, perhaps you are "
                        "outside of the range of validity for "
                        "the equation of state?"
                    )

        return opt.brentq(_delta_pressure, sol[0], sol[1], args=args)

    def pressure(self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume
        [Pa]
        """
        debye_T = self._debye_temperature(params["V_0"] / volume, params)
        gr = _grueneisen_parameter_slb(
            params["V_0"], volume, params["grueneisen_0"], params["q_0"]
        )
        # does not depend on pressure
        # thermal energy at temperature T
        E_th = debye.thermal_energy(temperature, debye_T, params["n"])
        # thermal energy at reference temperature
        E_th_ref = debye.thermal_energy(params["T_0"], debye_T, params["n"])

        b_iikk = 9.0 * params["K_0"]  # EQ 28
        b_iikkmm = 27.0 * params["K_0"] * (params["Kprime_0"] - 4.0)  # EQ 29
        f = 0.5 * (pow(params["V_0"] / volume, 2.0 / 3.0) - 1.0)  # EQ 24
        P = (1.0 / 3.0) * (pow(1.0 + 2.0 * f, 5.0 / 2.0)) * (
            (b_iikk * f) + (0.5 * b_iikkmm * pow(f, 2.0))
        ) + gr * (
            E_th - E_th_ref
        ) / volume  # EQ 21
        if self.conductive:
            P += el.pressure(
                temperature,
                volume,
                params["T_0"],
                params["V_0"],
                params["bel_0"],
                params["gel"],
            )
        return P

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        T_0 = params["T_0"]
        debye_T = self._debye_temperature(params["V_0"] / volume, params)
        gr = _grueneisen_parameter_slb(
            params["V_0"], volume, params["grueneisen_0"], params["q_0"]
        )

        # thermal energy at temperature T
        E_th = debye.thermal_energy(temperature, debye_T, params["n"])
        # thermal energy at reference temperature
        E_th_ref = debye.thermal_energy(T_0, debye_T, params["n"])

        # heat capacity at temperature T
        C_v = debye.molar_heat_capacity_v(temperature, debye_T, params["n"])
        # heat capacity at reference temperature
        C_v_ref = debye.molar_heat_capacity_v(T_0, debye_T, params["n"])

        q = self.volume_dependent_q(params["V_0"] / volume, params)

        K = (
            bm.bulk_modulus(volume, params)
            + (gr + 1.0 - q) * (gr / volume) * (E_th - E_th_ref)
            - (pow(gr, 2.0) / volume) * (C_v * temperature - C_v_ref * T_0)
        )

        if self.conductive:
            K += volume * el.KToverV(
                temperature,
                volume,
                params["T_0"],
                params["V_0"],
                params["bel_0"],
                params["gel"],
            )
        return K

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus. :math:`[Pa]`
        """
        T_0 = params["T_0"]
        debye_T = self._debye_temperature(params["V_0"] / volume, params)
        eta_s = self._isotropic_eta_s(params["V_0"] / volume, params)

        E_th = debye.thermal_energy(temperature, debye_T, params["n"])
        E_th_ref = debye.thermal_energy(T_0, debye_T, params["n"])

        if self.order == 2:
            return (
                bm.shear_modulus_second_order(volume, params)
                - eta_s * (E_th - E_th_ref) / volume
            )
        elif self.order == 3:
            return (
                bm.shear_modulus_third_order(volume, params)
                - eta_s * (E_th - E_th_ref) / volume
            )
        else:
            raise NotImplementedError("")

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        debye_T = self._debye_temperature(params["V_0"] / volume, params)
        C_v = debye.molar_heat_capacity_v(temperature, debye_T, params["n"])

        if self.conductive:
            C_v += temperature * el.CVoverT(
                volume, params["V_0"], params["bel_0"], params["gel"]
            )
        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        debye_T = self._debye_temperature(params["V_0"] / volume, params)
        C_v = debye.molar_heat_capacity_v(temperature, debye_T, params["n"])
        gr_slb = _grueneisen_parameter_slb(
            params["V_0"], volume, params["grueneisen_0"], params["q_0"]
        )
        K = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        alpha = gr_slb * C_v / volume / K

        if self.conductive:
            aKTel = el.aKT(
                temperature, volume, params["V_0"], params["bel_0"], params["gel"]
            )
            alpha += aKTel / K
        return alpha

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature
        of the mineral [J/K/mol]
        """
        Debye_T = self._debye_temperature(params["V_0"] / volume, params)
        S = debye.entropy(temperature, Debye_T, params["n"])

        if self.conductive:
            S += el.entropy(
                temperature, volume, params["V_0"], params["bel_0"], params["gel"]
            )
        return S

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature
        of the mineral [J/mol]
        """
        x = params["V_0"] / volume
        f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
        Debye_T = self._debye_temperature(params["V_0"] / volume, params)

        F_quasiharmonic = debye.helmholtz_free_energy(
            temperature, Debye_T, params["n"]
        ) - debye.helmholtz_free_energy(params["T_0"], Debye_T, params["n"])

        b_iikk = 9.0 * params["K_0"]  # EQ 28
        b_iikkmm = 27.0 * params["K_0"] * (params["Kprime_0"] - 4.0)  # EQ 29

        F = (
            params["F_0"]
            + 0.5 * b_iikk * f * f * params["V_0"]
            + (1.0 / 6.0) * params["V_0"] * b_iikkmm * f * f * f
            + F_quasiharmonic
        )

        if self.conductive:
            F += el.helmholtz(
                temperature,
                volume,
                params["T_0"],
                params["V_0"],
                params["bel_0"],
                params["gel"],
            )
        return F

    # Derived properties from here
    # These functions can use pressure as an argument,
    # as pressure should have been determined self-consistently
    # by the point at which these functions are called.
    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]`
        """
        if self.conductive:
            temperature = max(temperature, 1.0e-6)
            K_T = self.isothermal_bulk_modulus_reuss(
                pressure, temperature, volume, params
            )
            alpha = self.thermal_expansivity(pressure, temperature, volume, params)
            C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
            return alpha * K_T * volume / C_v
        else:
            return _grueneisen_parameter_slb(
                params["V_0"], volume, params["grueneisen_0"], params["q_0"]
            )

    def isentropic_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus. :math:`[Pa]`
        """
        K_T = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        return K_T * (1.0 + gr * alpha * temperature)

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v * (1.0 + gr * alpha * temperature)
        return C_p

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy at the pressure and temperature
        of the mineral [J/mol]
        """
        return self.helmholtz_free_energy(
            pressure, temperature, volume, params
        ) + temperature * self.entropy(pressure, temperature, volume, params)

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature
        of the mineral [J/mol]
        """
        G = (
            self.helmholtz_free_energy(pressure, temperature, volume, params)
            + pressure * volume
        )
        return G

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature
        of the mineral [J/mol]
        """

        return (
            self.helmholtz_free_energy(pressure, temperature, volume, params)
            + temperature * self.entropy(pressure, temperature, volume, params)
            + pressure * volume
        )

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if "T_0" not in params:
            params["T_0"] = 300.0

        # If eta_s_0 is not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without it,
        # so just add it to the dictionary as nan
        # The same goes for the standard state Helmholtz free energy
        if "eta_s_0" not in params:
            params["eta_s_0"] = float("nan")
        if "F_0" not in params:
            params["F_0"] = float("nan")

        # First, let's check the EoS parameters for Tref
        bm.BirchMurnaghanBase.validate_parameters(bm.BirchMurnaghanBase(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["molar_mass", "n", "Debye_0", "grueneisen_0", "q_0", "eta_s_0"]
        if self.conductive:
            expected_keys.extend(["bel_0", "gel"])

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
        if params["Debye_0"] < 1.0 or params["Debye_0"] > 10000.0:
            warnings.warn("Unusual value for Debye_0", stacklevel=2)
        if params["grueneisen_0"] < -1.0 or params["grueneisen_0"] > 10.0:
            warnings.warn("Unusual value for grueneisen_0", stacklevel=2)
        if params["q_0"] < -20.0 or params["q_0"] > 20.0:
            warnings.warn("Unusual value for q_0", stacklevel=2)
        if params["eta_s_0"] < -10.0 or params["eta_s_0"] > 10.0:
            warnings.warn("Unusual value for eta_s_0", stacklevel=2)


class SLB3(SLBBase):
    """
    SLB equation of state with third order finite strain expansion for the
    shear modulus (this should be preferred, as it is more thermodynamically
    consistent).
    """

    def __init__(self):
        self.order = 3
        self.conductive = False


class SLB2(SLBBase):
    """
    SLB equation of state with second order finite strain expansion for the
    shear modulus.  In general, this should not be used, but sometimes
    shear modulus data is fit to a second order equation of state.  In that
    case, you should use this.  The moral is, be careful!
    """

    def __init__(self):
        self.order = 2
        self.conductive = False


class SLB3Conductive(SLBBase):
    """
    SLB equation of state with third order finite strain expansion for the
    shear modulus and a contribution to the Helmholtz free energy
    that arises from the thermal excitation of electrons (Bukowinski, 1977).
    """

    def __init__(self):
        self.order = 3
        self.conductive = True
