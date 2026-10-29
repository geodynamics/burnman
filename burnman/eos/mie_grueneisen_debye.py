# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import scipy.optimize as opt
import warnings

from . import equation_of_state as eos
from . import birch_murnaghan as bm
from . import debye
from .. import constants
from ..utils.math import bracket


class MGDBase(eos.EquationOfState):
    """
    Base class for a generic finite-strain Mie-Grueneisen-Debye
    equation of state.  References for this can be found in many
    places, such as :cite:`Shim2002` and :cite:`Jackson1996`.
    Here we mostly follow the appendices of :cite:`Matas2007`.
    Of particular note is the thermal correction to the shear modulus, which
    was developed by :cite:`Hama1998`.
    """

    def _grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume (EQ B6)
        """
        return self._grueneisen_parameter(params["V_0"] / volume, params)

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ B7
        """
        T_0 = params["T_0"]
        func = (
            lambda x: bm.pressure_third_order(params["V_0"] / x, params)
            + self._thermal_pressure(temperature, x, params)
            - self._thermal_pressure(T_0, x, params)
            - pressure
        )
        try:
            sol = bracket(func, params["V_0"], 1.0e-2 * params["V_0"])
        except:
            raise ValueError(
                "Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?"
            )
        return opt.brentq(func, sol[0], sol[1])

    def isothermal_bulk_modulus_reuss(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B8
        """
        T_0 = params["T_0"]
        K_T = (
            bm.bulk_modulus_third_order(volume, params)
            + self._thermal_bulk_modulus(temperature, volume, params)
            - self._thermal_bulk_modulus(T_0, volume, params)
        )  # EQB13
        return K_T

    # calculate the mgd shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B11
        """
        T_0 = params["T_0"]
        if self.order == 2:
            return (
                bm.shear_modulus_second_order(volume, params)
                + self._thermal_shear_modulus(temperature, volume, params)
                - self._thermal_shear_modulus(T_0, volume, params)
            )  # EQ B11
        elif self.order == 3:
            return (
                bm.shear_modulus_third_order(volume, params)
                + self._thermal_shear_modulus(temperature, volume, params)
                - self._thermal_shear_modulus(T_0, volume, params)
            )  # EQ B11
        else:
            raise NotImplementedError("")

    # heat capacity at constant volume
    def _molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature, and volume [J/K/mol]
        """
        Debye_T = self._debye_temperature(params["V_0"] / volume, params)
        C_v = debye.molar_heat_capacity_v(temperature, Debye_T, params["n"])
        return C_v

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        """
        C_v = self._molar_heat_capacity_v(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(params["V_0"] / volume, params)
        K = self.isothermal_bulk_modulus_reuss(pressure, temperature, volume, params)
        alpha = gr * C_v / K / volume
        return alpha

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self._grueneisen_parameter(params["V_0"] / volume, params)
        C_v = self._molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v * (1.0 + gr * alpha * temperature)
        return C_p

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        T_0 = params["T_0"]
        return (
            bm.pressure_third_order(params["V_0"] / volume, params)
            + self._thermal_pressure(temperature, volume, params)
            - self._thermal_pressure(T_0, volume, params)
        )

    def gibbs_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        """
        G = (
            self._helmholtz_energy(pressure, temperature, volume, params)
            + pressure * volume
        )
        return G

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """
        Debye_T = self._debye_temperature(params["V_0"] / volume, params)
        S = debye.entropy(temperature, Debye_T, params["n"])
        return S

    def _helmholtz_energy(self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        """
        x = params["V_0"] / volume
        f = 1.0 / 2.0 * (pow(x, 2.0 / 3.0) - 1.0)
        b_iikk = 9.0 * params["K_0"]  # EQ 28, SLB2005
        b_iikkmm = 27.0 * params["K_0"] * (params["Kprime_0"] - 4.0)  # EQ 29, SLB2005

        F_pressure = (
            0.5 * b_iikk * f * f * params["V_0"]
            + (1.0 / 6.0) * params["V_0"] * b_iikkmm * f * f * f
        )

        Debye_T = self._debye_temperature(params["V_0"] / volume, params)
        F_thermal = debye.helmholtz_energy(
            temperature, Debye_T, params["n"]
        ) - debye.helmholtz_energy(params["T_0"], Debye_T, params["n"])

        return params["F_0"] + F_pressure + F_thermal

    # calculate the thermal correction to the shear modulus as a function of
    # V, T
    def _thermal_shear_modulus(self, T, V, params):
        if T > 1.0e-10:
            gr = self._grueneisen_parameter(params["V_0"] / V, params)
            Debye_T = self._debye_temperature(params["V_0"] / V, params)
            G_th = (
                3.0
                / 5.0
                * (
                    self._thermal_bulk_modulus(T, V, params)
                    - 6
                    * constants.gas_constant
                    * T
                    * params["n"]
                    / V
                    * gr
                    * debye.debye_fn(Debye_T / T)
                )
            )  # EQ B10
            return G_th
        else:
            return 0.0

    # compute the Debye temperature in K.  Takes the
    # parameter x, which is V_0/V (molar volumes).
    # Depends on the reference grueneisen parameter,
    # the reference Debye temperature, and the factor
    # q_0, see Matas eq B6
    def _debye_temperature(self, x, params):
        return params["Debye_0"] * np.exp(
            (params["grueneisen_0"] - self._grueneisen_parameter(x, params))
            / params["q_0"]
        )

    # compute the grueneisen parameter with depth, according
    # to q_0.  Takes x=V_0/V. See Matas eq B6
    def _grueneisen_parameter(self, x, params):
        return params["grueneisen_0"] * pow(1.0 / x, params["q_0"])

    # calculate isotropic thermal pressure, see
    # Matas et. al. (2007) eq B4
    def _thermal_pressure(self, T, V, params):
        Debye_T = self._debye_temperature(params["V_0"] / V, params)
        gr = self._grueneisen_parameter(params["V_0"] / V, params)
        P_th = gr * debye.thermal_energy(T, Debye_T, params["n"]) / V
        return P_th

    # calculate the thermal correction for the mgd
    # bulk modulus (see matas et al, 2007)
    def _thermal_bulk_modulus(self, T, V, params):
        if T > 1.0e-10:
            gr = self._grueneisen_parameter(params["V_0"] / V, params)
            Debye_T = self._debye_temperature(params["V_0"] / V, params)
            K_th = (
                3.0
                * params["n"]
                * constants.gas_constant
                * T
                / V
                * gr
                * (
                    (1.0 - params["q_0"] - 3.0 * gr) * debye.debye_fn(Debye_T / T)
                    + 3.0 * gr * (Debye_T / T) / (np.exp(Debye_T / T) - 1.0)
                )
            )  # EQ B5
            return K_th
        else:
            return 0.0

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if "T_0" not in params:
            params["T_0"] = 300.0
        if "F_0" not in params:
            params["F_0"] = 0.0

        # First, let's check the EoS parameters for Tref
        bm.BirchMurnaghanBase.validate_parameters(bm.BirchMurnaghanBase(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ["molar_mass", "n", "Debye_0", "grueneisen_0", "q_0"]
        for k in expected_keys:
            if k not in params:
                raise KeyError("params object missing parameter : " + k)

        # Finally, check that the values are reasonable.
        if params["T_0"] < 0.0:
            warnings.warn("Unusual value for T_0", stacklevel=2)
        if params["molar_mass"] < 0.001 or params["molar_mass"] > 1.0:
            warnings.warn("Unusual value for molar_mass", stacklevel=2)
        if params["n"] < 1.0 or params["n"] > 1000.0:
            warnings.warn("Unusual value for n", stacklevel=2)
        if params["Debye_0"] < 1.0 or params["Debye_0"] > 10000.0:
            warnings.warn("Unusual value for Debye_0", stacklevel=2)
        if params["grueneisen_0"] < 0.0 or params["grueneisen_0"] > 10.0:
            warnings.warn("Unusual value for grueneisen_0", stacklevel=2)
        if params["q_0"] < -10.0 or params["q_0"] > 10.0:
            warnings.warn("Unusual value for q_0", stacklevel=2)


class MGD3(MGDBase):
    """
    Mie-Grueneisen-Debye equation of state with third order finite strain
    expansion for the shear modulus.

    The Debye model for the Helmholtz free energy can be written as follows :cite:`Matas2007`

    .. math::
        \\mathcal{F} &= \\frac{9nRT}{V}\\frac{1}{x^3} \\int_{0}^{x} \\xi^2 \\ln (1-e^{-\\xi}) d\\xi, \\\\
        x &= \\theta / T, \\\\
        \\theta &= \\theta_0 \\exp \\left( \\frac{\\gamma_0-\\gamma}{q_0} \\right), \\\\
        \\gamma &= \\gamma_0 \\left( \\frac{V}{V_0} \\right)^{q_0}

    where :math:`\\theta` is the Debye temperature and :math:`\\gamma` is the Grüneisen parameter.

    Using thermodynamic relations we can derive equations for the thermal pressure and bulk modulus

    .. math::
        P_{th}(V,T) &= - \\frac{\\partial \\mathcal{F(V, T)}}{\\partial V}, \\\\
        &= \\frac{3 n \\gamma R T}{V} D(x), \\\\
        K_{th}(V,T) &= -V \\frac{\\partial P(V, T)}{\\partial V}, \\\\
        &= \\frac{3 n \\gamma R T}{V} \\gamma \\left[ (1-q_0 - 3 \\gamma) D(x) + 3\\gamma \\frac{x}{e^x - 1} \\right], \\\\
        D(x) &= \\frac{3}{x^3} \\int_{0}^{x} \\frac{\\xi^3}{e^{\\xi} - 1} d\\xi

    The thermal shear correction used in BurnMan was developed by :cite:`HS1998`

    .. math::
        G_{th}(V,T) = \\frac{3}{5} \\left[ K_{th} (V, T) - 2\\frac{3nRT}{V}\\gamma D(x) \\right]


    The total pressure, bulk and shear moduli can be calculated from the following sums

    .. math::
        P(V, T) &= P_{\\textrm{ref}}(V, T_0) + P_{th}(V, T) - P_{th}(V, T_0), \\\\
        K(V, T) &= K_{\\textrm{ref}}(V, T_0) + K_{th}(V, T) - K_{th}(V, T_0), \\\\
        G(V, T) &= G_{\\textrm{ref}}(V, T_0) + G_{th}(V, T) - G_{th}(V, T_0)

    This equation of state is substantially the same as that in :cite:`Stixrude2005`.
    The primary differences are in the thermal correction to the shear modulus and
    in the volume dependences of the Debye temperature and the Gruneisen parameter.

    .. list-table::
        :widths: 25 75 20
        :header-rows: 1

        * - Parameter
          - Description
          - Units
        * - ``F_0``
          - Reference Helmholtz free energy.
          - :math:`\\text{J/mol}`
        * - ``V_0``
          - Reference volume.
          - :math:`\\text{m}^3`
        * - ``K_0``
          - Reference bulk modulus.
          - :math:`\\text{Pa}`
        * - ``Kprime_0``
          - Pressure derivative of the bulk modulus at reference pressure.
          - Dimensionless
        * - ``molar_mass``
          - Molar mass of the mineral.
          - :math:`\\text{kg/mol}`
        * - ``n``
          - Number of atoms per formula unit.
          - Dimensionless
        * - ``Debye_0``
          - Debye temperature at reference state.
          - :math:`\\text{K}`
        * - ``grueneisen_0``
          - Grüneisen parameter at reference state.
          - Dimensionless
        * - ``q_0``
          - Volume dependence of the Grüneisen parameter.
          - Dimensionless
    """

    def __init__(self):
        self.order = 3


class MGD2(MGDBase):
    """
    MGD equation of state with second order finite strain expansion for the
    shear modulus.  In general, this should not be used, but sometimes
    shear modulus data is fit to a second order equation of state.  In that
    case, you should use this.  The moral is, be careful!
    """

    def __init__(self):
        self.order = 2
