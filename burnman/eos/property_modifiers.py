# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import scipy.optimize as opt
from ..constants import gas_constant
from . import debye, einstein
from collections import OrderedDict

"""
Functions for modifying the thermodynamic properties of minerals
via excess Gibbs energies that are added to the base energy calculated
via the equation of state. These excess Gibbs energies are a function of
pressure and temperature.
"""


def landau_excesses(pressure, temperature, params):
    """
    Applies a tricritical Landau correction to the properties
    of an endmember which undergoes a displacive phase transition. These
    transitions are not associated with an activation energy, and therefore
    they occur rapidly compared with seismic wave propagation.

    This correction follows :cite:`Putnis1992`, and is done relative to
    the completely *ordered* state (at 0 K).
    It therefore differs in implementation from both :cite:`Stixrude2011` and
    :cite:`HP2011`, who compute properties relative to
    the completely disordered state and standard states respectively.
    The current implementation is preferred, as the excess
    entropy (and heat capacity) terms are equal to zero at 0 K.

    .. math::
        Tc = Tc_0 + \\frac{V_D P}{S_D}

    If the temperature is above the critical temperature,
    Q (the order parameter) is equal to zero, and the Gibbs free energy is simply that of the disordered phase:

    .. math::
        \\mathcal{G}_{\\textrm{dis}} = -S_D \\left( \\left( T - Tc \\right) + \\frac{Tc_0}{3} \\right), \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P}_{\\textrm{dis}} = V_D, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T}_{\\textrm{dis}} = -S_D

    If temperature is below the critical temperature, Q is between 0 and 1.
    The Gibbs energy can be described thus:

    .. math::
        Q^2 = \\sqrt{\\left( 1 - \\frac{T}{Tc} \\right)}, \\\\
        \\mathcal{G} = S_D \\left((T - Tc) Q^2 + \\frac{Tc_0 Q^6}{3} \\right) + \\mathcal{G}_{\\textrm{dis}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = - V_D Q^2 \\left(1 + \\frac{T}{2 Tc} \\left(1. - \\frac{Tc_0}{Tc} \\right) \\right) + \\frac{\\partial \\mathcal{G}}{\\partial P}_{\\textrm{dis}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = S_D Q^2 \\left(\\frac{3}{2} - \\frac{Tc_0}{2 Tc} \\right) + \\frac{\\partial \\mathcal{G}}{\\partial T}_{\\textrm{dis}}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = V_D^2 \\frac{T}{S_D Tc^2 Q^2} \\left( \\frac{T}{4 Tc} \\left(1. + \\frac{Tc_0}{Tc} \\right) + Q^4 \\left(1. - \\frac{Tc_0}{Tc} \\right) - 1 \\right), \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = - \\frac{S_D}{Tc Q^2} \\left(\\frac{3}{4} - \\frac{Tc_0}{4 Tc} \\right), \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P \\partial T} = \\frac{V_D}{2 Tc Q^2} \\left(1 + \\left(\\frac{T}{2 Tc} - Q^4 \\right) \\left(1 - \\frac{Tc_0}{Tc} \\right) \\right)

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``Tc_0``
          - Critical temperature at reference pressure (K)
        * - ``S_D``
          - Entropy parameter (J/mol/K)
        * - ``V_D``
          - Volume parameter (m³/mol)
    """

    Tc = params["Tc_0"] + params["V_D"] * pressure / params["S_D"]

    G_disordered = -params["S_D"] * ((temperature - Tc) + params["Tc_0"] / 3.0)
    dGdT_disordered = -params["S_D"]
    dGdP_disordered = params["V_D"]

    if temperature < Tc:
        # Wolfram input to check partial differentials
        # x = T, y = P, a = S, c = Tc0, d = V
        # D[D[a ((x - c - d*y/a)*(1 - x/(c + d*y/a))^0.5 + c/3*(1 - x/(c +
        # d*y/a))^1.5), x], x]
        Q2 = np.sqrt(1.0 - temperature / Tc)
        G = (
            params["S_D"]
            * ((temperature - Tc) * Q2 + params["Tc_0"] * Q2 * Q2 * Q2 / 3.0)
            + G_disordered
        )
        dGdP = (
            -params["V_D"]
            * Q2
            * (1.0 + 0.5 * temperature / Tc * (1.0 - params["Tc_0"] / Tc))
            + dGdP_disordered
        )
        dGdT = params["S_D"] * Q2 * (1.5 - 0.5 * params["Tc_0"] / Tc) + dGdT_disordered
        d2GdP2 = (
            params["V_D"]
            * params["V_D"]
            * temperature
            / (params["S_D"] * Tc * Tc * Q2)
            * (
                temperature * (1.0 + params["Tc_0"] / Tc) / (4.0 * Tc)
                + Q2 * Q2 * (1.0 - params["Tc_0"] / Tc)
                - 1.0
            )
        )
        d2GdT2 = -params["S_D"] / (Tc * Q2) * (0.75 - 0.25 * params["Tc_0"] / Tc)
        d2GdPdT = (
            params["V_D"]
            / (2.0 * Tc * Q2)
            * (1.0 + (temperature / (2.0 * Tc) - Q2 * Q2) * (1.0 - params["Tc_0"] / Tc))
        )

    else:
        Q2 = 0.0
        G = G_disordered
        dGdT = dGdT_disordered
        dGdP = dGdP_disordered
        d2GdT2 = 0.0
        d2GdP2 = 0.0
        d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, {"Q": np.sqrt(Q2)})


def landau_slb_2022_excesses(pressure, temperature, params):
    """
    Applies a tricritical Landau correction to the properties
    of an endmember which undergoes a displacive phase transition.
    This correction follows :cite:`Stixrude2022`,
    and is done relative to the state with order parameter :math:`Q=1`.

    The order parameter of this formulation can exceed one,
    at odds with :cite:`Putnis1992`, but in better agreement with
    atomic intuition :cite:`Stixrude2022`.
    Nevertheless, this implementation is still not perfect,
    as the excess entropy (and heat capacity) terms are not equal
    to zero at 0 K. :math:`Q` is limited to values less than or equal to 2
    to avoid unrealistic stabilisation at ultrahigh pressure.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``Tc_0``
          - Critical temperature at reference pressure (K)
        * - ``S_D``
          - Entropy parameter (J/mol/K)
        * - ``V_D``
          - Volume parameter (m³/mol)

    """

    Tc = params["Tc_0"] + params["V_D"] * pressure / params["S_D"]

    G_disordered = -params["S_D"] * ((temperature - Tc) + params["Tc_0"] / 3.0)
    dGdT_disordered = -params["S_D"]
    dGdP_disordered = params["V_D"]

    if temperature < Tc:
        Q2 = np.sqrt((Tc - temperature) / params["Tc_0"])

        if Q2 < 4.0:
            # Wolfram input to check partial differentials
            # x = T, y = P, a = S, c = Tc0, d = V
            # D[D[a ((x - c - d*y/a)*((c + d*y/a - x)/c)^0.5
            # + c/3*((c + d*y/a - x)/c)^1.5), x], x]
            # where Q2 = ((c + d*y/a - x)/c)^0.5
            G = (
                params["S_D"]
                * ((temperature - Tc) * Q2 + params["Tc_0"] * Q2 * Q2 * Q2 / 3.0)
                + G_disordered
            )
            dGdP = -params["V_D"] * Q2 + dGdP_disordered
            dGdT = params["S_D"] * Q2 + dGdT_disordered
            d2GdP2 = (
                -0.5
                * params["V_D"]
                * params["V_D"]
                / (params["S_D"] * params["Tc_0"] * Q2)
            )
            d2GdT2 = -0.5 * params["S_D"] / (Tc * Q2)
            d2GdPdT = 0.5 * params["V_D"] / (Tc * Q2)
        else:
            # Wolfram input to check partial differentials
            # x = T, y = P, a = S, c = Tc0, d = V
            # D[D[a ((x - c - d*y/a)*4 + c/3*64), x], x]
            G = (
                params["S_D"] * ((temperature - Tc) * 4.0 + params["Tc_0"] * 64.0 / 3.0)
                + G_disordered
            )
            dGdP = -params["V_D"] * 4.0 + dGdP_disordered
            dGdT = params["S_D"] * 4.0 + dGdT_disordered
            d2GdP2 = 0.0
            d2GdT2 = 0.0
            d2GdPdT = 0.0

    else:
        Q2 = 0.0
        G = G_disordered
        dGdT = dGdT_disordered
        dGdP = dGdP_disordered
        d2GdT2 = 0.0
        d2GdP2 = 0.0
        d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, {"Q": np.sqrt(Q2)})


def landau_hp_excesses(pressure, temperature, params):
    """
    Applies a tricritical Landau correction similar to that described
    above. However, this implementation follows :cite:`HP2011`,
    who compute properties relative to the standard state.

    Includes the correction published within landaunote.pdf
    (Holland, pers. comm), which 'corrects' the terms involving
    the critical temperature Tc / Tc*.

    Note that this implementation allows the order parameter Q to be
    greater than one.

    .. math::
        Tc = Tc0 + \\frac{V_D P}{S_D}

    If the temperature is above the critical temperature,
    Q (the order parameter) is equal to zero. Otherwise

    .. math::
        Q^2 = \\sqrt{\\left( \\frac{Tc - T}{Tc0} \\right)}

    .. math::
        \\mathcal{G} = Tc_0 S_D \\left( Q_0^2 - \\frac{Q_0 ^ 6}{3} \\right) \\\\
            - S_D \\left( Tc Q^2 - Tc_0 \\frac{Q ^ 6}{3} \\right) \\\\
            - T S_D \\left( Q_0^2 - Q^2 \\right) + P V_D Q_0^2, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = -V_D \\left( Q^2 - Q_0^2
        \\right), \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = S_D \\left( Q^2 - Q_0^2
        \\right), \\\\

    The second derivatives of the Gibbs free energy are only non-zero if
    the order parameter exceeds zero. Then

    .. math::
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2}  = -\\frac{V_D^2}{2
        S_D Tc_0 Q^2}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2}  =
        -\\frac{S_D}{2 Tc_0 Q^2}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P \\partial T} =
        \\frac{V_D}{2 Tc_0 Q^2}

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``P_0``
          - Reference pressure (Pa)
        * - ``T_0``
          - Reference temperature (K)
        * - ``Tc_0``
          - Critical temperature at reference pressure (K)
        * - ``S_D``
          - Entropy parameter (J/mol/K)
        * - ``V_D``
          - Volume parameter (m³/mol)
    """

    P = pressure
    T = temperature

    if params["T_0"] < params["Tc_0"]:
        Q_0 = np.power((params["Tc_0"] - params["T_0"]) / params["Tc_0"], 0.25)
    else:
        Q_0 = 0.0

    Tc = params["Tc_0"] + params["V_D"] * (P - params["P_0"]) / params["S_D"]
    if T < Tc:
        Q = np.power((Tc - T) / params["Tc_0"], 0.25)
    else:
        Q = 0.0

    # Gibbs
    G = (
        params["Tc_0"] * params["S_D"] * (Q_0 * Q_0 - np.power(Q_0, 6.0) / 3.0)
        - params["S_D"] * (Tc * Q * Q - params["Tc_0"] * np.power(Q, 6.0) / 3.0)
        - T * params["S_D"] * (Q_0 * Q_0 - Q * Q)
        + (P - params["P_0"]) * params["V_D"] * Q_0 * Q_0
    )

    dGdT = params["S_D"] * (Q * Q - Q_0 * Q_0)
    dGdP = -params["V_D"] * (Q * Q - Q_0 * Q_0)

    if Q > 1.0e-12:
        d2GdT2 = -params["S_D"] / (2.0 * params["Tc_0"] * Q * Q)
        d2GdP2 = (
            -params["V_D"]
            * params["V_D"]
            / (2.0 * params["S_D"] * params["Tc_0"] * Q * Q)
        )
        d2GdPdT = params["V_D"] / (2.0 * params["Tc_0"] * Q * Q)
    else:
        d2GdT2 = 0.0
        d2GdP2 = 0.0
        d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, {"Q": Q})


def linear_excesses(pressure, temperature, params):
    """
    A simple linear correction in pressure and temperature.

    The excess Gibbs energy and its derivatives are given by:

    .. math::
        \\mathcal{G} = \\Delta \\mathcal{E} - T \\Delta S + P \\Delta V, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = - \\Delta S, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = \\Delta V, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T \\partial P} = 0

    This form of excess is extremely useful as a first order tweak
    to free energies (especially in solid solutions where data on
    all endmembers may not be available).

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``delta_E``
          - Energy correction (J/mol)
        * - ``delta_S``
          - Entropy correction (J/mol/K)
        * - ``delta_V``
          - Volume correction (m³/mol)

    """
    G = (
        params["delta_E"]
        - (temperature) * params["delta_S"]
        + (pressure) * params["delta_V"]
    )
    dGdT = -params["delta_S"]
    dGdP = params["delta_V"]
    d2GdT2 = 0.0
    d2GdP2 = 0.0
    d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


def bragg_williams_excesses(pressure, temperature, params):
    """
    The Bragg-Williams model is a symmetric solution
    model with an excess configurational entropy term
    determined by the specifics of order-disorder in the
    mineral, multiplied by some empirical factor. Expressions for the
    excess Gibbs free energy can be found in :cite:`HP1996`.

    Excess properties assume that order-disorder processes are
    rapid compared with the timescales of pressure and temperature
    changes (i.e., equilibrium is maintained).

    This may not be reasonable for order-disorder, especially
    for slow or coupled diffusers (Si-Al, for example).
    Properties for minerals that order-disorder slowly
    should be calculated using an explicit solid solution model.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``deltaH``
          - Enthalpy change (J/mol)
        * - ``deltaV``
          - Volume change (m³/mol)
        * - ``Wh``
          - Enthalpy interaction parameter (J/mol)
        * - ``Wv``
          - Volume interaction parameter (m³/mol)
        * - ``n``
          - Related to the number of available sites for ordering
        * - ``factor``
          - An empirical factor

    """

    R = gas_constant
    n = params["n"]
    if params["factor"] > 0.0:
        f = [params["factor"], params["factor"]]
    else:
        f = [1.0, -params["factor"]]

    # Equation A2-2
    def flnarxn(n, Q, f):
        return (
            n
            / (n + 1.0)
            * (
                f[0] * np.log(n * (1.0 - Q))
                + f[1] * np.log(1.0 - Q)
                - f[0] * np.log(1.0 + n * Q)
                - f[1] * np.log(n + Q)
            )
        )

    # Equation A2-4
    # Can be derived from A2-2 noting that
    # delta_H + f*R*T*lnarxn = delta_G + f*R*T*(lnadisord - lnaord)
    def reaction_bragg_williams(Q, delta_H, temperature, n, f, W):
        return delta_H + R * temperature * flnarxn(n, Q, f) + (2.0 * Q - 1.0) * W

    def order_gibbs(pressure, temperature, params):
        W = params["Wh"] + pressure * params["Wv"]
        H_disord = params["deltaH"] + pressure * params["deltaV"]

        # We can use brentq, but don't let the lower bracket = 0
        try:
            Q = opt.brentq(
                reaction_bragg_williams,
                1.0e-12,
                1.0 - 1.0e-12,
                args=(H_disord, temperature, n, f, W),
            )
        except ValueError:
            Q = 0.0

        S = (
            -R
            * (
                f[0]
                * (
                    (1.0 + n * Q) * np.log((1.0 + n * Q) / (n + 1.0))
                    + n * (1.0 - Q) * np.log(n * (1.0 - Q) / (n + 1.0))
                )
                + f[1]
                * (
                    n * (1.0 - Q) * np.log((1.0 - Q) / (n + 1.0))
                    + n * (n + Q) * np.log((n + Q) / (n + 1.0))
                )
            )
            / (n + 1.0)
        )

        G = (1.0 - Q) * H_disord + (1.0 - Q) * Q * W - temperature * S
        return Q, G

    # Calculating partial differentials with respect to P and T
    # are complicated by the fact that Q changes with P and T
    # TODO: Calculate the analytical derivatives of G via the chain rule.
    # For now, just use numerical differentiation.
    dT = 0.1
    dP = 1000.0

    Q, G = order_gibbs(pressure, temperature, params)
    Q, GsubPsubT = order_gibbs(pressure - dP, temperature - dT, params)
    Q, GsubPaddT = order_gibbs(pressure - dP, temperature + dT, params)
    Q, GaddPsubT = order_gibbs(pressure + dP, temperature - dT, params)
    Q, GaddPaddT = order_gibbs(pressure + dP, temperature + dT, params)
    Q, GsubP = order_gibbs(pressure - dP, temperature, params)
    Q, GaddP = order_gibbs(pressure + dP, temperature, params)
    Q, GsubT = order_gibbs(pressure, temperature - dT, params)
    Q, GaddT = order_gibbs(pressure, temperature + dT, params)

    dGdT = (GaddT - GsubT) / (2.0 * dT)
    dGdP = (GaddP - GsubP) / (2.0 * dP)
    d2GdT2 = (GaddT + GsubT - 2.0 * G) / (dT * dT)
    d2GdP2 = (GaddP + GsubP - 2.0 * G) / (dP * dP)
    d2GdPdT = (GaddPaddT - GsubPaddT - GaddPsubT + GsubPsubT) / (4.0 * dT * dP)

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, {"Q": Q})


def magnetic_excesses_chs(pressure, temperature, params):
    """
    This model approximates the excess energy due to magnetic ordering. It
    was originally described in :cite:`CHS1987`. The expressions used
    in this implementation can be found in :cite:`Sundman1991`.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``structural_parameter``
          - A dimensionless parameter related to the crystal structure
        * - ``curie_temperature``
          - A list of length 2: zero pressure Curie temperature (K)
            and pressure dependence (K/Pa)
        * - ``magnetic_moment``
          - A list of length 2: zero pressure magnetic moment (Bohr magnetons)
            and pressure dependence (Bohr magnetons/Pa)
    """

    structural_parameter = params["structural_parameter"]
    curie_temperature = (
        params["curie_temperature"][0] + pressure * params["curie_temperature"][1]
    )
    tau = temperature / curie_temperature
    dtaudT = 1.0 / curie_temperature
    dtaudP = -(temperature * params["curie_temperature"][1]) / (
        curie_temperature * curie_temperature
    )
    d2taudPdT = params["curie_temperature"][1] / (curie_temperature * curie_temperature)
    d2taudP2 = (
        2.0
        * temperature
        * params["curie_temperature"][1]
        * params["curie_temperature"][1]
        / (curie_temperature * curie_temperature * curie_temperature)
    )
    magnetic_moment = (
        params["magnetic_moment"][0] + pressure * params["magnetic_moment"][1]
    )
    dmagnetic_momentdP = params["magnetic_moment"][1]

    A = (518.0 / 1125.0) + (11692.0 / 15975.0) * ((1.0 / structural_parameter) - 1.0)
    if tau < 1:
        f = 1.0 - (1.0 / A) * (
            79.0 / (140.0 * structural_parameter * tau)
            + (474.0 / 497.0)
            * (1.0 / structural_parameter - 1.0)
            * (
                np.power(tau, 3.0) / 6.0
                + np.power(tau, 9.0) / 135.0
                + np.power(tau, 15.0) / 600.0
            )
        )
        dfdtau = -(1.0 / A) * (
            -79.0 / (140.0 * structural_parameter * tau * tau)
            + (474.0 / 497.0)
            * (1.0 / structural_parameter - 1.0)
            * (tau * tau / 2.0 + np.power(tau, 8.0) / 15.0 + np.power(tau, 14.0) / 40.0)
        )
        d2fdtau2 = -(1.0 / A) * (
            2.0 * 79.0 / (140.0 * structural_parameter * np.power(tau, 3.0))
            + (474.0 / 497.0)
            * (1.0 / structural_parameter - 1.0)
            * (
                tau
                + 8.0 * np.power(tau, 7.0) / 15.0
                + 14.0 * np.power(tau, 13.0) / 40.0
            )
        )

    else:
        f = -(1.0 / A) * (
            np.power(tau, -5.0) / 10.0
            + np.power(tau, -15.0) / 315.0
            + np.power(tau, -25.0) / 1500.0
        )
        dfdtau = (1.0 / A) * (
            np.power(tau, -6.0) / 2.0
            + np.power(tau, -16.0) / 21.0
            + np.power(tau, -26.0) / 60.0
        )
        d2fdtau2 = -(1.0 / A) * (
            6.0 * np.power(tau, -7.0) / 2.0
            + 16.0 * np.power(tau, -17.0) / 21.0
            + 26.0 * np.power(tau, -27.0) / 60.0
        )

    dfdT = dfdtau * dtaudT
    d2fdT2 = d2fdtau2 * dtaudT * dtaudT
    dfdP = dfdtau * dtaudP
    d2fdP2 = d2fdtau2 * dtaudP * dtaudP + dfdtau * d2taudP2
    d2fdPdT = d2fdtau2 * dtaudT * dtaudP - dfdtau * d2taudPdT

    G = gas_constant * temperature * np.log(magnetic_moment + 1.0) * f
    dGdT = gas_constant * np.log(magnetic_moment + 1.0) * (f + temperature * dfdT)
    d2GdT2 = (
        gas_constant
        * np.log(magnetic_moment + 1.0)
        * (2.0 * dfdT + temperature * d2fdT2)
    )

    dGdP = (
        gas_constant
        * temperature
        * (
            f * dmagnetic_momentdP / (magnetic_moment + 1.0)
            + dfdP * np.log(magnetic_moment + 1.0)
        )
    )
    d2GdP2 = (
        gas_constant
        * temperature
        * (
            -f * np.power(dmagnetic_momentdP / (magnetic_moment + 1.0), 2.0)
            + 2 * dfdP * dmagnetic_momentdP / (magnetic_moment + 1.0)
            + d2fdP2 * np.log(magnetic_moment + 1.0)
        )
    )
    d2GdPdT = dGdP / temperature + (
        gas_constant * temperature * np.log(magnetic_moment + 1.0) * d2fdPdT
        + gas_constant
        * temperature
        * dmagnetic_momentdP
        / (magnetic_moment + 1.0)
        * dfdT
    )

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


def debye_excesses(pressure, temperature, params):
    """
    Applies an excess contribution based on a Debye model.
    The excess Gibbs energy and its derivatives are given by:

    .. math::

        \\mathcal{G} = F_{\\textrm{Debye}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = - S_{\\textrm{Debye}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = - \\frac{C_{V,\\textrm{Debye}}}{T}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T \\partial P} = 0

    The excess heat capacity tends toward a constant value at high temperature.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``Cv_inf``
          - Heat capacity at infinite temperature (J/mol/K).
        * - ``Theta_0``
          - Debye temperature (K).

    """
    f = params["Cv_inf"] / 3.0 / gas_constant
    theta = params["Theta_0"]

    G = debye.helmholtz_energy(temperature, theta, f)
    dGdT = -debye.entropy(temperature, theta, f)
    dGdP = 0.0
    if temperature > 1.0e-20:
        d2GdT2 = -debye.molar_heat_capacity_v(temperature, theta, f) / temperature
    else:
        d2GdT2 = 0.0
    d2GdP2 = 0.0
    d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


def debye_delta_excesses(pressure, temperature, params):
    """
    Applies an excess contribution based on a Debye model.
    The excess Gibbs energy and its derivatives are given by:

    .. math::

        \\mathcal{G} = - U_{\\textrm{Debye}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = - C_{V,\\textrm{Debye}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = - \\frac{dC_{V,\\textrm{Debye}}}{dT}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T \\partial P} = 0

    The excess entropy tends toward a constant value at high temperature
    and behaves like the heat capacity of a Debye model at finite temperature.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``S_inf``
          - Entropy at infinite temperature (J/mol/K).
        * - ``Theta_0``
          - Einstein temperature (K).

    """
    f = params["S_inf"] / 3.0 / gas_constant
    theta = params["Theta_0"]

    G = -debye.thermal_energy(temperature, theta, f)
    dGdT = -debye.molar_heat_capacity_v(temperature, theta, f)
    dGdP = 0.0
    d2GdT2 = -debye.dmolar_heat_capacity_v_dT(temperature, theta, f)
    d2GdP2 = 0.0
    d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


def einstein_excesses(pressure, temperature, params):
    """
    Applies an excess contribution based an Einstein model.
    The excess Gibbs energy and its derivatives are given by:

    .. math::

        \\mathcal{G} = F_{\\textrm{Einstein}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = - S_{\\textrm{Einstein}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = - \\frac{C_{V,\\textrm{Einstein}}}{T}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T \\partial P} = 0

    The excess heat capacity tends toward a constant value at high temperature.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``Cv_inf``
          - Heat capacity at infinite temperature (J/mol/K).
        * - ``Theta_0``
          - Einstein temperature (K).

    """
    f = params["Cv_inf"] / 3.0 / gas_constant
    theta = params["Theta_0"]

    G = einstein.helmholtz_energy(temperature, theta, f)
    dGdT = -einstein.entropy(temperature, theta, f)
    dGdP = 0.0
    if temperature > 1.0e-20:
        d2GdT2 = -einstein.molar_heat_capacity_v(temperature, theta, f) / temperature
    else:
        d2GdT2 = 0.0
    d2GdP2 = 0.0
    d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


def einstein_delta_excesses(pressure, temperature, params):
    """
    Applies an excess contribution based an Einstein model.
    The excess Gibbs energy and its derivatives are given by:

    .. math::

        \\mathcal{G} = - U_{\\textrm{Einstein}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial T} = - C_{V,\\textrm{Einstein}}, \\\\
        \\frac{\\partial \\mathcal{G}}{\\partial P} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T^2} = - \\frac{dC_{V,\\textrm{Einstein}}}{dT}, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial P^2} = 0, \\\\
        \\frac{\\partial^2 \\mathcal{G}}{\\partial T \\partial P} = 0

    The excess entropy tends toward a
    constant value at high temperature
    and behaves like the heat capacity of an Einstein model
    at finite temperature.

    .. list-table:: Parameters
        :widths: 25 75
        :header-rows: 1

        * - Parameter
          - Description
        * - ``S_inf``
          - Entropy at infinite temperature (J/mol/K).
        * - ``Theta_0``
          - Einstein temperature (K).

    """
    f = params["S_inf"] / 3.0 / gas_constant
    theta = params["Theta_0"]

    G = -einstein.thermal_energy(temperature, theta, f)
    dGdT = -einstein.molar_heat_capacity_v(temperature, theta, f)
    dGdP = 0.0
    d2GdT2 = -einstein.dmolar_heat_capacity_v_dT(temperature, theta, f)
    d2GdP2 = 0.0
    d2GdPdT = 0.0

    excesses = {
        "G": G,
        "dGdT": dGdT,
        "dGdP": dGdP,
        "d2GdT2": d2GdT2,
        "d2GdP2": d2GdP2,
        "d2GdPdT": d2GdPdT,
    }

    return (excesses, None)


modifier_names = OrderedDict()
modifier_names["landau"] = "Landau tricritical model :cite:`Putnis1992`"
modifier_names["landau_slb_2022"] = "Landau tricritical model :cite:`Stixrude2022`"
modifier_names["landau_hp"] = "Landau tricritical model :cite:`Holland1998`"
modifier_names["linear"] = "Linear in P and T"
modifier_names["bragg_williams"] = "Bragg-Williams model"
modifier_names["magnetic_chs"] = "Magnetic ordering :cite:`CHS1987`"
modifier_names["debye"] = "Debye model excess (Helmholtz)"
modifier_names["debye_delta"] = "Debye model excess (internal energy)"
modifier_names["einstein"] = "Einstein model excess (Helmholtz)"
modifier_names["einstein_delta"] = "Einstein model excess (internal energy)"

modifier_functions = OrderedDict()
modifier_functions["landau"] = landau_excesses
modifier_functions["landau_slb_2022"] = landau_slb_2022_excesses
modifier_functions["landau_hp"] = landau_hp_excesses
modifier_functions["linear"] = linear_excesses
modifier_functions["bragg_williams"] = bragg_williams_excesses
modifier_functions["magnetic_chs"] = magnetic_excesses_chs
modifier_functions["debye"] = debye_excesses
modifier_functions["debye_delta"] = debye_delta_excesses
modifier_functions["einstein"] = einstein_excesses
modifier_functions["einstein_delta"] = einstein_delta_excesses


def calculate_property_modifications(mineral):
    """
    Sums the excesses from all the modifiers.

    To calculate thermodynamic properties from the outputs,
    the following functions should be used
    (the _o suffix stands for original value):

    gibbs = gibbs_o + excesses['G']
    S = S_o - excesses['dGdT']
    V = V_o + excesses['dGdP']
    K_T = V / ((V_o / K_T_o) - excesses['d2GdP2'])
    C_p = C_p_o - temperature*excesses['d2GdT2']
    alpha = ((alpha_o*V_o) + excesses['d2GdPdT']) / V

    H = gibbs + temperature*S
    helmholtz = gibbs - pressure*V
    C_v = C_p - V*temperature*alpha*alpha*K_T
    gr = alpha*K_T*V/C_v
    K_S = K_T*C_p/C_v
    """
    excesses = {
        "G": 0.0,
        "dGdT": 0.0,
        "dGdP": 0.0,
        "d2GdT2": 0.0,
        "d2GdP2": 0.0,
        "d2GdPdT": 0.0,
    }
    mineral.property_modifier_properties = []
    for modifier in mineral.property_modifiers:
        if modifier[0] in modifier_functions:
            xs_function = modifier_functions[modifier[0]]
        else:
            raise Exception(
                f"Property modifier label for {mineral.name} ({modifier[0]}) not recognised."
            )

        xs_component, properties = xs_function(
            mineral.pressure, mineral.temperature, modifier[1]
        )
        mineral.property_modifier_properties.append(properties)
        for key in xs_component:
            excesses[key] += xs_component[key]

    return excesses
