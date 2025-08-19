# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np

try:
    import os

    if "NUMBA_DISABLE_JIT" in os.environ and int(os.environ["NUMBA_DISABLE_JIT"]) == 1:
        raise ImportError("NOOOO!")
    from numba import jit
except ImportError:

    def jit(nopython=True):
        def decorator(fn):
            return fn

        return decorator


# Pade coefficients for approximation of the anharmonic contribution to the
# Helmholtz energy. Based on the Taylor expansion of the integral of the
# Debye thermal energy at x_0 = 0.4 divided by x^4,
# followed by a 3-5 Pade approximation multiplied by x^4.
p = np.array([0.0, 0.0, 0.0, 0.0, 0.235256039, -0.00744162069, 23.0722431, 366.827130])
q = np.array([1.0, 4.1784544, 46.00154792, 168.5910059, 568.90929887, 737.13224108])

# First derivatives
p1 = np.polyder(p[::-1])[::-1]
q1 = np.polyder(q[::-1])[::-1]

# Second derivatives
p2 = np.polyder(p1[::-1])[::-1]
q2 = np.polyder(q1[::-1])[::-1]


@jit(nopython=True)
def _helmholtz_pade_pq(t, to_nth_derivative=0):
    """
    Evaluate the Pade approximant to the nondimensional Helmholtz energy
    at a given nondimensional temperature. See the documentation of
    `helmholtz_energy` for details on the model.

    :param t: Nondimensional temperature, defined as
    :math:`T / T_{D}`, where :math:`T_{D}` is the
    Debye temperature.
    :param to_nth_derivative: If 0, return p and q. If 1, also return the
    first derivatives of p and q. If 2, also return the second derivatives
    of p and q.
    :return: A tuple containing the coefficients of the Pade approximant in
    the form ([p, ...], [q, ...]) up to the desired derivative.
    :rtype: tuple of (list, list)
    """
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t
    t5 = t4 * t
    t6 = t5 * t
    t7 = t6 * t

    xp = [p[4] * t4 + p[5] * t5 + p[6] * t6 + p[7] * t7]
    xq = [q[0] + q[1] * t + q[2] * t2 + q[3] * t3 + q[4] * t4 + q[5] * t5]

    if to_nth_derivative > 0:
        xp.append(p1[3] * t3 + p1[4] * t4 + p1[5] * t5 + p1[6] * t6)
        xq.append(q1[0] + q1[1] * t + q1[2] * t2 + q1[3] * t3 + q1[4] * t4)

    if to_nth_derivative > 1:
        xp.append(p2[2] * t2 + p2[3] * t3 + p2[4] * t4 + p2[5] * t5)
        xq.append(q2[0] + q2[1] * t + q2[2] * t2 + q2[3] * t3)

    return xp, xq


@jit(nopython=True)
def _helmholtz_pade(t):
    """
    Evaluate the Pade approximant to the nondimensional Helmholtz energy
    at a given nondimensional temperature. See the documentation of
    `helmholtz_energy` for details on the model.

    :param t: Nondimensional temperature, defined as
    :math:`T / T_{D}`, where :math:`T_{D}` is the
    Debye temperature.
    :return: Nondimensional Helmholtz energy.
    :rtype: float
    """
    xp, xq = _helmholtz_pade_pq(t, to_nth_derivative=0)
    return xp[0] / xq[0]


@jit(nopython=True)
def _dhelmholtzdt_pade(t):
    """
    Evaluate the first derivative of the Pade approximant to the
    nondimensional Helmholtz energy with respect to nondimensional temperature.
    See the documentation of `helmholtz_energy` for details on the model.

    :param t: Nondimensional temperature, defined as
    :math:`T / T_{D}`, where :math:`T_{D}` is the
    Debye temperature.
    :return: float
    """
    xp, xq = _helmholtz_pade_pq(t, to_nth_derivative=1)
    return (xp[1] * xq[0] - xp[0] * xq[1]) / (xq[0] * xq[0])


@jit(nopython=True)
def _d2helmholtzdt2_pade(t):
    """
    Evaluate the second derivative of the Pade approximant to the
    nondimensional Helmholtz energy with respect to nondimensional temperature.
    See the documentation of `helmholtz_energy` for details on the model.

    :param t: Nondimensional temperature, defined as
    :math:`T / T_{D}`, where :math:`T_{D}` is the
    Debye temperature.
    :return: float
    """
    xp, xq = _helmholtz_pade_pq(t, to_nth_derivative=2)
    return (
        xp[2] * xq[0] * xq[0]
        - xp[0] * xq[2] * xq[0]
        - 2.0 * (xp[1] * xq[0] - xp[0] * xq[1]) * xq[1]
    ) / (xq[0] * xq[0] * xq[0])


@jit(nopython=True)
def _nondimensional_helmholtz_energy(T, debye_T):
    """
    Helmholtz free energy of the anharmonic contribution
    at a given temperature.

    This model is based on a 3-5 Pade approximation to the following
    expression:
    :math:`\\int_0^x (E_{D}/3nR) dt / x^{4}`, which is then
    post-multiplied by :math:`x^{4}` to yield the Helmholtz energy.
    The :math:`E_{D}` term inside the integral is the thermal energy of a
    Debye solid per mole of atoms. This expression is chosen because it
    matches the behaviour of the anharmonic contribution to the entropy
    at low and high temperatures - i.e., it is equal to zero at low temperature
    (with all derivatives also equal to zero) and linear at high temperature.
    See Figure 3 in
    Oganov and Dorogokupets (2004; dx.doi.org/10.1088/0953-8984/16/8/018).

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin, which is used to
    nondimensionalise the temperature.
    :type debye_T: float
    :return: Helmholtz energy
    :rtype: float
    """
    t = T / debye_T
    return _helmholtz_pade(t)


@jit(nopython=True)
def _nondimensional_entropy(T, debye_T):
    """
    Entropy of the anharmonic contribution. See the documentation of
    `helmholtz_energy` for details on the model.

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin.
    :type debye_T: float
    :return: Entropy
    :rtype: float
    """
    t = T / debye_T
    return -_dhelmholtzdt_pade(t) / debye_T


@jit(nopython=True)
def _nondimensional_heat_capacity(T, debye_T):
    """
    Heat capacity of the anharmonic contribution. See the documentation of
    `helmholtz_energy` for details on the model.

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin.
    :type debye_T: float
    :return: Heat capacity
    :rtype: float
    """
    t = T / debye_T
    return -t * _d2helmholtzdt2_pade(t) / debye_T


@jit(nopython=True)
def _nondimensional_dhelmholtz_dTheta(T, debye_T):
    """
    Derivative of the anharmonic contribution to the Helmholtz energy
    with respect to the Debye temperature. See the documentation of
    `helmholtz_energy` for details on the model.

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin.
    :type debye_T: float
    :return: Derivative of Helmholtz energy with respect to Debye temperature
    :rtype: float
    """
    t = T / debye_T
    return -_dhelmholtzdt_pade(t) * t / debye_T


@jit(nopython=True)
def _nondimensional_d2helmholtz_dTheta2(T, debye_T):
    """
    Second derivative of the anharmonic contribution to the Helmholtz energy
    with respect to the Debye temperature. See the documentation of
    `helmholtz_energy` for details on the model.

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin.
    :type debye_T: float
    :return: Second derivative of Helmholtz energy with respect to Debye temperature
    :rtype: float
    """
    t = T / debye_T
    return t * (t * _d2helmholtzdt2_pade(t) + 2.0 * _dhelmholtzdt_pade(t)) / debye_T**2


@jit(nopython=True)
def _nondimensional_dentropy_dTheta(T, debye_T):
    """
    Derivative of the anharmonic contribution to the entropy
    with respect to the Debye temperature. See the documentation of
    `helmholtz_energy` for details on the model.

    :param T: Temperature in Kelvin.
    :type T: float
    :param debye_T: Debye temperature in Kelvin.
    :type debye_T: float
    :return: Derivative of entropy with respect to Debye temperature
    :rtype: float
    """
    t = T / debye_T
    return (_d2helmholtzdt2_pade(t) * t + _dhelmholtzdt_pade(t)) / debye_T**2


class AnharmonicDebyePade:
    """
    Class providing methods to compute the anharmonic contribution to the
    Helmholtz free energy, pressure, entropy, isochoric heat capacity,
    isothermal bulk modulus and thermal expansion coefficient
    multiplied by the isothermal bulk modulus.

    The full Helmholtz free energy (relative to the reference isotherm) is
    given by :math:`F = A * (F_a - F_a(T_0))`, where
    :math:`A = a_{anh} * (V/V_0)^{m_{anh}}`, with both :math:`a_{anh}`
    and :math:`m_{anh}` being parameters of the model.
    The term :math:`F_a` is calculated using the 3-5 Pade approximant to
    the function: :math:`\\int_0^x (E_{D}/3nR) dt / x^{4}`, then
    post-multiplied by :math:`x^{4}`.

    The :math:`E_{D}` term inside the integral is the thermal energy of a
    Debye solid per mole of atoms. This expression is chosen because it
    matches the behaviour of the anharmonic contribution to the entropy
    at low and high temperatures - i.e., it is equal to zero at low temperature
    (with all derivatives also equal to zero) and linear at high temperature.
    See Figure 3 in Oganov and Dorogokupets
    (2004; dx.doi.org/10.1088/0953-8984/16/8/018).

    :return: _description_
    :rtype: _type_
    """

    @staticmethod
    def helmholtz_energy(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)
        F_a = _nondimensional_helmholtz_energy(temperature, debye_T)
        F_a0 = _nondimensional_helmholtz_energy(params["T_0"], debye_T)
        return A * (F_a - F_a0)

    @staticmethod
    def entropy(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)
        S_a = _nondimensional_entropy(temperature, debye_T)
        return A * S_a

    @staticmethod
    def heat_capacity_v(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)
        Cv_a = _nondimensional_heat_capacity(temperature, debye_T)
        return A * Cv_a

    @staticmethod
    def pressure(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)
        F_a = _nondimensional_helmholtz_energy(temperature, debye_T)
        F_a0 = _nondimensional_helmholtz_energy(params["T_0"], debye_T)
        F_ad = _nondimensional_dhelmholtz_dTheta(temperature, debye_T)
        F_ad0 = _nondimensional_dhelmholtz_dTheta(params["T_0"], debye_T)
        return -A * (
            (params["m_anh"] / volume) * (F_a - F_a0)
            + (theta_model.dVrel(x, params) / params["V_0"]) * (F_ad - F_ad0)
        )

    @staticmethod
    def isothermal_bulk_modulus(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)

        F_a = _nondimensional_helmholtz_energy(temperature, debye_T)
        F_a0 = _nondimensional_helmholtz_energy(params["T_0"], debye_T)
        F_ad = _nondimensional_dhelmholtz_dTheta(temperature, debye_T)
        F_ad0 = _nondimensional_dhelmholtz_dTheta(params["T_0"], debye_T)
        F_add = _nondimensional_d2helmholtz_dTheta2(temperature, debye_T)
        F_add0 = _nondimensional_d2helmholtz_dTheta2(params["T_0"], debye_T)

        return (
            A
            * volume
            * (
                params["m_anh"] * (params["m_anh"] - 1.0) / volume**2 * (F_a - F_a0)
                + 2
                * params["m_anh"]
                / volume
                * (F_ad - F_ad0)
                * theta_model.dVrel(x, params)
                / params["V_0"]
                + (F_add - F_add0) * (theta_model.dVrel(x, params) / params["V_0"]) ** 2
                + (F_ad - F_ad0) * theta_model.d2dVrel2(x, params) / params["V_0"] ** 2
            )
        )

    @staticmethod
    def dSdV(temperature, volume, params):
        x = volume / params["V_0"]
        A = params["a_anh"] * np.power(x, params["m_anh"])
        theta_model = params["debye_temperature_model"]
        debye_T = theta_model(x, params)

        S_a = _nondimensional_entropy(temperature, debye_T)
        S_ad = _nondimensional_dentropy_dTheta(temperature, debye_T)

        aK_T = A * (
            (params["m_anh"] / volume) * S_a
            + (theta_model.dVrel(x, params) / params["V_0"]) * S_ad
        )

        return aK_T
