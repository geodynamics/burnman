# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
from scipy.special import erf

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


def singleton(cls):
    instances = {}

    def get_instance(*args, **kwargs):
        if cls not in instances:
            instances[cls] = cls(*args, **kwargs)
        return instances[cls]

    return get_instance


class AnharmonicThermalModelBase(object):
    """Abstract base class for anharmonic thermal models."""

    def nondimensional_helmholtz_energy(self, T, debye_T, params):
        """
        Anharmonic contribution to the Helmholtz free energy
        at a given temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_helmholtz_energy() in derived class!"
        )

    def nondimensional_entropy(self, T, debye_T, params):
        """
        Anharmonic contribution to the entropy
        at a given temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_entropy() in derived class!"
        )

    def nondimensional_heat_capacity(self, T, debye_T, params):
        """
        Anharmonic contribution to the heat capacity
        at a given temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_heat_capacity() in derived class!"
        )

    def nondimensional_dhelmholtz_dTheta(self, T, debye_T, params):
        """
        Anharmonic contribution to the derivative of the Helmholtz energy
        with respect to the Debye temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_dhelmholtz_dTheta() in derived class!"
        )

    def nondimensional_d2helmholtz_dTheta2(self, T, debye_T, params):
        """
        Anharmonic contribution to the second derivative of the Helmholtz energy
        with respect to the Debye temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_d2helmholtz_dTheta2() in derived class!"
        )

    def nondimensional_dentropy_dTheta(self, T, debye_T, params):
        """
        Anharmonic contribution to the derivative of the entropy
        with respect to the Debye temperature.
        """
        raise NotImplementedError(
            "need to implement nondimensional_dentropy_dTheta() in derived class!"
        )


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


@singleton
class Pade(object):
    def nondimensional_helmholtz_energy(self, T, debye_T, params=None):
        """
        Anharmonic contribution to the Helmholtz free energy
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
        return _helmholtz_pade(t) * debye_T

    def nondimensional_entropy(self, T, debye_T, params=None):
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
        return -_dhelmholtzdt_pade(t)

    def nondimensional_heat_capacity(self, T, debye_T, params=None):
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
        return -t * _d2helmholtzdt2_pade(t)

    def nondimensional_dhelmholtz_dTheta(self, T, debye_T, params=None):
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
        return _helmholtz_pade(t) - _dhelmholtzdt_pade(t) * t

    def nondimensional_d2helmholtz_dTheta2(self, T, debye_T, params=None):
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
        return t * t * _d2helmholtzdt2_pade(t) / debye_T

    def nondimensional_dentropy_dTheta(self, T, debye_T, params=None):
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
        return _d2helmholtzdt2_pade(t) * t / debye_T


@singleton
class LogNormal(object):
    def _helmholtz(self, x, mu, sigma):
        u = np.log(x)
        u_minus_mu = u - mu
        sigmasqr = sigma**2
        z = u_minus_mu / (sigma * np.sqrt(2))

        term1 = 0.5 * x * ((u_minus_mu - 1) * erf(z) + (u - 1))
        term2 = (
            x * sigma / np.sqrt(2 * np.pi) * np.exp(-(u_minus_mu**2) / (2 * sigmasqr))
        )
        term3 = (
            0.5
            * np.exp(mu + sigmasqr / 2)
            * erf((u_minus_mu - sigmasqr) / (sigma * np.sqrt(2)))
        )

        return term1 + term2 + term3

    def _dhelmholtzdt(self, x, mu, sigma):
        u = np.log(x)
        u_minus_mu = u - mu
        sigmasqr = sigma**2
        z = u_minus_mu / (sigma * np.sqrt(2))

        term1 = 0.5 * (u_minus_mu * erf(z) + u)
        term2 = sigma / np.sqrt(2 * np.pi) * np.exp(-(u_minus_mu**2) / (2 * sigmasqr))

        return term1 + term2

    def _d2helmholtzdt2(self, x, mu, sigma):
        u = np.log(x)
        return 0.5 * (1 + erf((u - mu) / (sigma * np.sqrt(2)))) / x

    def nondimensional_helmholtz_energy(self, T, debye_T, params):
        """
        Anharmonic contribution to the Helmholtz free energy
        at a given temperature.

        This model is based on treating the volumetric heat capacity as
        a scaled version of the cumulative distribution function (CDF) of the
        log-normal distribution, which has the form:

        :math:`F(x) = 0.5 * (1 + erf((\\log(x) - \\mu) / (\\sigma \\sqrt{2})))`
        where :math:`\\mu` and :math:`\\sigma` are the mean and standard deviation
        of the log-normal distribution, respectively, and x is the ratio of the
        temperature to the Debye temperature.

        The nondimensional Helmholtz free energy is then obtained by
        dividing :math:`F(x)` by :math:`x`, integrating twice
        with respect to :math:`x`, and then multiplying by the Debye temperature.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin, which is used to
        nondimensionalise the temperature.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Helmholtz energy
        :rtype: float
        """
        t = T / debye_T
        return self._helmholtz(t, params["mu_anh"], params["sigma_anh"]) * debye_T

    def nondimensional_entropy(self, T, debye_T, params):
        """
        Entropy of the anharmonic contribution. See the documentation of
        `helmholtz_energy` for details on the model.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Entropy
        :rtype: float
        """
        t = T / debye_T
        return -self._dhelmholtzdt(t, params["mu_anh"], params["sigma_anh"])

    def nondimensional_heat_capacity(self, T, debye_T, params):
        """
        Heat capacity of the anharmonic contribution. See the documentation of
        `helmholtz_energy` for details on the model.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Heat capacity
        :rtype: float
        """
        t = T / debye_T
        return -t * self._d2helmholtzdt2(t, params["mu_anh"], params["sigma_anh"])

    def nondimensional_dhelmholtz_dTheta(self, T, debye_T, params):
        """
        Derivative of the anharmonic contribution to the Helmholtz energy
        with respect to the Debye temperature. See the documentation of
        `helmholtz_energy` for details on the model.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Derivative of Helmholtz energy with respect to Debye temperature
        :rtype: float
        """
        t = T / debye_T
        return (
            self._helmholtz(t, params["mu_anh"], params["sigma_anh"])
            - self._dhelmholtzdt(t, params["mu_anh"], params["sigma_anh"]) * t
        )

    def nondimensional_d2helmholtz_dTheta2(self, T, debye_T, params):
        """
        Second derivative of the anharmonic contribution to the Helmholtz energy
        with respect to the Debye temperature. See the documentation of
        `helmholtz_energy` for details on the model.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Second derivative of Helmholtz energy with respect to Debye temperature
        :rtype: float
        """
        t = T / debye_T
        return (
            t
            * t
            * self._d2helmholtzdt2(t, params["mu_anh"], params["sigma_anh"])
            / debye_T
        )

    def nondimensional_dentropy_dTheta(self, T, debye_T, params):
        """
        Derivative of the anharmonic contribution to the entropy
        with respect to the Debye temperature. See the documentation of
        `helmholtz_energy` for details on the model.

        :param T: Temperature in Kelvin.
        :type T: float
        :param debye_T: Debye temperature in Kelvin.
        :type debye_T: float
        :param params: Parameters for the anharmonic model, which should contain
            - mu_anh: Mean of the log-normal distribution
            - sigma_anh: Standard deviation of the log-normal distribution
        :type params: dict
        :return: Derivative of entropy with respect to Debye temperature
        :rtype: float
        """
        t = T / debye_T
        return (
            self._d2helmholtzdt2(t, params["mu_anh"], params["sigma_anh"]) * t / debye_T
        )
