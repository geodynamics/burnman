# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import importlib
import numpy as np
from ..utils.chemistry import process_solution_chemistry
from .. import constants
import warnings
from .material import material_property, cached_property
from dataclasses import make_dataclass

Interaction = make_dataclass(
    "Interaction", ["inds", "expts", "f_r", "m_jr", "f_rs", "m_jrs"]
)

try:
    ag = importlib.import_module("autograd")
except ImportError as err:
    print(
        f"Warning: {err}. "
        "For full functionality of BurnMan, please install autograd."
    )


def _ideal_activities_fct(
    molar_fractions,
    endmember_noccupancies,
    n_endmembers,
    n_occupancies,
    site_multiplicities,
    endmember_configurational_entropies,
):
    site_noccupancies = np.dot(molar_fractions, endmember_noccupancies)
    site_multiplicities = np.einsum("i, ij", molar_fractions, site_multiplicities)
    site_occupancies = site_noccupancies * inverseish(site_multiplicities)

    a = np.power(site_occupancies, endmember_noccupancies).prod(-1)
    normalisation_constants = np.exp(
        endmember_configurational_entropies / constants.gas_constant
    )
    return normalisation_constants * a


def _non_ideal_hessian_fct(phi, molar_fractions, n_endmembers, alpha, W):
    q = np.eye(n_endmembers) - phi * np.ones((n_endmembers, n_endmembers))
    sum_pa = np.dot(molar_fractions, alpha)
    hess = np.einsum("m, i, ij, jk, mk->im", -alpha / sum_pa, -alpha, q, W, q)
    hess += hess.T
    return hess


def _non_ideal_interactions_fct(phi, molar_fractions, n_endmembers, alpha, W):
    # -sum(sum(qi.qj.Wij*)
    # equation (2) of Holland and Powell 2003
    q = np.eye(n_endmembers) - phi * np.ones((n_endmembers, n_endmembers))
    # The following are equivalent to
    # np.einsum('i, ij, jk, ik->i', -self.alphas, q, self.Wx, q)
    Wint = -alpha * (q.dot(W) * q).sum(-1)
    return Wint


def _non_ideal_hessian_subreg(p, n_endmembers, Wijk):
    Id = np.identity(n_endmembers)
    IIp = np.einsum("il, jm, k->ijklm", Id, Id, p)
    Ipp = np.einsum("il, j, k->ijkl", Id, p, p)
    ppp = np.einsum("i, j, k->ijk", p, p, p)

    A = (
        IIp
        + np.transpose(IIp, axes=[0, 2, 1, 3, 4])
        + np.transpose(IIp, axes=[1, 0, 2, 3, 4])
        + np.transpose(IIp, axes=[1, 2, 0, 3, 4])
        + np.transpose(IIp, axes=[2, 1, 0, 3, 4])
        + np.transpose(IIp, axes=[2, 0, 1, 3, 4])
    )
    B = 2.0 * (
        Ipp
        + np.transpose(Ipp, axes=[1, 0, 2, 3])
        + np.transpose(Ipp, axes=[2, 1, 0, 3])
    )

    Asum = (
        A - B[:, :, :, :, None] - B[:, :, :, None, :] + 6.0 * ppp[:, :, :, None, None]
    )
    hess = np.einsum("ijklm, ijk->lm", Asum, Wijk)
    return hess


def _non_ideal_interactions_subreg(p, n_endmembers, Wijk):
    Aijkl = np.einsum("li, j, k->ijkl", np.identity(n_endmembers), p, p)
    ppp = np.einsum("i, j, k->ijk", p, p, p)

    Asum = (
        Aijkl
        + np.transpose(Aijkl, axes=[1, 0, 2, 3])
        + np.transpose(Aijkl, axes=[1, 2, 0, 3])
        - 2 * ppp[:, :, :, None]
    )

    Wint = np.einsum("ijk, ijkl->l", Wijk, Asum)
    return Wint


def logish(x, eps=1.0e-7):
    """
    2nd order series expansion of log(x) about eps:
    log(eps) - sum_k=1^infty (f_eps)^k / k
    Prevents infinities at x=0
    """
    f_eps = 1.0 - x / eps
    mask = x > eps
    ln = np.where(x <= eps, np.log(eps) - f_eps - f_eps * f_eps / 2.0, 0.0)
    ln[mask] = np.log(x[mask])
    return ln


def inverseish(x, eps=1.0e-5):
    """
    1st order series expansion of 1/x about eps: 2/eps - x/eps/eps
    Prevents infinities at x=0
    """
    mask = x > eps
    oneoverx = np.where(x <= eps, 2.0 / eps - x / eps / eps, 0.0)
    oneoverx[mask] = 1.0 / x[mask]
    return oneoverx


def dpdn(molar_amounts, n, ones, eye):
    """
    The partial derivative of endmember proportions
    with respect to endmember amounts.

    :param molar_amounts: molar amounts of independent endmembers
    :type molar_fractions: 1D numpy array
    :param n: sum of the endmember amounts (usually equal to one)
    :type n: float
    :param ones: a vector of ones of length equal to the number of endmembers
    :type ones: 1D numpy array
    :param eye: the identity matrix of size equal to the number of endmembers
    :type eye: 2D numpy array
    """
    return (eye - np.einsum("k, m->km", molar_amounts, ones)) / n


def d2pdndn(molar_amounts, nsqr, ones, eyeones):
    """
    The second partial derivative of endmember proportions
    with respect to endmember amounts.

    :param molar_amounts: molar amounts of independent endmembers
    :type molar_fractions: 1D numpy array
    :param nsqr: square of the sum of the endmember amounts (usually equal to one)
    :type n: float
    :param ones: a vector of ones of length equal to the number of endmembers
    :type ones: 1D numpy array
    :param eyeones: delta_ij 1_k + delta_ik 1_j
    :type eyeones: 3D numpy array
    """

    return (2.0 * np.einsum("n, m, k->knm", ones, ones, molar_amounts) - eyeones) / nsqr


class SolutionModel(object):
    """
    This is the base class for a solution model,  intended for use
    in defining solutions and performing thermodynamic calculations
    on them.  All minerals of type :class:`burnman.Solution` use
    a solution model for defining how the endmembers in the solution
    interact.

    A user wanting a new solution model should define the functions included
    in the base class. All of the functions in the base class return zero,
    so if the user-defined solution model does not implement them,
    they essentially have no effect, and the Gibbs free energy and molar
    volume of a solution will be equal to the weighted arithmetic
    averages of the different endmember values.
    """

    def __init__(self):
        """
        Does nothing.
        """
        pass

    def excess_gibbs_free_energy(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Gibbs free energy of the solution.
        The base class implementation assumes that the excess gibbs
        free energy is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess Gibbs energy.
        :rtype: float
        """
        return np.dot(
            np.array(molar_fractions),
            self.excess_partial_gibbs_free_energies(
                pressure, temperature, molar_fractions
            ),
        )

    def excess_volume(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess volume of the solution.
        The base class implementation assumes that the excess volume is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess volume of the solution.
        :rtype: float
        """
        return np.dot(
            molar_fractions,
            self.excess_partial_volumes(pressure, temperature, molar_fractions),
        )

    def excess_entropy(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy of the solution.
        The base class implementation assumes that the excess entropy is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess entropy of the solution.
        :rtype: float
        """
        return np.dot(
            molar_fractions,
            self.excess_partial_entropies(pressure, temperature, molar_fractions),
        )

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess enthalpy of the solution.
        The base class implementation assumes that the excess enthalpy is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess enthalpy of the solution.
        :rtype: float
        """
        return self.excess_gibbs_free_energy(
            pressure, temperature, molar_fractions
        ) + temperature * self.excess_entropy(pressure, temperature, molar_fractions)

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        """
        Given a list of molar fractions of different phases,
        compute the excess Gibbs free energy for each endmember
        of the solution.
        The base class implementation assumes that the excess gibbs
        free energy is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess partial Gibbs free energy of each endmember.
        :rtype: numpy.array
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy for each endmember of the solution.
        The base class implementation assumes that the excess entropy
        is zero (true for mechanical solutions).

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess partial entropy of each endmember.
        :rtype: numpy.array
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess volume for each endmember of the solution.
        The base class implementation assumes that the excess volume
        is zero.

        :param pressure: Pressure at which to evaluate the solution model [Pa].
        :type pressure: float
        :param temperature: Temperature at which to evaluate the solution model [K].
        :type temperature: float
        :param molar_fractions: List of molar fractions of the
            different independent endmembers in the solution model.
        :type molar_fractions: list of floats

        :returns: The excess partial volume of each endmember.
        :rtype: numpy.array
        """
        return np.zeros_like(np.array(molar_fractions))

    def Cp_excess(self):
        """
        Returns the excess heat capacity of the solution model
        at its current state
        """
        return 0.0

    def alphaV_excess(self):
        """
        Returns the excess alpha*V of the solution model
        at its current state
        """
        return 0.0

    def VoverKT_excess(self):
        """
        Returns the excess V/K_T of the solution model
        at its current state
        """
        return 0.0


class MechanicalSolution(SolutionModel):
    """
    An extremely simple class representing a mechanical solution model.
    A mechanical solution experiences no interaction between endmembers.
    Therefore, unlike ideal solutions there is no entropy of mixing;
    the total gibbs free energy of the solution is equal to the
    dot product of the molar gibbs free energies and molar fractions
    of the constituent materials.
    """

    def __init__(self, endmembers):
        self.endmembers = endmembers
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

    def excess_gibbs_free_energy(self, pressure, temperature, molar_fractions):
        return 0.0

    def excess_volume(self, pressure, temperature, molar_fractions):
        return 0.0

    def excess_entropy(self, pressure, temperature, molar_fractions):
        return 0.0

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        return 0.0

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        return np.zeros_like(molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def activities(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)


class IdealSolution(SolutionModel):
    """
    A class representing an ideal solution model.
    Calculates the excess gibbs free energy and entropy due to configurational
    entropy. Excess internal energy and volume are equal to zero.

    The multiplicity of each type of site in the structure is allowed to
    change linearly as a function of endmember proportions. This class
    is therefore equivalent to the entropic part of
    a Temkin-type model :cite:`Temkin1945`.
    """

    def __init__(self, endmembers):
        self.endmembers = endmembers
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

        # Process solution chemistry
        process_solution_chemistry(self)

        self._calculate_endmember_configurational_entropies()

    def _calculate_endmember_configurational_entropies(self):
        S_conf = -(
            constants.gas_constant
            * (
                self.endmember_noccupancies
                * (
                    logish(self.endmember_noccupancies)
                    - logish(self.site_multiplicities)
                )
            ).sum(-1)
        )
        self.endmember_configurational_entropies = S_conf

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        return self._ideal_excess_partial_gibbs(temperature, molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        return self._ideal_excess_partial_entropies(temperature, molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return np.zeros((self.n_endmembers))

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return -temperature * hess_S

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return hess_S

    def volume_hessian(self, pressure, temperature, molar_fractions):
        return np.zeros((len(molar_fractions), len(molar_fractions)))

    def _configurational_entropy(self, molar_fractions):
        site_noccupancies = np.einsum(
            "i, ij", molar_fractions, self.endmember_noccupancies
        )
        site_multiplicities = np.einsum(
            "i, ij", molar_fractions, self.site_multiplicities
        )
        site_occupancies = site_noccupancies * inverseish(site_multiplicities)
        conf_entropy = -(
            constants.gas_constant
            * (site_noccupancies * logish(site_occupancies)).sum(-1)
        )
        return conf_entropy

    def _ideal_excess_partial_gibbs(self, temperature, molar_fractions):
        return -(
            temperature
            * self._ideal_excess_partial_entropies(temperature, molar_fractions)
        )

    def _ideal_excess_partial_entropies(self, temperature, molar_fractions):
        return -(constants.gas_constant * self._log_ideal_activities(molar_fractions))

    def _ideal_entropy_hessian(self, temperature, molar_fractions):
        hessian = -(
            constants.gas_constant
            * self._log_ideal_activity_derivatives(molar_fractions)
        )
        return hessian

    def _log_ideal_activities(self, molar_fractions):
        site_noccupancies = np.einsum(
            "i, ij", molar_fractions, self.endmember_noccupancies
        )
        site_multiplicities = np.einsum(
            "i, ij", molar_fractions, self.site_multiplicities
        )

        lna = np.einsum(
            "ij, j->i",
            self.endmember_noccupancies,
            logish(site_noccupancies) - logish(site_multiplicities),
        )
        normalisation_constants = (
            self.endmember_configurational_entropies / constants.gas_constant
        )
        return lna + normalisation_constants

    def _log_ideal_activity_derivatives(self, molar_fractions):
        site_noccupancies = np.einsum(
            "i, ij", molar_fractions, self.endmember_noccupancies
        )
        site_multiplicities = np.einsum(
            "i, ij", molar_fractions, self.site_multiplicities
        )

        dlnadp = np.einsum(
            "pj, qj, j->pq",
            self.endmember_noccupancies,
            self.endmember_noccupancies,
            inverseish(site_noccupancies),
        ) - np.einsum(
            "pj, qj, j->pq",
            self.endmember_noccupancies,
            self.site_multiplicities,
            inverseish(site_multiplicities),
        )

        return dlnadp

    def _ideal_activities(self, molar_fractions):
        return _ideal_activities_fct(
            molar_fractions,
            self.endmember_noccupancies,
            self.n_endmembers,
            self.n_occupancies,
            self.site_multiplicities,
            self.endmember_configurational_entropies,
        )

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def activities(self, pressure, temperature, molar_fractions):
        return self._ideal_activities(molar_fractions)

    @cached_property
    def ones(self):
        """
        A vector of ones with length equal to the number of endmembers
        :return: ones
        :rtype: 1D numpy array
        """
        return np.ones(self.n_endmembers)

    @cached_property
    def eye(self):
        """
        An identity matrix with size equal to the number of endmembers
        :return: eye
        :rtype: 2D numpy array
        """
        return np.eye(self.n_endmembers)

    @cached_property
    def eyeones(self):
        """
        A convenience function consisting of two concatenations
        of an identity matrix and ones vector with size
        equal to the number of endmembers.
        :return: delta_ij 1_k + delta_ik 1_j
        :rtype: 3D numpy array
        """
        return np.einsum("km, n->kmn", self.eye, self.ones) + np.einsum(
            "kn, m->kmn", self.eye, self.ones
        )


class AsymmetricRegularSolution(IdealSolution):
    """
    Solution model implementing the asymmetric regular solution model
    formulation as described in :cite:`HP2003`.

    The excess nonconfigurational Gibbs energy is given by the
    expression:

    .. math::
        \\mathcal{G}_{\\textrm{excess}} = \\alpha^T p (\\phi^T W \\phi)

    :math:`\\alpha` is a vector of van Laar parameters governing asymmetry
    in the excess properties.

    .. math::
        \\phi_i = \\frac{\\alpha_i p_i}{\\sum_{k=1}^{n} \\alpha_k p_k},
        W_{ij} = \\frac{2 w_{ij}}{\\alpha_i + \\alpha_j} \\textrm{for i<j}
    """

    def __init__(
        self,
        endmembers,
        alphas,
        energy_interaction,
        volume_interaction=None,
        entropy_interaction=None,
    ):
        self.n_endmembers = len(endmembers)

        # Create array of van Laar parameters
        self.alphas = np.array(alphas)

        # Create 2D arrays of interaction parameters
        self.We = np.triu(2.0 / (self.alphas[:, np.newaxis] + self.alphas), 1)
        self.We[np.triu_indices(self.n_endmembers, 1)] *= np.array(
            [i for row in energy_interaction for i in row]
        )

        if entropy_interaction is not None:
            self.Ws = np.triu(2.0 / (self.alphas[:, np.newaxis] + self.alphas), 1)
            self.Ws[np.triu_indices(self.n_endmembers, 1)] *= np.array(
                [i for row in entropy_interaction for i in row]
            )
        else:
            self.Ws = np.zeros((self.n_endmembers, self.n_endmembers))

        if volume_interaction is not None:
            self.Wv = np.triu(2.0 / (self.alphas[:, np.newaxis] + self.alphas), 1)
            self.Wv[np.triu_indices(self.n_endmembers, 1)] *= np.array(
                [i for row in volume_interaction for i in row]
            )
        else:
            self.Wv = np.zeros((self.n_endmembers, self.n_endmembers))

        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def _phi(self, molar_fractions):
        phi = self.alphas * molar_fractions
        phi = np.divide(phi, np.sum(phi))
        return phi

    def _non_ideal_interactions(self, W, molar_fractions):
        # -sum(sum(qi.qj.Wij*)
        # equation (2) of Holland and Powell 2003
        phi = self._phi(molar_fractions)
        return _non_ideal_interactions_fct(
            phi, np.array(molar_fractions), self.n_endmembers, self.alphas, W
        )

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions):
        Eint = self._non_ideal_interactions(self.We, molar_fractions)
        Sint = self._non_ideal_interactions(self.Ws, molar_fractions)
        Vint = self._non_ideal_interactions(self.Wv, molar_fractions)
        return Eint - temperature * Sint + pressure * Vint

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions
        )
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions
        )
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_interactions(self.Ws, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return self._non_ideal_interactions(self.Wv, molar_fractions)

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        phi = self._phi(molar_fractions)
        nonideal_gibbs_hessian = _non_ideal_hessian_fct(
            phi,
            molar_fractions,
            self.n_endmembers,
            self.alphas,
            self.We - temperature * self.Ws + pressure * self.Wv,
        )

        return nonideal_gibbs_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        phi = self._phi(molar_fractions)
        nonideal_entropy_hessian = _non_ideal_hessian_fct(
            phi, molar_fractions, self.n_endmembers, self.alphas, self.Ws
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def volume_hessian(self, pressure, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        return _non_ideal_hessian_fct(
            phi, molar_fractions, self.n_endmembers, self.alphas, self.Wv
        )

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.0e-10:
            return np.exp(
                self._non_ideal_excess_partial_gibbs(
                    pressure, temperature, molar_fractions
                )
                / (constants.gas_constant * temperature)
            )
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(
            self, pressure, temperature, molar_fractions
        ) * self.activity_coefficients(pressure, temperature, molar_fractions)


class SymmetricRegularSolution(AsymmetricRegularSolution):
    """
    Solution model implementing the symmetric regular solution model.
    This is a special case of the
    :class:`burnman.solutionmodel.AsymmetricRegularSolution` class.
    """

    def __init__(
        self,
        endmembers,
        energy_interaction,
        volume_interaction=None,
        entropy_interaction=None,
    ):
        alphas = np.ones(len(endmembers))
        AsymmetricRegularSolution.__init__(
            self,
            endmembers,
            alphas,
            energy_interaction,
            volume_interaction,
            entropy_interaction,
        )


class SubregularSolution(IdealSolution):
    """
    Solution model implementing the subregular solution model formulation
    as described in :cite:`HW1989`. The excess nonconfigurational
    Gibbs energy is given by the expression:

    .. math::
        \\mathcal{G}_{\\textrm{excess}} = \\sum_i \\sum_{j > i} (p_i p_j^2
        W_{ij} + p_j p_i^2 W_{ji} + \\sum_{k > j > i} p_i p_j p_k W_{ijk})

    Interaction parameters are inserted into a 3D interaction matrix during
    initialization to make use of numpy vector algebra.

    :param endmembers: A list of all the independent endmembers in the solution.
        The first item of each list gives the Mineral object corresponding
        to the endmember. The second item gives the site-species formula.
    :type endmembers: list of lists
    :param energy_interaction: The binary endmember interaction energies.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).
    :type energy_interaction: list of list of lists
    :param volume_interaction: The binary endmember interaction volumes.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).
    :type volume_interaction: list of list of lists
    :param entropy_interaction: The binary endmember interaction entropies.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).
    :type entropy_interaction: list of list of lists
    :param energy_ternary_terms: The ternary interaction energies.
        Each list should contain four entries:
        the indices i, j, k and the value of the interaction.
    :type energy_ternary_terms: list of lists
    :param volume_ternary_terms: The ternary interaction volumes.
        Each list should contain four entries:
        the indices i, j, k and the value of the interaction.
    :type volume_ternary_terms: list of lists
    :param entropy_ternary_terms: The ternary interaction entropies.
        Each list should contain four entries:
        the indices i, j, k and the value of the interaction.
    :type entropy_ternary_terms: list of lists
    """

    def __init__(
        self,
        endmembers,
        energy_interaction,
        volume_interaction=None,
        entropy_interaction=None,
        energy_ternary_terms=None,
        volume_ternary_terms=None,
        entropy_ternary_terms=None,
    ):
        """
        Initialization function for the SubregularSolution class.
        """

        self.n_endmembers = len(endmembers)

        # Create 3D arrays of interaction parameters
        self.Wijke = np.zeros(
            shape=(self.n_endmembers, self.n_endmembers, self.n_endmembers)
        )
        self.Wijks = np.zeros_like(self.Wijke)
        self.Wijkv = np.zeros_like(self.Wijke)

        # setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i + 1, self.n_endmembers):
                w0 = energy_interaction[i][j - i - 1][0] / 2.0
                w1 = energy_interaction[i][j - i - 1][1] / 2.0
                self.Wijke[:, i, j] += w0
                self.Wijke[:, j, i] += w1

                self.Wijke[i, j, j] += w0
                self.Wijke[j, i, i] += w1

                self.Wijke[i, j, i] -= w0
                self.Wijke[j, i, j] -= w1

        if energy_ternary_terms is not None:
            for i, j, k, v in energy_ternary_terms:
                self.Wijke[i, j, k] += v

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    w0 = entropy_interaction[i][j - i - 1][0] / 2.0
                    w1 = entropy_interaction[i][j - i - 1][1] / 2.0
                    self.Wijks[:, i, j] += w0
                    self.Wijks[:, j, i] += w1

                    self.Wijks[i, j, j] += w0
                    self.Wijks[j, i, i] += w1

                    self.Wijks[i, j, i] -= w0
                    self.Wijks[j, i, j] -= w1

        if entropy_ternary_terms is not None:
            for i, j, k, v in entropy_ternary_terms:
                self.Wijks[i, j, k] += v

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    w0 = volume_interaction[i][j - i - 1][0] / 2.0
                    w1 = volume_interaction[i][j - i - 1][1] / 2.0
                    self.Wijkv[:, i, j] += w0
                    self.Wijkv[:, j, i] += w1

                    self.Wijkv[i, j, j] += w0
                    self.Wijkv[j, i, i] += w1

                    self.Wijkv[i, j, i] -= w0
                    self.Wijkv[j, i, j] -= w1

        if volume_ternary_terms is not None:
            for i, j, k, v in volume_ternary_terms:
                self.Wijkv[i, j, k] += v

        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def _non_ideal_function(self, Wijk, molar_fractions):
        n = len(molar_fractions)
        return _non_ideal_interactions_subreg(molar_fractions, n, Wijk)

    def _non_ideal_interactions(self, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        Eint = self._non_ideal_function(self.Wijke, molar_fractions)
        Sint = self._non_ideal_function(self.Wijks, molar_fractions)
        Vint = self._non_ideal_function(self.Wijkv, molar_fractions)
        return Eint, Sint, Vint

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions):
        Eint, Sint, Vint = self._non_ideal_interactions(molar_fractions)
        return Eint - temperature * Sint + pressure * Vint

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions
        )
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions
        )
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_function(self.Wijks, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        non_ideal_volumes = self._non_ideal_function(self.Wijkv, molar_fractions)
        return non_ideal_volumes

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_gibbs_hessian = _non_ideal_hessian_subreg(
            molar_fractions,
            n,
            self.Wijke - temperature * self.Wijks + pressure * self.Wijkv,
        )

        return nonideal_gibbs_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_entropy_hessian = _non_ideal_hessian_subreg(
            molar_fractions, n, self.Wijks
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def volume_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        return _non_ideal_hessian_subreg(molar_fractions, n, self.Wijkv)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.0e-10:
            return np.exp(
                self._non_ideal_excess_partial_gibbs(
                    pressure, temperature, molar_fractions
                )
                / (constants.gas_constant * temperature)
            )
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(
            self, pressure, temperature, molar_fractions
        ) * self.activity_coefficients(pressure, temperature, molar_fractions)


class FunctionSolution(IdealSolution):
    """
    Solution model implementing a generalized solution model.
    The extensive excess nonconfigurational Gibbs energy is
    provided as a function by the user.

    Derivatives are calculated using the autograd module,
    and so the user-defined excess Gibbs energy function
    should be defined using autograd-friendly expressions.

    :param endmembers: A list of all the independent endmembers in the solution.
        The first item of each list gives the Mineral object corresponding
        to the endmember. The second item gives the site-species formula.
    :type endmembers: list of lists

    :param excess_gibbs_function: The nonconfigurational Gibbs energy function
        with arguments pressure, temperature and molar_amounts, in that order.
        Note that the function must be extensive; if the molar amounts
        are doubled, the Gibbs energy must also double.
    :type excess_gibbs_function: function
    """

    def __init__(self, endmembers, excess_gibbs_function):
        """
        Initialization function for the GeneralSolution class.
        """

        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

        self.n_endmembers = len(endmembers)
        self._excess_gibbs_function = excess_gibbs_function

        self._non_ideal_excess_partial_gibbs = ag.jacobian(
            excess_gibbs_function, argnum=2
        )

        def partial_entropies(pressure, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return -ag.jacobian(self._non_ideal_excess_partial_gibbs, argnum=1)(
                    pressure, temperature, molar_amounts
                )

        self._non_ideal_excess_partial_entropies = partial_entropies

        def partial_volumes(pressure, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return ag.jacobian(self._non_ideal_excess_partial_gibbs, argnum=0)(
                    pressure, temperature, molar_amounts
                )

        self.excess_partial_volumes = partial_volumes

        self._non_ideal_gibbs_hessian = ag.jacobian(
            self._non_ideal_excess_partial_gibbs, argnum=2
        )

        def entropy_hess(pressure, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return ag.jacobian(partial_entropies, argnum=2)(
                    pressure, temperature, molar_amounts
                )

        self._non_ideal_entropy_hessian = entropy_hess

        def volume_hess(pressure, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return ag.jacobian(partial_volumes, argnum=2)(
                    pressure, temperature, molar_amounts
                )

        self.volume_hessian = volume_hess

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions
        )
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions
        )
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_excess_partial_entropies(
            pressure, temperature, molar_fractions
        )
        return ideal_entropies + non_ideal_entropies

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_gibbs_hessian = self._non_ideal_gibbs_hessian(
            pressure, temperature, molar_fractions
        )

        return nonideal_gibbs_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_entropy_hessian = self._non_ideal_entropy_hessian(
            pressure, temperature, molar_fractions
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.0e-10:
            return np.exp(
                self._non_ideal_excess_partial_gibbs(
                    pressure, temperature, molar_fractions
                )
                / (constants.gas_constant * temperature)
            )
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(
            self, pressure, temperature, molar_fractions
        ) * self.activity_coefficients(pressure, temperature, molar_fractions)


class PolynomialSolution(IdealSolution):
    """
    Solution model implementing a general polynomial solution model.

    :param endmembers: A list of all the independent endmembers in the solution.
        The first item of each list gives the Mineral object corresponding
        to the endmember. The second item gives the site-species formula.
    :type endmembers: list of lists
    :param ESV_interactions: A list containing lists where the first three elements are
        energy, entropy and volume interactions and the rest of the elements
        are indices of the transformed endmembers to which those
        interactions correspond, and the proportion exponents.
        For example, [2., 0., 0., 0, 4, 1, 1] would correspond to an interaction
        of 2*p'[0]^4*p'[1]^1.
    :type ESV_interactions: list of lists
    :param interaction_endmembers: A list of minerals involved in the interaction terms.
    :type interaction_endmembers: list of :class:`burnman.Mineral` objects
    :param endmember_coefficients_and_interactions: list of lists
        A list containing lists where the first n elements are
        coefficients for each of the interaction_endmembers and the
        rest of the elements are indices of the transformed
        endmembers to which those interactions correspond.
        For example, [1., 0., -1., 0, 4, 1, 1] would correspond to an interaction
        of (mbr[0].gibbs - mbr[2].gibbs)*p'[0]^4*p'[1]^1.
    :type endmember_coefficients_and_interactions: list of lists
    :param transformation_matrix: The interactions for a given solution may
        be most compactly expressed not as a polynomial function of the
        proportions of the endmembers, but a polynomial function of a
        linearly transformed set. This parameter is a square numpy array A,
        where p'i = A_ij p_j
    :type transformation_matrix: 2D numpy array
    """

    def __init__(
        self,
        endmembers,
        ESV_interactions=None,
        interaction_endmembers=[],
        endmember_coefficients_and_interactions=None,
        transformation_matrix=None,
    ):
        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

        self.n_endmembers = len(endmembers)
        self.endmembers = endmembers

        self.W_ESV = None
        self.c_ESV = None
        self.n_W_ESV = 0
        if ESV_interactions is not None:
            self.W_ESV, self.c_ESV = self._make_interaction_arrays(ESV_interactions, 3)
            self.n_W_ESV = len(ESV_interactions)

        self.n_interaction_endmembers = len(interaction_endmembers)
        self.interaction_endmembers = interaction_endmembers

        self.W_mbr = None
        self.c_mbr = None
        self.n_W_mbr = 0
        if self.n_interaction_endmembers > 0:
            self.W_mbr, self.c_mbr = self._make_interaction_arrays(
                endmember_coefficients_and_interactions, self.n_interaction_endmembers
            )
            self.n_W_mbr = len(endmember_coefficients_and_interactions)

        if transformation_matrix is None:
            self.dqdp = np.eye(self.n_endmembers)
        else:
            self.dqdp = transformation_matrix
        self.n_transformed_endmembers = len(self.dqdp)
        self.reset()

    def reset(self):
        """
        Resets all cached material properties.
        It is typically not required for the user to call this function.
        """
        self._cached = {}

    def set_composition(self, molar_fractions):
        """
        Resets all cached material properties.
        It is typically not required for the user to call this function.
        """
        self.reset()
        self.molar_fractions = molar_fractions
        self.trans_fractions = np.einsum("ij, j->i", self.dqdp, molar_fractions)

    def set_state(self, pressure, temperature):
        """
        Sets the states for the interaction endmembers.
        It is typically not required for the user to call this function.
        """
        for mbr in self.interaction_endmembers:
            mbr.set_state(pressure, temperature)

    def _make_prefactors_and_exponents(self, a):
        """
        This hidden method makes an Interaction object
        containing convenience arrays to facilitate rapid
        calculation of the first and second compositional
        derivatives of the Gibbs, entropy and volume
        excesses.

        :param a: A list of polynomial indices and exponents.
        :type a: list
        :return: An interaction object containing the
            indices (inds), exponents (expts),
            first derivative prefactors (j_r) and exponents (m_jr) and
            second derivative prefactors (f_rs) and exponents (m_jrs).
        :rtype: Interaction dataclass instance
        """
        # assign exponent matrix for first derivatives
        if len(a) % 2 != 0:
            raise Exception(
                "Interaction must have paired indices " f"and exponents, currently {a}"
            )
        indices = a[::2]
        exponents = np.array(a[1::2])
        n_indices = len(indices)
        m_jr = np.zeros((n_indices, self.n_endmembers))
        m_jr[:, :] = exponents[:, np.newaxis]
        m_jr[list(range(len(m_jr))), indices] -= np.ones(len(indices))

        # assign factors for first derivatives
        f_jr = np.zeros((n_indices, self.n_endmembers))
        f_jr[:, indices] = 1.0
        f_jr[range(len(m_jr)), indices] = exponents
        f_r = np.prod(f_jr, axis=0)

        # assign exponent matrix for second derivatives
        m_jrs = np.zeros((n_indices, self.n_endmembers, self.n_endmembers))
        m_jrs[:, :, :] = exponents[:, np.newaxis, np.newaxis]

        m_jrs -= (
            np.einsum(
                "ij, k->ijk", np.eye(self.n_endmembers), np.ones(self.n_endmembers)
            )
            + np.einsum(
                "ik, j->ijk", np.eye(self.n_endmembers), np.ones(self.n_endmembers)
            )
        )[indices]

        # assign factors for second derivatives
        f_jrs = np.zeros((n_indices, self.n_endmembers, self.n_endmembers))
        f_jrs[np.ix_(range(n_indices), indices, indices)] = 1.0

        ival = np.array([indices, exponents]).T
        for k, (i, iexp) in enumerate(ival):
            for l, (j, jexp) in enumerate(ival[k:]):
                if i == j:
                    f_jrs[k, i, i] = iexp * (iexp - 1)
                else:
                    f_jrs[k, i, j] = iexp
                    f_jrs[k, j, i] = iexp
                    f_jrs[l, j, i] = jexp
                    f_jrs[l, i, j] = jexp

        f_rs = np.prod(f_jrs, axis=0)
        return Interaction(indices, exponents, f_r, m_jr, f_rs, m_jrs)

    def _make_interaction_arrays(self, Ws, n):
        """
        A hidden convenience function that splits each
        excess term into _summary_

        :param Ws: List of interactions in form input by the user.
        :type Ws: list of lists
        :param n: Number of interaction values (e.g. 3 for ESV interactions)
        :type n: integer
        :return: The processed interaction values, indices and exponents.
        :rtype: tuple of 2D numpy array and list of Interaction objects
        """
        values = np.empty((len(Ws), n))
        interactions = []
        for i, W in enumerate(Ws):
            values[i] = W[:n]
            interactions.append(self._make_prefactors_and_exponents(W[n:]))
        return values, interactions

    @material_property
    def _c_xs(self):
        """
        Calculates the product of molar fractions raised to the exponents in
        each polynomial interaction.
        :return: Interaction postfactors
            (not including the values of the interactions themselves)
        :rtype: 1D numpy array
        """
        c = np.empty(self.n_W_ESV + self.n_W_mbr)
        for i in range(self.n_W_ESV):
            c[i] = np.prod(
                np.power(self.trans_fractions[self.c_ESV[i].inds], self.c_ESV[i].expts),
                axis=0,
            )
        for i in range(self.n_W_mbr):
            j = i + self.n_W_ESV
            c[j] = np.prod(
                np.power(self.trans_fractions[self.c_mbr[i].inds], self.c_mbr[i].expts),
                axis=0,
            )

        return c

    @material_property
    def _dc_xsdq(self):
        """
        Calculates the first derivative of the product of
        molar fractions raised to the exponents in each polynomial interaction
        with respect to the transformed endmember proportions.
        :return: Interaction postfactors
            (not including the values of the interactions themselves)
        :rtype: 2D numpy array
        """
        c = np.empty((self.n_W_ESV + self.n_W_mbr, self.n_transformed_endmembers))
        for i in range(self.n_W_ESV):
            c[i] = self.c_ESV[i].f_r * np.prod(
                np.power(
                    self.trans_fractions[self.c_ESV[i].inds], self.c_ESV[i].m_jr.T
                ),
                axis=1,
            )
        for i in range(self.n_W_mbr):
            j = i + self.n_W_ESV
            c[j] = self.c_mbr[i].f_r * np.prod(
                np.power(
                    self.trans_fractions[self.c_mbr[i].inds], self.c_mbr[i].m_jr.T
                ),
                axis=1,
            )
        return c

    @material_property
    def _d2c_xsdqdq(self):
        """
        Calculates the second derivative of the product of
        molar fractions raised to the exponents in each polynomial interaction
        with respect to the transformed endmember proportions.
        :return: Interaction postfactors
            (not including the values of the interactions themselves)
        :rtype: 3D numpy array
        """
        c = np.empty(
            (
                self.n_W_ESV + self.n_W_mbr,
                self.n_transformed_endmembers,
                self.n_transformed_endmembers,
            )
        )
        for i in range(self.n_W_ESV):
            c[i] = self.c_ESV[i].f_rs * np.prod(
                np.power(
                    self.trans_fractions[self.c_ESV[i].inds], self.c_ESV[i].m_jrs.T
                ),
                axis=2,
            )
        for i in range(self.n_W_mbr):
            j = i + self.n_W_ESV
            c[j] = self.c_mbr[i].f_rs * np.prod(
                np.power(
                    self.trans_fractions[self.c_mbr[i].inds], self.c_mbr[i].m_jrs.T
                ),
                axis=2,
            )
        return c

    @material_property
    def _dqdn(self):
        """
        The first derivative of the transformed endmember proportions
        with respect to the original endmember amounts.

        :return: dqdn
        :rtype: 2D numpy array
        """
        return np.einsum(
            "ik, km->im",
            self.dqdp,
            dpdn(self.molar_fractions, 1.0, self.ones, self.eye),
        )

    @material_property
    def _dc_xsdn(self):
        """
        Calculates the first derivative of the product of
        molar fractions raised to the exponents in each polynomial interaction
        with respect to the endmember amounts.
        :return: Interaction postfactors
            (not including the values of the interactions themselves)
        :rtype: 2D numpy array
        """
        return np.einsum("ir, rm->im", self._dc_xsdq, self._dqdn) + np.einsum(
            "i, m->im", self._c_xs, self.ones
        )

    @material_property
    def _d2c_xsdndn(self):
        """
        Calculates the second derivative of the product of
        molar fractions raised to the exponents in each polynomial interaction
        with respect to the endmember amounts.
        :return: Interaction postfactors
            (not including the values of the interactions themselves)
        :rtype: 3D numpy array
        """
        n = 1.0
        a1 = np.einsum("irs, sn, rm", self._d2c_xsdqdq, self._dqdn, self._dqdn)
        a2a = 1.0 / n * (np.einsum("ir, rm, n", self._dc_xsdq, self._dqdn, self.ones))
        a2 = a2a + np.einsum("imn->inm", a2a)
        a3 = np.einsum(
            "ir, rk, knm->inm",
            self._dc_xsdq,
            self.dqdp,
            d2pdndn(self.molar_fractions, n * n, self.ones, self.eyeones),
        )
        return a1 + a2 + a3

    @material_property
    def _ESV_gradient_list(self):
        """
        Calculates the first derivative of the
        internal energy, entropy and volume
        with respect to the endmember amounts.

        :return: dEdn, dSdn and dVdn.
        :rtype: 2D numpy array
        """
        if self.W_ESV is None:
            return np.zeros((3, self.n_endmembers))
        else:
            return np.einsum("ij, ik->jk", self.W_ESV, self._dc_xsdn[: self.n_W_ESV])

    @material_property
    def _ESV_hessian_list(self):
        """
        Calculates the second derivative of the
        internal energy, entropy and volume
        with respect to the endmember amounts.

        :return: d2Edndn, d2Sdndn and d2Vdndn.
        :rtype: 3D numpy array
        """
        if self.W_ESV is None:
            return np.zeros((3, self.n_endmembers, self.n_endmembers))
        else:
            return np.einsum(
                "ij, ikl->jkl", self.W_ESV, self._d2c_xsdndn[: self.n_W_ESV]
            )

    @material_property
    def mbr_scalar_list(self):
        """
        :return: The value with which to multiply
        each interaction endmember property to get the
        total excess property.
        :rtype: 1D numpy array
        """
        if self.W_mbr is None:
            return np.zeros((0))
        else:
            return np.einsum("ij, i->j", self.W_mbr, self._c_xs[self.n_W_ESV :])

    @material_property
    def _mbr_gradient_list(self):
        """
        :return: The values with which to multiply
        each interaction endmember property to get the
        first derivatives of the total excess property
        with respect to the endmember amounts.
        :rtype: 2D numpy array
        """
        if self.W_mbr is None:
            return np.zeros((0, self.n_endmembers))
        else:
            return np.einsum("ij, ik->jk", self.W_mbr, self._dc_xsdn[self.n_W_ESV :])

    @material_property
    def _mbr_hessian_list(self):
        """
        :return: The values with which to multiply
        each interaction endmember property to get the
        second derivatives of the total excess property
        with respect to the endmember amounts.
        :rtype: 3D numpy array
        """
        if self.W_mbr is None:
            return np.zeros((0, self.n_endmembers, self.n_endmembers))
        else:
            return np.einsum(
                "ij, ikl->jkl", self.W_mbr, self._d2c_xsdndn[self.n_W_ESV :]
            )

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions):
        dEdn, dSdn, dVdn = self._ESV_gradient_list
        mbr_gradients = self._mbr_gradient_list
        mbr_gibbs = np.array([mbr.gibbs for mbr in self.interaction_endmembers])
        gibbs = dEdn - temperature * dSdn + pressure * dVdn
        gibbs += np.einsum("i, ij->j", mbr_gibbs, mbr_gradients)
        return gibbs

    def excess_partial_gibbs_free_energies(
        self, pressure, temperature, molar_fractions
    ):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions
        )
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions
        )
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )

        dSdn = self._ESV_gradient_list[1]
        mbr_gradients = self._mbr_gradient_list
        mbr_entropies = np.array([mbr.S for mbr in self.interaction_endmembers])

        non_ideal_entropies = dSdn + np.einsum("i, ij->j", mbr_entropies, mbr_gradients)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        dVdn = self._ESV_gradient_list[2]
        mbr_gradients = self._mbr_gradient_list
        mbr_volumes = np.array([mbr.V for mbr in self.interaction_endmembers])

        return dVdn + np.einsum("i, ij->j", mbr_volumes, mbr_gradients)

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )

        d2Edndn, d2Sdndn, d2Vdndn = self._ESV_hessian_list
        mbr_hessian = self._mbr_hessian_list
        mbr_gibbs = np.array([mbr.gibbs for mbr in self.interaction_endmembers])

        d2Gdndn = d2Edndn - temperature * d2Sdndn + pressure * d2Vdndn
        d2Gdndn += np.einsum("i, ijk->jk", mbr_gibbs, mbr_hessian)

        return d2Gdndn - temperature * ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )

        d2Sdndn = self._ESV_hessian_list[1]
        mbr_hessian = self._mbr_hessian_list
        mbr_entropies = np.array([mbr.S for mbr in self.interaction_endmembers])

        d2Sdndn += np.einsum("i, ijk->jk", mbr_entropies, mbr_hessian)

        return d2Sdndn + ideal_entropy_hessian

    def volume_hessian(self, pressure, temperature, molar_fractions):
        d2Vdndn = self._ESV_hessian_list[2]
        mbr_hessian = self._mbr_hessian_list
        mbr_volumes = np.array([mbr.V for mbr in self.interaction_endmembers])

        d2Vdndn += np.einsum("i, ijk->jk", mbr_volumes, mbr_hessian)

        return d2Vdndn

    def Cp_excess(self):
        mbr_scalar = self.mbr_scalar_list
        mbr_Cp = np.array(
            [mbr.molar_heat_capacity_p for mbr in self.interaction_endmembers]
        )
        return np.einsum("i, i", mbr_scalar, mbr_Cp)

    def alphaV_excess(self):
        mbr_scalar = self.mbr_scalar_list
        mbr_d2gibbsdpdt = np.array(
            [mbr.alpha * mbr.V for mbr in self.interaction_endmembers]
        )
        return np.einsum("i, i", mbr_scalar, mbr_d2gibbsdpdt)

    def VoverKT_excess(self):
        mbr_scalar = self.mbr_scalar_list
        mbr_d2gibbsdpdp = np.array(
            [
                mbr.V / mbr.isothermal_bulk_modulus_reuss
                for mbr in self.interaction_endmembers
            ]
        )
        return np.einsum("i, i", mbr_scalar, mbr_d2gibbsdpdp)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.0e-10:
            return np.exp(
                self._non_ideal_excess_partial_gibbs(
                    pressure, temperature, molar_fractions
                )
                / (constants.gas_constant * temperature)
            )
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(
            self, pressure, temperature, molar_fractions
        ) * self.activity_coefficients(pressure, temperature, molar_fractions)
