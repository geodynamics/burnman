# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from ..utils.chemistry import process_solution_chemistry
from .. import constants


def _ideal_activities_fct(molar_fractions, endmember_noccupancies, n_endmembers,
                          n_occupancies, site_multiplicities,
                          endmember_configurational_entropies):
    site_noccupancies = np.dot(molar_fractions, endmember_noccupancies)
    site_multiplicities = np.einsum('i, ij', molar_fractions,
                                    site_multiplicities)
    site_occupancies = site_noccupancies * inverseish(site_multiplicities)

    a = np.power(site_occupancies, endmember_noccupancies).prod(-1)
    normalisation_constants = np.exp(endmember_configurational_entropies
                                     / constants.gas_constant)
    return normalisation_constants * a


def _non_ideal_hessian_fct(phi, molar_fractions, n_endmembers, alpha, W):
    q = np.eye(n_endmembers) - phi*np.ones((n_endmembers, n_endmembers))
    sum_pa = np.dot(molar_fractions, alpha)
    hess = np.einsum('m, i, ij, jk, mk->im', -alpha/sum_pa, -alpha, q, W, q)
    hess += hess.T
    return hess


def _non_ideal_interactions_fct(phi, molar_fractions, n_endmembers, alpha, W):
    # -sum(sum(qi.qj.Wij*)
    # equation (2) of Holland and Powell 2003
    q = np.eye(n_endmembers) - phi*np.ones((n_endmembers, n_endmembers))
    # The following are equivalent to
    # np.einsum('i, ij, jk, ik->i', -self.alphas, q, self.Wx, q)
    Wint = -alpha * (q.dot(W)*q).sum(-1)
    return Wint


def _non_ideal_hessian_subreg(p, n_endmembers, Wijk):
    Id = np.identity(n_endmembers)
    IIp = np.einsum('il, jm, k->ijklm', Id, Id, p)
    Ipp = np.einsum('il, j, k->ijkl', Id, p, p)
    ppp = np.einsum('i, j, k->ijk', p, p, p)

    A = (IIp
         + np.transpose(IIp, axes=[0, 2, 1, 3, 4])
         + np.transpose(IIp, axes=[1, 0, 2, 3, 4])
         + np.transpose(IIp, axes=[1, 2, 0, 3, 4])
         + np.transpose(IIp, axes=[2, 1, 0, 3, 4])
         + np.transpose(IIp, axes=[2, 0, 1, 3, 4]))
    B = 2.*(Ipp
            + np.transpose(Ipp, axes=[1, 0, 2, 3])
            + np.transpose(Ipp, axes=[2, 1, 0, 3]))

    Asum = (A
            - B[:, :, :, :, None]
            - B[:, :, :, None, :]
            + 6.*ppp[:, :, :, None, None])
    hess = np.einsum('ijklm, ijk->lm', Asum, Wijk)
    return hess


def _non_ideal_interactions_subreg(p, n_endmembers, Wijk):
    Aijkl = np.einsum('li, j, k->ijkl', np.identity(n_endmembers), p, p)
    ppp = np.einsum('i, j, k->ijk', p, p, p)

    Asum = (Aijkl
            + np.transpose(Aijkl, axes=[1, 0, 2, 3])
            + np.transpose(Aijkl, axes=[1, 2, 0, 3])
            - 2*ppp[:, :, :, None])

    Wint = np.einsum('ijk, ijkl->l', Wijk, Asum)
    return Wint


def logish(x, eps=1.e-5):
    """
    2nd order series expansion of log(x) about eps:
    log(eps) - sum_k=1^infty (f_eps)^k / k
    Prevents infinities at x=0
    """
    f_eps = 1. - x/eps
    mask = x > eps
    ln = np.where(x <= eps, np.log(eps) - f_eps - f_eps*f_eps/2., 0.)
    ln[mask] = np.log(x[mask])
    return ln


def inverseish(x, eps=1.e-5):
    """
    1st order series expansion of 1/x about eps: 2/eps - x/eps/eps
    Prevents infinities at x=0
    """
    mask = x > eps
    oneoverx = np.where(x <= eps, 2./eps - x/eps/eps, 0.)
    oneoverx[mask] = 1./x[mask]
    return oneoverx


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

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        G_excess : float
            The excess Gibbs free energy
        """
        return np.dot(np.array(molar_fractions), self.excess_partial_gibbs_free_energies(pressure, temperature, molar_fractions))

    def excess_volume(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess volume of the solution.
        The base class implementation assumes that the excess volume is zero.

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        V_excess : float
            The excess volume of the solution
        """
        return np.dot(molar_fractions,
                      self.excess_partial_volumes(pressure, temperature, molar_fractions))

    def excess_entropy(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy of the solution.
        The base class implementation assumes that the excess entropy is zero.

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        S_excess : float
            The excess entropy of the solution
        """
        return np.dot(molar_fractions,
                      self.excess_partial_entropies(pressure, temperature, molar_fractions))

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess enthalpy of the solution.
        The base class implementation assumes that the excess enthalpy is zero.

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        H_excess : float
            The excess enthalpy of the solution
        """
        return (self.excess_gibbs_free_energy(pressure, temperature, molar_fractions)
                + temperature*self.excess_entropy(pressure, temperature, molar_fractions))

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Gibbs free energy for each endmember of the solution.
        The base class implementation assumes that the excess gibbs
        free energy is zero.

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        partial_G_excess : numpy array
            The excess Gibbs free energy of each endmember
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy for each endmember of the solution.
        The base class implementation assumes that the excess entropy
        is zero (true for mechanical solutions).

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        partial_S_excess : numpy array
            The excess entropy of each endmember
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess volume for each endmember of the solution.
        The base class implementation assumes that the excess volume
        is zero.

        Parameters
        ----------
        pressure : float
            Pressure at which to evaluate the solution model. [Pa]

        temperature : float
            Temperature at which to evaluate the solution. [K]

        molar_fractions : list of floats
            List of molar fractions of the different endmembers in solution

        Returns
        -------
        partial_V_excess : numpy array
            The excess volume of each endmember
        """
        return np.zeros_like(np.array(molar_fractions))


class MechanicalSolution (SolutionModel):

    """
    An extremely simple class representing a mechanical solution model.
    A mechanical solution experiences no interaction between endmembers.
    Therefore, unlike ideal solutions there is no entropy of mixing;
    the total gibbs free energy of the solution is equal to the
    dot product of the molar gibbs free energies and molar fractions
    of the constituent materials.
    """

    def __init__(self, endmembers):
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

    def excess_gibbs_free_energy(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_volume(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_entropy(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def activities(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)


class IdealSolution (SolutionModel):

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
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

        # Process solution chemistry
        process_solution_chemistry(self)

        self._calculate_endmember_configurational_entropies()

    def _calculate_endmember_configurational_entropies(self):
        S_conf = -(constants.gas_constant
                   * (self.endmember_noccupancies
                      * logish(self.endmember_occupancies)).sum(-1))
        self.endmember_configurational_entropies = S_conf

    def excess_partial_gibbs_free_energies(self, pressure, temperature,
                                           molar_fractions):
        return self._ideal_excess_partial_gibbs(temperature, molar_fractions)

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        return self._ideal_excess_partial_entropies(temperature,
                                                    molar_fractions)

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return np.zeros((self.n_endmembers))

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return -temperature*hess_S

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return hess_S

    def volume_hessian(self, pressure, temperature, molar_fractions):
        return np.zeros((len(molar_fractions), len(molar_fractions)))

    def _configurational_entropy(self, molar_fractions):
        site_noccupancies = np.einsum('i, ij', molar_fractions,
                                      self.endmember_noccupancies)
        site_multiplicities = np.einsum('i, ij', molar_fractions,
                                        self.site_multiplicities)
        site_occupancies = (site_noccupancies
                            * inverseish(site_multiplicities))
        conf_entropy = -(constants.gas_constant
                         * (site_noccupancies
                            * logish(site_occupancies)).sum(-1))
        return conf_entropy

    def _ideal_excess_partial_gibbs(self, temperature, molar_fractions):
        return -temperature*self._ideal_excess_partial_entropies(temperature, molar_fractions)

    def _ideal_excess_partial_entropies(self, temperature, molar_fractions):
        return -constants.gas_constant * self._log_ideal_activities(molar_fractions)

    def _ideal_entropy_hessian(self, temperature, molar_fractions):
        hessian = -constants.gas_constant * self._log_ideal_activity_derivatives(molar_fractions)
        return hessian

    def _log_ideal_activities(self, molar_fractions):
        site_noccupancies = np.einsum('i, ij', molar_fractions,
                                      self.endmember_noccupancies)
        site_multiplicities = np.einsum('i, ij', molar_fractions,
                                        self.site_multiplicities)

        lna = np.einsum('ij, j->i', self.endmember_noccupancies,
                        logish(site_noccupancies) - logish(site_multiplicities))

        normalisation_constants = (self.endmember_configurational_entropies
                                   / constants.gas_constant)
        return lna + normalisation_constants

    def _log_ideal_activity_derivatives(self, molar_fractions):
        site_noccupancies = np.einsum('i, ij', molar_fractions,
                                      self.endmember_noccupancies)
        site_multiplicities = np.einsum('i, ij', molar_fractions,
                                        self.site_multiplicities)

        dlnadp = (np.einsum('pj, qj, j->pq', self.endmember_noccupancies,
                            self.endmember_noccupancies,
                            inverseish(site_noccupancies))
                  - np.einsum('pj, qj, j->pq', self.endmember_noccupancies,
                              self.site_multiplicities,
                              inverseish(site_multiplicities)))

        return dlnadp

    def _ideal_activities(self, molar_fractions):
        return _ideal_activities_fct(molar_fractions,
                                     self.endmember_noccupancies,
                                     self.n_endmembers,
                                     self.n_occupancies,
                                     self.site_multiplicities,
                                     self.endmember_configurational_entropies)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def activities(self, pressure, temperature, molar_fractions):
        return self._ideal_activities(molar_fractions)


class AsymmetricRegularSolution (IdealSolution):

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

    def __init__(self, endmembers, alphas, energy_interaction,
                 volume_interaction=None, entropy_interaction=None):

        self.n_endmembers = len(endmembers)

        # Create array of van Laar parameters
        self.alphas = np.array(alphas)

        # Create 2D arrays of interaction parameters
        self.We = np.triu(2. / (self.alphas[:, np.newaxis] + self.alphas), 1)
        self.We[np.triu_indices(self.n_endmembers, 1)] *= np.array([i for row in energy_interaction
                                                                    for i in row])

        if entropy_interaction is not None:
            self.Ws = np.triu(2. / (self.alphas[:, np.newaxis] + self.alphas), 1)
            self.Ws[np.triu_indices(self.n_endmembers, 1)] *= np.array([i for row in entropy_interaction
                                                                        for i in row])
        else:
            self.Ws = np.zeros((self.n_endmembers, self.n_endmembers))

        if volume_interaction is not None:
            self.Wv = np.triu(2. / (self.alphas[:, np.newaxis] + self.alphas), 1)
            self.Wv[np.triu_indices(self.n_endmembers, 1)] *= np.array([i for row in volume_interaction
                                                                        for i in row])
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
        return _non_ideal_interactions_fct(phi, np.array(molar_fractions), self.n_endmembers, self.alphas, W)

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions):
        Eint = self._non_ideal_interactions(self.We, molar_fractions)
        Sint = self._non_ideal_interactions(self.Ws, molar_fractions)
        Vint = self._non_ideal_interactions(self.Wv, molar_fractions)
        return Eint - temperature * Sint + pressure * Vint

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions)
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions)
        non_ideal_entropies = self._non_ideal_interactions(self.Ws, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        return self._non_ideal_interactions(self.Wv, molar_fractions)

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(self, temperature, molar_fractions)
        phi = self._phi(molar_fractions)
        nonideal_gibbs_hessian = _non_ideal_hessian_fct(phi, molar_fractions,
                                                        self.n_endmembers, self.alphas,
                                                        self.We - temperature*self.Ws + pressure*self.Wv)

        return nonideal_gibbs_hessian - temperature*ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(self, temperature, molar_fractions)
        phi = self._phi(molar_fractions)
        nonideal_entropy_hessian = _non_ideal_hessian_fct(phi, molar_fractions,
                                                          self.n_endmembers, self.alphas,
                                                          self.Ws)
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def volume_hessian(self, pressure, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        return _non_ideal_hessian_fct(phi, molar_fractions,
                                      self.n_endmembers, self.alphas,
                                      self.Wv)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.e-10:
            return np.exp(self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions) / (constants.gas_constant * temperature))
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(self, pressure, temperature, molar_fractions) * self.activity_coefficients(pressure, temperature, molar_fractions)


class SymmetricRegularSolution (AsymmetricRegularSolution):

    """
    Solution model implementing the symmetric regular solution model.
    This is a special case of the
    :class:`burnman.solutionmodel.AsymmetricRegularSolution` class.
    """

    def __init__(self, endmembers, energy_interaction, volume_interaction=None, entropy_interaction=None):
        alphas = np.ones(len(endmembers))
        AsymmetricRegularSolution.__init__(
            self, endmembers, alphas, energy_interaction, volume_interaction, entropy_interaction)


class SubregularSolution (IdealSolution):

    """
    Solution model implementing the subregular solution model formulation
    as described in :cite:`HW1989`. The excess conconfigurational
    Gibbs energy is given by the expression:

    .. math::
        \\mathcal{G}_{\\textrm{excess}} = \\sum_i \\sum_{j > i} (p_i p_j^2
        W_{ij} + p_j p_i^2 W_{ji} + \\sum_{k > j > i} p_i p_j p_k W_{ijk})

    Interaction parameters are inserted into a 3D interaction matrix during
    initialization to make use of numpy vector algebra.

    Parameters
    ----------
    endmembers : list of lists
        A list of all the independent endmembers in the solution.
        The first item of each list gives the Mineral object corresponding
        to the endmember. The second item gives the site-species formula.

    energy_interaction : list of list of lists
        The binary endmember interaction energies.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).

    volume_interaction : list of list of lists
        The binary endmember interaction volumes.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).

    entropy_interaction : list of list of lists
        The binary endmember interaction entropies.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).

    energy_ternary_terms : list of lists
        The ternary interaction energies. Each list should contain
        four entries: the indices i, j, k and the value of the interaction.

    volume_ternary_terms : list of lists
        The ternary interaction volumes. Each list should contain
        four entries: the indices i, j, k and the value of the interaction.

    entropy_ternary_terms : list of lists
        The ternary interaction entropies. Each list should contain
        four entries: the indices i, j, k and the value of the interaction.
    """

    def __init__(self, endmembers, energy_interaction,
                 volume_interaction=None, entropy_interaction=None,
                 energy_ternary_terms=None,
                 volume_ternary_terms=None, entropy_ternary_terms=None):
        """
        Initialization function for the SubregularSolution class.
        """

        self.n_endmembers = len(endmembers)

        # Create 3D arrays of interaction parameters
        self.Wijke = np.zeros(shape=(self.n_endmembers,
                                     self.n_endmembers,
                                     self.n_endmembers))
        self.Wijks = np.zeros_like(self.Wijke)
        self.Wijkv = np.zeros_like(self.Wijke)

        # setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i + 1, self.n_endmembers):
                w0 = energy_interaction[i][j - i - 1][0]/2.
                w1 = energy_interaction[i][j - i - 1][1]/2.
                self.Wijke[:, i, j] += w0
                self.Wijke[:, j, i] += w1

                self.Wijke[i, j, j] += w0
                self.Wijke[j, i, i] += w1

                self.Wijke[i, j, i] -= w0
                self.Wijke[j, i, j] -= w1

        if energy_ternary_terms is not None:
            for (i, j, k, v) in energy_ternary_terms:
                self.Wijke[i, j, k] += v

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    w0 = entropy_interaction[i][j - i - 1][0]/2.
                    w1 = entropy_interaction[i][j - i - 1][1]/2.
                    self.Wijks[:, i, j] += w0
                    self.Wijks[:, j, i] += w1

                    self.Wijks[i, j, j] += w0
                    self.Wijks[j, i, i] += w1

                    self.Wijks[i, j, i] -= w0
                    self.Wijks[j, i, j] -= w1

        if entropy_ternary_terms is not None:
            for (i, j, k, v) in entropy_ternary_terms:
                self.Wijks[i, j, k] += v

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    w0 = volume_interaction[i][j - i - 1][0]/2.
                    w1 = volume_interaction[i][j - i - 1][1]/2.
                    self.Wijkv[:, i, j] += w0
                    self.Wijkv[:, j, i] += w1

                    self.Wijkv[i, j, j] += w0
                    self.Wijkv[j, i, i] += w1

                    self.Wijkv[i, j, i] -= w0
                    self.Wijkv[j, i, j] -= w1

        if volume_ternary_terms is not None:
            for (i, j, k, v) in volume_ternary_terms:
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

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions)
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_partial_entropies(self, pressure, temperature, molar_fractions):
        ideal_entropies = IdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions)
        non_ideal_entropies = self._non_ideal_function(self.Wijks, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_volumes(self, pressure, temperature, molar_fractions):
        non_ideal_volumes = self._non_ideal_function(self.Wijkv, molar_fractions)
        return non_ideal_volumes

    def gibbs_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(self, temperature, molar_fractions)
        nonideal_gibbs_hessian = _non_ideal_hessian_subreg(molar_fractions, n,
                                                           self.Wijke - temperature*self.Wijks
                                                           + pressure*self.Wijkv)

        return nonideal_gibbs_hessian - temperature*ideal_entropy_hessian

    def entropy_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = IdealSolution._ideal_entropy_hessian(self, temperature, molar_fractions)
        nonideal_entropy_hessian = _non_ideal_hessian_subreg(molar_fractions, n,
                                                             self.Wijks)
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def volume_hessian(self, pressure, temperature, molar_fractions):
        n = len(molar_fractions)
        return _non_ideal_hessian_subreg(molar_fractions, n, self.Wijkv)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.e-10:
            return np.exp(self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions) / (constants.gas_constant * temperature))
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(self, pressure, temperature, molar_fractions) * self.activity_coefficients(pressure, temperature, molar_fractions)
