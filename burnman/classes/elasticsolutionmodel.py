# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import importlib
import warnings
import numpy as np
from ..utils.chemistry import process_solution_chemistry
from .solutionmodel import _ideal_activities_fct
from .solutionmodel import _non_ideal_hessian_fct, _non_ideal_interactions_fct
from .solutionmodel import _non_ideal_hessian_subreg
from .solutionmodel import _non_ideal_interactions_subreg
from .solutionmodel import logish, inverseish
from .. import constants

try:
    ag = importlib.import_module("autograd")
except ImportError as err:
    print(
        f"Warning: {err}. "
        "For full functionality of BurnMan, please install autograd."
    )


class ElasticSolutionModel(object):
    """
    This is the base class for an Elastic solution model, intended for use
    in defining solutions and performing thermodynamic calculations
    on them.  All minerals of type :class:`burnman.Solution` use
    a solution model for defining how the endmembers in the solution
    interact.

    A user wanting a new solution model should define the functions included
    in the base class. All of the functions in the base class return zero,
    so if the user-defined solution model does not implement them,
    they essentially have no effect, and the Helmholtz energy and
    pressure of a solution will be equal to the weighted arithmetic
    averages of the different endmember values.
    """

    def __init__(self):
        """
        Does nothing.
        """
        pass

    def excess_helmholtz_energy(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Helmholtz free energy of the solution.
        The base class implementation assumes that the excess Helmholtz
        energy is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess Helmholtz energy.
        :rtype: float
        """
        return np.dot(
            np.array(molar_fractions),
            self.excess_partial_helmholtz_energies(
                volume, temperature, molar_fractions
            ),
        )

    def excess_pressure(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess pressure of the solution.
        The base class implementation assumes that the excess pressure is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess pressure of the solution.
        :rtype: float
        """
        return np.dot(
            molar_fractions,
            self.excess_partial_pressures(volume, temperature, molar_fractions),
        )

    def excess_entropy(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy of the solution.
        The base class implementation assumes that the excess entropy is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess entropy of the solution.
        :rtype: float
        """
        return np.dot(
            molar_fractions,
            self.excess_partial_entropies(volume, temperature, molar_fractions),
        )

    def excess_enthalpy(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess enthalpy of the solution.
        The base class implementation assumes that the excess enthalpy is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess enthalpy of the solution.
        :rtype: float
        """
        return (
            self.excess_helmholtz_energy(volume, temperature, molar_fractions)
            + temperature * self.excess_entropy(volume, temperature, molar_fractions)
            - volume * self.excess_pressure(volume, temperature, molar_fractions)
        )

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Helmholtz energy for each endmember of the solution.
        The base class implementation assumes that the excess Helmholtz energy
        is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess Helmholtz energy of each endmember
        :rtype: numpy.array
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess entropy for each endmember of the solution.
        The base class implementation assumes that the excess entropy
        is zero (true for mechanical solutions).

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess entropy of each endmember.
        :rtype: numpy.array
        """
        return np.zeros_like(molar_fractions)

    def excess_partial_pressures(self, volume, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess pressure for each endmember of the solution.
        The base class implementation assumes that the excess pressure
        is zero.

        :param volume: Volume at which to evaluate the solution model. [m^3/mol]
        :type volume: float

        :param temperature: Temperature at which to evaluate the solution. [K]
        :type temperature: float

        :param molar_fractions: List of molar fractions of the different
            endmembers in solution.
        :type molar_fractions: list of floats

        :returns: The excess pressure of each endmember.
        :rtype: numpy.array
        """
        return np.zeros_like(np.array(molar_fractions))


class ElasticMechanicalSolution(ElasticSolutionModel):
    """
    An extremely simple class representing a mechanical solution model.
    A mechanical solution experiences no interaction between endmembers.
    Therefore, unlike ideal solutions there is no entropy of mixing;
    the total Helmholtz energy of the solution is equal to the
    dot product of the molar Helmholtz energies and molar fractions
    of the constituent materials.
    """

    def __init__(self, endmembers):
        self.endmembers = endmembers
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

    def excess_helmholtz_energy(self, volume, temperature, molar_fractions):
        return 0.0

    def excess_pressure(self, volume, temperature, molar_fractions):
        return 0.0

    def excess_entropy(self, volume, temperature, molar_fractions):
        return 0.0

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def excess_partial_pressures(self, volume, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)


class ElasticIdealSolution(ElasticSolutionModel):
    """
    A class representing an ideal solution model.
    Calculates the excess Helmholtz energy and entropy due to configurational
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
            * (self.endmember_noccupancies * logish(self.endmember_occupancies)).sum(-1)
        )
        self.endmember_configurational_entropies = S_conf

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        return self._ideal_excess_partial_helmholtz(temperature, molar_fractions)

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        return self._ideal_excess_partial_entropies(temperature, molar_fractions)

    def excess_partial_pressures(self, volume, temperature, molar_fractions):
        return np.zeros((self.n_endmembers))

    def helmholtz_hessian(self, volume, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return -temperature * hess_S

    def entropy_hessian(self, volume, temperature, molar_fractions):
        hess_S = self._ideal_entropy_hessian(temperature, molar_fractions)
        return hess_S

    def pressure_hessian(self, volume, temperature, molar_fractions):
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

    def _ideal_excess_partial_helmholtz(self, temperature, molar_fractions):
        return -(
            temperature
            * self._ideal_excess_partial_entropies(temperature, molar_fractions)
        )

    def _ideal_excess_partial_entropies(self, temperature, molar_fractions):
        return -(constants.gas_constant * self._log_ideal_activities(molar_fractions))

    def _ideal_entropy_hessian(self, temperature, molar_fractions):
        hessian = -constants.gas_constant * self._log_ideal_activity_derivatives(
            molar_fractions
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


class ElasticAsymmetricRegularSolution(ElasticIdealSolution):
    """
    Solution model implementing the asymmetric regular solution model
    formulation as described in :cite:`HP2003`.

    The excess nonconfigurational Helmholtz energy is given by the
    expression:

    .. math::
        \\mathcal{F}_{\\textrm{excess}} = \\alpha^T p (\\phi^T W \\phi)

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
        pressure_interaction=None,
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

        if pressure_interaction is not None:
            self.Wp = np.triu(2.0 / (self.alphas[:, np.newaxis] + self.alphas), 1)
            self.Wp[np.triu_indices(self.n_endmembers, 1)] *= np.array(
                [i for row in pressure_interaction for i in row]
            )
        else:
            self.Wp = np.zeros((self.n_endmembers, self.n_endmembers))

        # initialize ideal solution model
        ElasticIdealSolution.__init__(self, endmembers)

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

    def _non_ideal_excess_partial_helmholtz(self, volume, temperature, molar_fractions):
        Eint = self._non_ideal_interactions(self.We, molar_fractions)
        Sint = self._non_ideal_interactions(self.Ws, molar_fractions)
        Pint = self._non_ideal_interactions(self.Wp, molar_fractions)
        return Eint - temperature * Sint - volume * Pint

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        ideal_helmholtz = ElasticIdealSolution._ideal_excess_partial_helmholtz(
            self, temperature, molar_fractions
        )
        non_ideal_helmholtz = self._non_ideal_excess_partial_helmholtz(
            volume, temperature, molar_fractions
        )
        return ideal_helmholtz + non_ideal_helmholtz

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        ideal_entropies = ElasticIdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_interactions(self.Ws, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_pressures(self, volume, temperature, molar_fractions):
        return self._non_ideal_interactions(self.Wp, molar_fractions)

    def helmholtz_hessian(self, volume, temperature, molar_fractions):
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        phi = self._phi(molar_fractions)
        nonideal_helmholtz_hessian = _non_ideal_hessian_fct(
            phi,
            molar_fractions,
            self.n_endmembers,
            self.alphas,
            self.We - temperature * self.Ws - volume * self.Wp,
        )

        return nonideal_helmholtz_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, volume, temperature, molar_fractions):
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        phi = self._phi(molar_fractions)
        nonideal_entropy_hessian = _non_ideal_hessian_fct(
            phi, molar_fractions, self.n_endmembers, self.alphas, self.Ws
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def pressure_hessian(self, volume, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        return _non_ideal_hessian_fct(
            phi, molar_fractions, self.n_endmembers, self.alphas, self.Wp
        )


class ElasticSymmetricRegularSolution(ElasticAsymmetricRegularSolution):
    """
    Solution model implementing the symmetric regular solution model.
    This is a special case of the
    :class:`burnman.solutionmodel.AsymmetricRegularSolution` class.
    """

    def __init__(
        self,
        endmembers,
        energy_interaction,
        pressure_interaction=None,
        entropy_interaction=None,
    ):
        alphas = np.ones(len(endmembers))
        ElasticAsymmetricRegularSolution.__init__(
            self,
            endmembers,
            alphas,
            energy_interaction,
            pressure_interaction,
            entropy_interaction,
        )


class ElasticSubregularSolution(ElasticIdealSolution):
    """
    Solution model implementing the subregular solution model formulation
    as described in :cite:`HW1989`. The excess conconfigurational
    Helmholtz energy is given by the expression:

    .. math::
        \\mathcal{F}_{\\textrm{excess}} = \\sum_i \\sum_{j > i} (p_i p_j^2
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
    :param pressure_interaction: The binary endmember interaction pressures.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).
    :type pressure_interaction: list of list of lists
    :param entropy_interaction: The binary endmember interaction entropies.
        Each interaction[i, j-i-1, 0] corresponds to W(i,j), while
        interaction[i, j-i-1, 1] corresponds to W(j,i).
    :type entropy_interaction: list of list of lists
    :param energy_ternary_terms: The ternary interaction energies.
        Each list should contain four entries: the indices i, j, k
        and the value of the interaction.
    :type energy_ternary_terms: list of lists
    :param pressure_ternary_terms: The ternary interaction pressures.
        Each list should contain four entries: the indices i, j, k
        and the value of the interaction.
    :type pressure_ternary_terms: list of lists
    :param entropy_ternary_terms: The ternary interaction entropies.
        Each list should contain four entries:
        the indices i, j, k and the value of the interaction.
    :type entropy_ternary_terms: list of lists
    """

    def __init__(
        self,
        endmembers,
        energy_interaction,
        pressure_interaction=None,
        entropy_interaction=None,
        energy_ternary_terms=None,
        pressure_ternary_terms=None,
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
        self.Wijkp = np.zeros_like(self.Wijke)

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

        if pressure_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    w0 = pressure_interaction[i][j - i - 1][0] / 2.0
                    w1 = pressure_interaction[i][j - i - 1][1] / 2.0
                    self.Wijkp[:, i, j] += w0
                    self.Wijkp[:, j, i] += w1

                    self.Wijkp[i, j, j] += w0
                    self.Wijkp[j, i, i] += w1

                    self.Wijkp[i, j, i] -= w0
                    self.Wijkp[j, i, j] -= w1

        if pressure_ternary_terms is not None:
            for i, j, k, v in pressure_ternary_terms:
                self.Wijkv[i, j, k] += v

        # initialize ideal solution model
        ElasticIdealSolution.__init__(self, endmembers)

    def _non_ideal_function(self, Wijk, molar_fractions):
        n = len(molar_fractions)
        return _non_ideal_interactions_subreg(molar_fractions, n, Wijk)

    def _non_ideal_interactions(self, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        Eint = self._non_ideal_function(self.Wijke, molar_fractions)
        Sint = self._non_ideal_function(self.Wijks, molar_fractions)
        Pint = self._non_ideal_function(self.Wijkp, molar_fractions)
        return Eint, Sint, Pint

    def _non_ideal_excess_partial_helmholtz(self, volume, temperature, molar_fractions):
        Eint, Sint, Pint = self._non_ideal_interactions(molar_fractions)
        return Eint - temperature * Sint - volume * Pint

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        ideal_helmholtz = ElasticIdealSolution._ideal_excess_partial_helmholtz(
            self, temperature, molar_fractions
        )
        non_ideal_helmholtz = self._non_ideal_excess_partial_helmholtz(
            volume, temperature, molar_fractions
        )
        return ideal_helmholtz + non_ideal_helmholtz

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        ideal_entropies = ElasticIdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_function(self.Wijks, molar_fractions)
        return ideal_entropies + non_ideal_entropies

    def excess_partial_pressures(self, volume, temperature, molar_fractions):
        non_ideal_pressures = self._non_ideal_function(self.Wijkp, molar_fractions)
        return non_ideal_pressures

    def helmholtz_hessian(self, volume, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_helmholtz_hessian = _non_ideal_hessian_subreg(
            molar_fractions,
            n,
            self.Wijke - temperature * self.Wijks - volume * self.Wijkp,
        )

        return nonideal_helmholtz_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, volume, temperature, molar_fractions):
        n = len(molar_fractions)
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_entropy_hessian = _non_ideal_hessian_subreg(
            molar_fractions, n, self.Wijks
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian

    def pressure_hessian(self, volume, temperature, molar_fractions):
        n = len(molar_fractions)
        return _non_ideal_hessian_subreg(molar_fractions, n, self.Wijkp)


class ElasticFunctionSolution(ElasticIdealSolution):
    """
    Solution model implementing a generalized elastic solution model.
    The extensive excess nonconfigurational Helmholtz energy is
    provided as a function by the user.

    Derivatives are calculated using the autograd module,
    and so the user-defined excess Helmholtz energy function
    should be defined using autograd-friendly expressions.

    :param endmembers: A list of all the independent endmembers in the solution.
        The first item of each list gives the Mineral object corresponding
        to the endmember. The second item gives the site-species formula.
    :type endmembers: list of lists

    :param excess_helmholtz_function: The nonconfigurational
        Helmholtz energy function with arguments volume,
        temperature and molar_amounts, in that order.
        Note that the function must be extensive; if the molar amounts
        are doubled, the Helmholtz energy must also double.
    :type excess_helmholtz_function: function
    """

    def __init__(self, endmembers, excess_helmholtz_function):
        """
        Initialization function for the GeneralSolution class.
        """

        # initialize ideal solution model
        ElasticIdealSolution.__init__(self, endmembers)

        self.n_endmembers = len(endmembers)
        self._excess_helmholtz_function = excess_helmholtz_function

        self._non_ideal_excess_partial_helmholtz = ag.jacobian(
            excess_helmholtz_function, argnum=2
        )

        def partial_entropies(volume, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return -ag.jacobian(self._non_ideal_excess_partial_helmholtz, argnum=1)(
                    volume, temperature, molar_amounts
                )

        self._non_ideal_excess_partial_entropies = partial_entropies

        def partial_pressures(volume, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return -ag.jacobian(self._non_ideal_excess_partial_helmholtz, argnum=0)(
                    volume, temperature, molar_amounts
                )

        self.excess_partial_pressures = partial_pressures

        self._non_ideal_helmholtz_hessian = ag.jacobian(
            self._non_ideal_excess_partial_helmholtz, argnum=2
        )

        def entropy_hess(volume, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return ag.jacobian(partial_entropies, argnum=2)(
                    volume, temperature, molar_amounts
                )

        self._non_ideal_entropy_hessian = entropy_hess

        def pressure_hess(volume, temperature, molar_amounts):
            with warnings.catch_warnings(record=True):
                warnings.simplefilter("always")
                return ag.jacobian(partial_pressures, argnum=2)(
                    volume, temperature, molar_amounts
                )

        self.pressure_hessian = pressure_hess

    def excess_partial_helmholtz_energies(self, volume, temperature, molar_fractions):
        ideal_helmholtz = ElasticIdealSolution._ideal_excess_partial_helmholtz(
            self, temperature, molar_fractions
        )
        non_ideal_helmholtz = self._non_ideal_excess_partial_helmholtz(
            volume, temperature, molar_fractions
        )
        return ideal_helmholtz + non_ideal_helmholtz

    def excess_partial_entropies(self, volume, temperature, molar_fractions):
        ideal_entropies = ElasticIdealSolution._ideal_excess_partial_entropies(
            self, temperature, molar_fractions
        )
        non_ideal_entropies = self._non_ideal_excess_partial_entropies(
            volume, temperature, molar_fractions
        )
        return ideal_entropies + non_ideal_entropies

    def helmholtz_hessian(self, volume, temperature, molar_fractions):
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_helmholtz_hessian = self._non_ideal_helmholtz_hessian(
            volume, temperature, molar_fractions
        )

        return nonideal_helmholtz_hessian - temperature * ideal_entropy_hessian

    def entropy_hessian(self, volume, temperature, molar_fractions):
        ideal_entropy_hessian = ElasticIdealSolution._ideal_entropy_hessian(
            self, temperature, molar_fractions
        )
        nonideal_entropy_hessian = self._non_ideal_entropy_hessian(
            volume, temperature, molar_fractions
        )
        return ideal_entropy_hessian + nonideal_entropy_hessian
