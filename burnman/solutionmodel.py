# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from . import processchemistry
from . import constants

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    import os
    if 'NUMBA_DISABLE_JIT' in os.environ and int(os.environ['NUMBA_DISABLE_JIT']) == 1:
        raise ImportError("NOOOO!")
    from numba import jit
except ImportError:
    def jit(fn):
        return fn

@jit
def _ideal_activities_fct(molar_fractions, endmember_occupancies, n_endmembers, n_occupancies, site_multiplicities, endmember_configurational_entropies):
    site_occupancies = np.dot(molar_fractions, endmember_occupancies)
    activities = np.empty(shape=(n_endmembers))
    
    for e in range(n_endmembers):
        a = 1.0
        for occ in range(n_occupancies):
            if endmember_occupancies[e][occ] > 1e-10:
                a *= np.power(
                        site_occupancies[occ], endmember_occupancies[e][occ] * site_multiplicities[occ])
        normalisation_constant = np.exp(
            endmember_configurational_entropies[e] / constants.gas_constant)
        activities[e] = normalisation_constant * a
    return activities

@jit
def _non_ideal_interactions_fct(phi, molar_fractions, n_endmembers, alpha, We, Ws, Wv):
    # -sum(sum(qi.qj.Wij*)
    # equation (2) of Holland and Powell 2003
    Eint = np.zeros(len(molar_fractions))
    Sint = np.zeros(len(molar_fractions))
    Vint = np.zeros(len(molar_fractions))
    
    for l in range(n_endmembers):
        q = -phi
        q[l] += 1.0

        Eint[l] = 0. - alpha[l] * np.dot(q, np.dot(We, q))
        Sint[l] = 0. - alpha[l] * np.dot(q, np.dot(Ws, q))
        Vint[l] = 0. - alpha[l] * np.dot(q, np.dot(Wv, q))
        
    return Eint, Sint, Vint


class SolutionModel(object):

    """
    This is the base class for a solution model,  intended for use
    in defining solid solutions and performing thermodynamic calculations
    on them.  All minerals of type :class:`burnman.SolidSolution` use
    a solution model for defining how the endmembers in the solid solution
    interact.

    A user wanting a new solution model should define the functions below.
    In the base class all of these return zero, so if the solution model
    does not implement them, they essentially have no effect, and
    then the Gibbs free energy and molar volume of a solid solution are
    just the weighted arithmetic averages of the different endmember values.
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
        return np.empty_like(np.array(molar_fractions))

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
        return 0.0

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
        return 0.0

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
        return 0.0

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
    
    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        return np.zeros_like(molar_fractions)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def activities(self, pressure, temperature, molar_fractions):
        return np.ones_like(molar_fractions)

    def excess_volume(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_entropy(self, pressure, temperature, molar_fractions):
        return 0.

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        return 0.

    
class IdealSolution (SolutionModel):

    """
    A very simple class representing an ideal solution model.
    Calculate the excess gibbs free energy due to configurational
    entropy, all the other excess terms return zero.
    """

    def __init__(self, endmembers):
        self.n_endmembers = len(endmembers)
        self.formulas = [e[1] for e in endmembers]

        # Process solid solution chemistry
        self.solution_formulae, self.n_sites, self.sites, self.n_occupancies, self.endmember_occupancies, self.site_multiplicities = \
            processchemistry.process_solution_chemistry(self.formulas)

        self._calculate_endmember_configurational_entropies()

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        return self._ideal_excess_partial_gibbs(temperature, molar_fractions)

    def _calculate_endmember_configurational_entropies(self):
        self.endmember_configurational_entropies = np.zeros(
            shape=(self.n_endmembers))
        for idx, endmember_occupancy in enumerate(self.endmember_occupancies):
            for occ in range(self.n_occupancies):
                if endmember_occupancy[occ] > 1e-10:
                    self.endmember_configurational_entropies[idx] = \
                        self.endmember_configurational_entropies[idx] - \
                        constants.gas_constant * self.site_multiplicities[
                            occ] * endmember_occupancy[occ] * np.log(endmember_occupancy[occ])

    def _endmember_configurational_entropy_contribution(self, molar_fractions):
        return np.dot(molar_fractions, self.endmember_configurational_entropies)

    def _configurational_entropy(self, molar_fractions):
        site_occupancies = np.dot(molar_fractions, self.endmember_occupancies)
        conf_entropy = 0
        for idx, occupancy in enumerate(site_occupancies):
            if occupancy > 1e-10:
                conf_entropy = conf_entropy - constants.gas_constant * \
                    occupancy * \
                    self.site_multiplicities[idx] * np.log(occupancy)

        return conf_entropy

    def _ideal_excess_partial_gibbs(self, temperature, molar_fractions):
        return constants.gas_constant * temperature * self._log_ideal_activities(molar_fractions)

    def _log_ideal_activities(self, molar_fractions):
        site_occupancies = np.dot(molar_fractions, self.endmember_occupancies)
        lna = np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            lna[e] = 0.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] > 1e-10 and site_occupancies[occ] > 1e-10:
                    lna[e] = lna[e] + self.endmember_occupancies[e][occ] * \
                        self.site_multiplicities[
                            occ] * np.log(site_occupancies[occ])

            normalisation_constant = self.endmember_configurational_entropies[
                e] / constants.gas_constant
            lna[e] = lna[e] + self.endmember_configurational_entropies[
                e] / constants.gas_constant
        return lna


    def _ideal_activities(self, molar_fractions):
        return _ideal_activities_fct(molar_fractions, 
                                     self.endmember_occupancies, 
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
    Solution model implementing the asymmetric regular solution model formulation (Holland and Powell, 2003)
    """

    def __init__(self, endmembers, alphas, energy_interaction, volume_interaction=None, entropy_interaction=None):

        self.n_endmembers = len(endmembers)

        # Create array of van Laar parameters
        self.alpha = np.array(alphas)

        # Create 2D arrays of interaction parameters
        self.We = np.zeros(shape=(self.n_endmembers, self.n_endmembers))
        self.Ws = np.zeros(shape=(self.n_endmembers, self.n_endmembers))
        self.Wv = np.zeros(shape=(self.n_endmembers, self.n_endmembers))

        # setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i + 1, self.n_endmembers):
                self.We[i][j] = 2. * energy_interaction[
                    i][j - i - 1] / (self.alpha[i] + self.alpha[j])

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    self.Ws[i][j] = 2. * entropy_interaction[
                        i][j - i - 1] / (self.alpha[i] + self.alpha[j])

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    self.Wv[i][j] = 2. * volume_interaction[
                        i][j - i - 1] / (self.alpha[i] + self.alpha[j])

        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def _phi(self, molar_fractions):
        phi = np.array([self.alpha[i] * molar_fractions[i]
                       for i in range(self.n_endmembers)])
        phi = np.divide(phi, np.sum(phi))
        return phi


    def _non_ideal_interactions(self, molar_fractions):
        # -sum(sum(qi.qj.Wij*)
        # equation (2) of Holland and Powell 2003
        phi = self._phi(molar_fractions)
        return _non_ideal_interactions_fct(phi, molar_fractions, self.n_endmembers, self.alpha, self.We, self.Ws, self.Wv)

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions):
        Eint, Sint, Vint = self._non_ideal_interactions(molar_fractions)
        return Eint - temperature * Sint + pressure * Vint

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs(
            self, temperature, molar_fractions)
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(
            pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_volume(self, pressure, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        V_excess = np.dot(self.alpha.T, molar_fractions) * np.dot(
            phi.T, np.dot(self.Wv, phi))
        return V_excess

    def excess_entropy(self, pressure, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        S_conf = -constants.gas_constant * \
            np.dot(IdealSolution._log_ideal_activities(
                self, molar_fractions), molar_fractions)
        S_excess = np.dot(self.alpha.T, molar_fractions) * np.dot(
            phi.T, np.dot(self.Ws, phi))
        return S_conf + S_excess

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        phi = self._phi(molar_fractions)
        E_excess = np.dot(self.alpha.T, molar_fractions) * np.dot(
            phi.T, np.dot(self.We, phi))
        return E_excess + pressure * self.excess_volume(pressure, temperature, molar_fractions)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.e-10:
            return np.exp(self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions) / (constants.gas_constant * temperature))
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(self, pressure, temperature, molar_fractions) * self.activity_coefficients(pressure, temperature, molar_fractions)


class SymmetricRegularSolution (AsymmetricRegularSolution):

    """
    Solution model implementing the symmetric regular solution model
    """

    def __init__(self, endmembers, energy_interaction, volume_interaction=None, entropy_interaction=None):
        alphas = np.ones(len(endmembers))
        AsymmetricRegularSolution.__init__(
            self, endmembers, alphas, energy_interaction, volume_interaction, entropy_interaction)


class SubregularSolution (IdealSolution):

    """
    Solution model implementing the subregular solution model formulation (Helffrich and Wood, 1989)
    """

    def __init__(self, endmembers, energy_interaction, volume_interaction=None, entropy_interaction=None):

        self.n_endmembers = len(endmembers)

        # Create 2D arrays of interaction parameters
        self.We = np.zeros(shape=(self.n_endmembers, self.n_endmembers))
        self.Ws = np.zeros(shape=(self.n_endmembers, self.n_endmembers))
        self.Wv = np.zeros(shape=(self.n_endmembers, self.n_endmembers))

        # setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i + 1, self.n_endmembers):
                self.We[i][j] = energy_interaction[i][j - i - 1][0]
                self.We[j][i] = energy_interaction[i][j - i - 1][1]

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    self.Ws[i][j] = entropy_interaction[i][j - i - 1][0]
                    self.Ws[j][i] = entropy_interaction[i][j - i - 1][1]

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i + 1, self.n_endmembers):
                    self.Wv[i][j] = volume_interaction[i][j - i - 1][0]
                    self.Wv[j][i] = volume_interaction[i][j - i - 1][1]

        # initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def _non_ideal_function(self, W, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        n = len(molar_fractions)
        RTlny = np.zeros(n)
        for l in range(n):
            val = 0.
            for i in range(n):
                if i != l:
                    val += 0.5 * molar_fractions[i] * (W[l][i] * (1 - molar_fractions[l] + molar_fractions[i] + 2. * molar_fractions[l] * (molar_fractions[l] - molar_fractions[i] - 1)) + W[
                                                       i][l] * (1. - molar_fractions[l] - molar_fractions[i] - 2. * molar_fractions[l] * (molar_fractions[l] - molar_fractions[i] - 1)))
                    for j in range(i + 1, n):
                        if j != l:
                            val += molar_fractions[i] * molar_fractions[j] * (
                                W[i][j] * (molar_fractions[i] - molar_fractions[j] - 0.5) + W[j][i] * (molar_fractions[j] - molar_fractions[i] - 0.5))
            RTlny[l] = val
        return RTlny

    def _non_ideal_interactions(self, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        Eint = self._non_ideal_function(self.We, molar_fractions)
        Sint = self._non_ideal_function(self.Ws, molar_fractions)
        Vint = self._non_ideal_function(self.Wv, molar_fractions)
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

    def excess_volume(self, pressure, temperature, molar_fractions):
        V_excess = np.dot(
            molar_fractions, self._non_ideal_function(self.Wv, molar_fractions))
        return V_excess

    def excess_entropy(self, pressure, temperature, molar_fractions):
        S_conf = -constants.gas_constant * \
            np.dot(IdealSolution._log_ideal_activities(
                self, molar_fractions), molar_fractions)
        S_excess = np.dot(
            molar_fractions, self._non_ideal_function(self.Ws, molar_fractions))
        return S_conf + S_excess

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        E_excess = np.dot(
            molar_fractions, self._non_ideal_function(self.We, molar_fractions))
        return E_excess + pressure * self.excess_volume(pressure, temperature, molar_fractions)

    def activity_coefficients(self, pressure, temperature, molar_fractions):
        if temperature > 1.e-10:
            return np.exp(self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions) / (constants.gas_constant * temperature))
        else:
            raise Exception("Activity coefficients not defined at 0 K.")

    def activities(self, pressure, temperature, molar_fractions):
        return IdealSolution.activities(self, pressure, temperature, molar_fractions) * self.activity_coefficients(pressure, temperature, molar_fractions)
