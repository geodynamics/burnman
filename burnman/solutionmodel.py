# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import

import numpy as np
import warnings
from .processchemistry import *
from . import constants

"""
kronecker delta function for integers
"""
kd = lambda x, y: 1 if x == y else 0


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
        return np.dot(np.array(molar_fractions), self.excess_partial_gibbs_free_energies( pressure, temperature, molar_fractions))

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
        return np.empty_like( np.array(molar_fractions) )

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

    def excess_dVdP(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess d(volume)/d(pressure) of the solution.
        The base class implementation assumes that this excess is zero.

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
        dVdP_excess : float
            The pressure gradient of the volume excess of the solution
        """
        return 0.0

    def excess_dVdT(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess d(volume)/d(temperature) of the solution.
        The base class implementation assumes that this excess is zero.

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
        dVdT_excess : float
            The temperature gradient of the volume excess of the solution
        """
        return 0.0

    def excess_dSdT(self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess d(entropy)/d(temperature) of the solution.
        The base class implementation assumes that this excess is zero.

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
        dSdT_excess : float
            The temperature gradient of the entropy excess of the solution
        """
        return 0.0

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
            process_solution_chemistry(self.formulas)

        self._calculate_endmember_configurational_entropies()

    def excess_partial_gibbs_free_energies( self, pressure, temperature, molar_fractions ):
        return self._ideal_excess_partial_gibbs( temperature, molar_fractions )

    def _calculate_endmember_configurational_entropies(self):
        self.endmember_configurational_entropies=np.zeros(shape=(self.n_endmembers))
        for idx, endmember_occupancy in enumerate(self.endmember_occupancies):
            for occ in range(self.n_occupancies):
                if endmember_occupancy[occ] > 1e-10:
                    self.endmember_configurational_entropies[idx] = \
                        self.endmember_configurational_entropies[idx] - \
                        constants.gas_constant*self.site_multiplicities[occ]*endmember_occupancy[occ]*np.log(endmember_occupancy[occ])

    def _endmember_configurational_entropy_contribution(self, molar_fractions):
        return np.dot(molar_fractions, self.endmember_configurational_entropies)

    def _configurational_entropy (self, molar_fractions):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        conf_entropy=0
        for idx, occupancy in enumerate(site_occupancies):
            if occupancy > 1e-10:
                conf_entropy=conf_entropy-constants.gas_constant*occupancy*self.site_multiplicities[idx]*np.log(occupancy)

        return conf_entropy


    def _ideal_excess_partial_gibbs(self, temperature, molar_fractions):
        return  constants.gas_constant*temperature * self._log_ideal_activities(molar_fractions)

    def _log_ideal_activities (self, molar_fractions):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        lna=np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            lna[e]=0.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] > 1e-10 and site_occupancies[occ] > 1e-10:
                    lna[e]=lna[e] + self.endmember_occupancies[e][occ]*self.site_multiplicities[occ]*np.log(site_occupancies[occ])

            normalisation_constant=self.endmember_configurational_entropies[e]/constants.gas_constant
            lna[e]=lna[e] + self.endmember_configurational_entropies[e]/constants.gas_constant
        return lna

    def _ideal_activities (self, molar_fractions):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        activities=np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            activities[e]=1.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] > 1e-10:
                    activities[e]=activities[e]*np.power(site_occupancies[occ],self.endmember_occupancies[e][occ]*self.site_multiplicities[occ])
            normalisation_constant=np.exp(self.endmember_configurational_entropies[e]/constants.gas_constant)
            activities[e]=normalisation_constant*activities[e]
        return activities

 
class AsymmetricRegularSolution (IdealSolution):
    """
    Solution model implementing the asymmetric regular solution model formulation (Holland and Powell, 2003)
    """

    def __init__( self, endmembers, alphas, enthalpy_interaction, volume_interaction = None, entropy_interaction = None ): 

        self.n_endmembers = len(endmembers)

        # Create array of van Laar parameters
        self.alpha=np.array(alphas)

        # Create 2D arrays of interaction parameters
        self.Wh=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Ws=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wv=np.zeros(shape=(self.n_endmembers,self.n_endmembers))

        #setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i+1, self.n_endmembers):
                self.Wh[i][j]=2.*enthalpy_interaction[i][j-i-1]/(self.alpha[i]+self.alpha[j])

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i+1, self.n_endmembers):
                    self.Ws[i][j]=2.*entropy_interaction[i][j-i-1]/(self.alpha[i]+self.alpha[j])

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i+1, self.n_endmembers):
                    self.Wv[i][j]=2.*volume_interaction[i][j-i-1]/(self.alpha[i]+self.alpha[j])

        #initialize ideal solution model
        IdealSolution.__init__(self, endmembers)
        
    def _phi( self, molar_fractions):
        phi=np.array([self.alpha[i]*molar_fractions[i] for i in range(self.n_endmembers)])
        phi=np.divide(phi, np.sum(phi))
        return phi

    def _non_ideal_interactions(self, molar_fractions):
        # -sum(sum(qi.qj.Wij*)
        # equation (2) of Holland and Powell 2003
        phi=self._phi(molar_fractions)

        q=np.zeros(len(molar_fractions))
        Hint=np.zeros(len(molar_fractions))
        Sint=np.zeros(len(molar_fractions))
        Vint=np.zeros(len(molar_fractions))

        for l in range(self.n_endmembers):
            q=np.array([kd(i,l)-phi[i] for i in range(self.n_endmembers)])

            Hint[l]=0.-self.alpha[l]*np.dot(q,np.dot(self.Wh,q))
            Sint[l]=0.-self.alpha[l]*np.dot(q,np.dot(self.Ws,q))
            Vint[l]=0.-self.alpha[l]*np.dot(q,np.dot(self.Wv,q))
     
        return Hint, Sint, Vint

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions) :
        Hint, Sint, Vint = self._non_ideal_interactions( molar_fractions )
        return Hint - temperature*Sint + pressure*Vint

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs (self, temperature, molar_fractions )
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_volume (self, pressure, temperature, molar_fractions):
        phi=self._phi(molar_fractions)
        V_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Wv,phi))
        return V_excess

    def excess_entropy(self, pressure, temperature, molar_fractions):
        phi=self._phi(molar_fractions)
        S_conf=-constants.gas_constant*np.dot(IdealSolution._log_ideal_activities(self, molar_fractions), molar_fractions)
        S_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Ws,phi))
        return S_conf + S_excess

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        phi=self._phi(molar_fractions)
        H_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Wh,phi))
        return H_excess + pressure*self.excess_volume ( pressure, temperature, molar_fractions )


class SymmetricRegularSolution (AsymmetricRegularSolution):
    """
    Solution model implementing the symmetric regular solution model
    """
    def __init__( self, endmembers, enthalpy_interaction, volume_interaction = None, entropy_interaction = None ):
        alphas = np.ones( len(endmembers) )
        AsymmetricRegularSolution.__init__(self, endmembers, alphas, enthalpy_interaction, volume_interaction, entropy_interaction )

class SubregularSolution (IdealSolution):
    """
    Solution model implementing the subregular solution model formulation (Helffrich and Wood, 1989)
    """

    def __init__( self, endmembers, enthalpy_interaction, volume_interaction = None, entropy_interaction = None ): 

        self.n_endmembers = len(endmembers)

        # Create 2D arrays of interaction parameters
        self.Wh=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Ws=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wv=np.zeros(shape=(self.n_endmembers,self.n_endmembers))

        #setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i+1, self.n_endmembers):
                self.Wh[i][j]=enthalpy_interaction[i][j-i-1][0]
                self.Wh[j][i]=enthalpy_interaction[i][j-i-1][1]

        if entropy_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i+1, self.n_endmembers):
                    self.Ws[i][j]=entropy_interaction[i][j-i-1][0]
                    self.Ws[j][i]=entropy_interaction[i][j-i-1][1]

        if volume_interaction is not None:
            for i in range(self.n_endmembers):
                for j in range(i+1, self.n_endmembers):
                    self.Wv[i][j]=volume_interaction[i][j-i-1][0]
                    self.Wv[j][i]=volume_interaction[i][j-i-1][1]

        #initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def _non_ideal_function(self, W, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        n=len(molar_fractions)
        RTlny=np.zeros(n)
        for l in range(n):
            val=0.
            for i in range(n):
                if i != l:
                    val+=0.5*molar_fractions[i]*(W[l][i]*(1-molar_fractions[l] + molar_fractions[i] + 2.*molar_fractions[l]*(molar_fractions[l] - molar_fractions[i] - 1)) + W[i][l]*(1.-molar_fractions[l] - molar_fractions[i] - 2.*molar_fractions[l]*(molar_fractions[l] - molar_fractions[i] - 1)))
                    for j in range(i+1,n):
                        if j != l:
                            val+=molar_fractions[i]*molar_fractions[j]*(W[i][j]*(molar_fractions[i] - molar_fractions[j] - 0.5) + W[j][i]*(molar_fractions[j] - molar_fractions[i] - 0.5))
            RTlny[l] = val
        return RTlny

    def _non_ideal_interactions(self, molar_fractions):
        # equation (6') of Helffrich and Wood, 1989
        Hint=self._non_ideal_function(self.Wh, molar_fractions)
        Sint=self._non_ideal_function(self.Ws, molar_fractions)
        Vint=self._non_ideal_function(self.Wv, molar_fractions)
        return Hint, Sint, Vint

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions) :
        Hint, Sint, Vint = self._non_ideal_interactions(molar_fractions)
        return Hint - temperature*Sint + pressure*Vint

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs (self, temperature, molar_fractions)
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_volume (self, pressure, temperature, molar_fractions):
        V_excess=np.dot(molar_fractions, self._non_ideal_function(self.Wv, molar_fractions))
        return V_excess

    def excess_entropy(self, pressure, temperature, molar_fractions):
        S_conf=-constants.gas_constant*np.dot(IdealSolution._log_ideal_activities(self, molar_fractions), molar_fractions)
        S_excess=np.dot(molar_fractions, self._non_ideal_function(self.Ws, molar_fractions))
        return S_conf + S_excess

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        H_excess=np.dot(molar_fractions, self._non_ideal_function(self.Wh, molar_fractions))
        return H_excess + pressure*self.excess_volume (pressure, temperature, molar_fractions)

class FullSubregularSolution (SubregularSolution):
    """
    Solution model implementing the subregular solution model formulation (Helffrich and Wood, 1989)
    Intermediates are defined to describe each of the interaction terms as a function of P and T
    """

    def __init__( self, endmembers, intermediates): 
        self.endmembers = endmembers
        self.intermediates = intermediates
        #initialize ideal solution model
        IdealSolution.__init__(self, endmembers)

    def set_interaction_terms(self, pressure, temperature):
        for i in range(self.n_endmembers-1):
            for j in range(self.n_endmembers-i-1):
                self.intermediates[i][j][0].set_state(pressure, temperature)
                self.intermediates[i][j][1].set_state(pressure, temperature)


        # Create 2D arrays of interaction parameters
        self.enthalpy_interaction = [[[[] for k in xrange(2)] for j in xrange(self.n_endmembers - i - 1)] for i in xrange(self.n_endmembers - 1)]
        self.entropy_interaction = [[[[] for k in xrange(2)] for j in xrange(self.n_endmembers - i - 1)] for i in xrange(self.n_endmembers - 1)]
        self.volume_interaction = [[[[] for k in xrange(2)] for j in xrange(self.n_endmembers - i - 1)] for i in xrange(self.n_endmembers - 1)]

        self.Wg=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wh=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Ws=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wv=np.zeros(shape=(self.n_endmembers,self.n_endmembers))

        self.Wdvdp=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wdvdt=np.zeros(shape=(self.n_endmembers,self.n_endmembers))
        self.Wdsdt=np.zeros(shape=(self.n_endmembers,self.n_endmembers))

        #setup excess enthalpy interaction matrix
        for i in range(self.n_endmembers):
            for j in range(i+1, self.n_endmembers):
                fractions = []
                for mbr in range(self.n_endmembers):
                    if mbr == i or mbr == j:
                        fractions.append(0.5)
                    else:
                        fractions.append(0.0)
                Sconf = -constants.gas_constant*np.dot(IdealSolution._log_ideal_activities(self, fractions), fractions)


                self.Wg[i][j]=4.*(self.intermediates[i][j-i-1][0].gibbs + temperature*Sconf) \
                               - 2.*(self.endmembers[i][0].gibbs + self.endmembers[j][0].gibbs)
                self.Wg[j][i]=4.*(self.intermediates[i][j-i-1][1].gibbs + temperature*Sconf) \
                               - 2.*(self.endmembers[i][0].gibbs + self.endmembers[j][0].gibbs) 

                self.Ws[i][j] = 4.*(self.intermediates[i][j-i-1][0].S - Sconf) \
                               - 2.*(self.endmembers[i][0].S + self.endmembers[j][0].S) 
                self.Ws[j][i] = 4.*(self.intermediates[i][j-i-1][1].S - Sconf) \
                                - 2.*(self.endmembers[i][0].S + self.endmembers[j][0].S)
                self.entropy_interaction[i][j-i-1][0] = self.Ws[i][j]
                self.entropy_interaction[i][j-i-1][1] = self.Ws[j][i]

                self.Wv[i][j] = 4.*(self.intermediates[i][j-i-1][0].V) \
                                - 2.*(self.endmembers[i][0].V + self.endmembers[j][0].V)
                self.Wv[j][i] = 4.*(self.intermediates[i][j-i-1][1].V) \
                                - 2.*(self.endmembers[i][0].V + self.endmembers[j][0].V)
                self.volume_interaction[i][j-i-1][0] = self.Wv[i][j]
                self.volume_interaction[i][j-i-1][1] = self.Wv[j][i]

                self.Wdvdp[i][j] = -4.*(self.intermediates[i][j-i-1][0].V \
                                        / self.intermediates[i][j-i-1][0].K_T) \
                    + 2.*(self.endmembers[i][0].V/self.endmembers[i][0].K_T \
                          + self.endmembers[j][0].V/self.endmembers[j][0].K_T)
                self.Wdvdp[j][i] = -4.*(self.intermediates[i][j-i-1][1].V \
                                        / self.intermediates[i][j-i-1][1].K_T) \
                    + 2.*(self.endmembers[i][0].V/self.endmembers[i][0].K_T \
                          + self.endmembers[j][0].V/self.endmembers[j][0].K_T)
                
                self.Wdvdt[i][j] = 4.*(self.intermediates[i][j-i-1][0].alpha \
                                       * self.intermediates[i][j-i-1][0].V) \
                    - 2.*(self.endmembers[i][0].alpha*self.endmembers[i][0].V \
                          + self.endmembers[j][0].alpha*self.endmembers[j][0].V)        
                self.Wdvdt[j][i] = 4.*(self.intermediates[i][j-i-1][1].alpha \
                                       * self.intermediates[i][j-i-1][1].V) \
                    - 2.*(self.endmembers[i][0].alpha*self.endmembers[i][0].V \
                          + self.endmembers[j][0].alpha*self.endmembers[j][0].V)

                self.Wdsdt[i][j] = (1./temperature)*(4.*(self.intermediates[i][j-i-1][0].C_p)
                                    - 2.*(self.endmembers[i][0].C_p + self.endmembers[j][0].C_p))
                self.Wdsdt[j][i] = (1./temperature)*(4.*(self.intermediates[i][j-i-1][1].C_p)
                                    - 2.*(self.endmembers[i][0].C_p + self.endmembers[j][0].C_p))

    def _non_ideal_function(self, W, molar_fractions):
        return SubregularSolution._non_ideal_function(self, W, molar_fractions )

    def _non_ideal_excess_partial_gibbs(self, pressure, temperature, molar_fractions) :
        return self._non_ideal_function(self.Wg, molar_fractions)

    def excess_partial_gibbs_free_energies(self, pressure, temperature, molar_fractions):
        ideal_gibbs = IdealSolution._ideal_excess_partial_gibbs (self, temperature, molar_fractions)
        non_ideal_gibbs = self._non_ideal_excess_partial_gibbs(pressure, temperature, molar_fractions)
        return ideal_gibbs + non_ideal_gibbs

    def excess_volume (self, pressure, temperature, molar_fractions):
        return SubregularSolution.excess_volume(self, pressure, temperature, molar_fractions )

    def excess_entropy(self, pressure, temperature, molar_fractions):
        return SubregularSolution.excess_entropy(self, pressure, temperature, molar_fractions )

    def excess_enthalpy(self, pressure, temperature, molar_fractions):
        return self.excess_gibbs_free_energy(pressure, temperature, molar_fractions) \
            + temperature*self.excess_entropy(pressure, temperature, molar_fractions)

    def excess_dVdP (self, pressure, temperature, molar_fractions):
        dVdP_excess=np.dot(molar_fractions, self._non_ideal_function(self.Wdvdp, molar_fractions))
        return dVdP_excess

    def excess_dVdT (self, pressure, temperature, molar_fractions):
        dVdT_excess=np.dot(molar_fractions, self._non_ideal_function(self.Wdvdt, molar_fractions))
        return dVdT_excess

    def excess_dSdT (self, pressure, temperature, molar_fractions):
        dSdT_excess=np.dot(molar_fractions, self._non_ideal_function(self.Wdsdt, molar_fractions))
        return dSdT_excess
