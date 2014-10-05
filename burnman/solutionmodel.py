# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings
import burnman
from burnman.processchemistry import *

R = 8.31446 # J/K/mol
kd = lambda x,y : 1 if x==y else 0
class SolutionModel:
    """
    This is the base class for a solution model,  intended for use
    in defining solid solutions and performing thermodynamic calculations
    on them.  

    A user wanting a new solution model should define the functions for
    excess_gibbs_free_energy and excess_volume.  In the base class these
    return zero, so if this class is used then the Gibbs free energy and molar
    volume of a solution are just the weighted arithmetic averages of the
    different endmember values.
    """
    
    def excess_gibbs_free_energy( self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Gibbs free energy of the solution,
        relative to an ideal model.
 
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
        return 0.0

    def excess_volume( self, pressure, temperature, molar_fractions):
        """
        Given a list of molar fractions of different phases,
        compute the excess Gibbs free energy of the solution,
        relative to an ideal model.  The base class implementation
        assumes that the excess volume is zero.
 
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


class IdealSolution ( SolutionModel ):
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
            ProcessSolidSolutionChemistry(self.formulas)

        self.endmember_configurational_entropies=np.zeros(shape=(self.n_endmembers))
        for idx, endmember_occupancy in enumerate(self.endmember_occupancies):
            for occ in range(self.n_occupancies):
                if endmember_occupancy[occ] > 1e-10: 
                    self.endmember_configurational_entropies[idx] = \
                        self.endmember_configurational_entropies[idx] - \
                        R*self.site_multiplicities[occ]*endmember_occupancy[occ]*np.log(endmember_occupancy[occ])



    def excess_gibbs_free_energy( self, pressure, temperature, molar_fractions ):
        return self.ideal_gibbs_excess( temperature, molar_fractions )

    def configurational_entropy (self, molar_fractions):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        conf_entropy=0
        for idx, occupancy in enumerate(site_occupancies):
            if occupancy > 1e-10:
                conf_entropy=conf_entropy-R*occupancy*self.site_multiplicities[idx]*np.log(occupancy)

        # Alternative entropy construction
        #activities=self.ideal_activities(molar_fractions)
        #conf_entropy=0.
        #for i, activity in enumerate(activities):
        #    if activity > 1e-10:
        #        conf_entropy=conf_entropy - R*molar_fractions[i]*np.log(activity) + molar_fractions[i]*self.endmember_configurational_entropies[i]

        return conf_entropy

    def endmember_configurational_entropy_contribution(self, molar_fractions):
        return np.dot(molar_fractions, self.endmember_configurational_entropies)

    def ideal_gibbs_excess( self, temperature, molar_fractions ): 
        return 0.0-temperature*(self.configurational_entropy(molar_fractions)  - self.endmember_configurational_entropy_contribution(molar_fractions))

    def ln_ideal_activities ( self, molar_fractions ):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        lna=np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            lna[e]=0.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] > 1e-10 and site_occupancies[occ] > 1e-10:
                    lna[e]=lna[e] + self.endmember_occupancies[e][occ]*self.site_multiplicities[occ]*np.log(site_occupancies[occ])

            normalisation_constant=self.endmember_configurational_entropies[e]/R
            lna[e]=lna[e] + self.endmember_configurational_entropies[e]/R
        return lna


    def ideal_activities ( self, molar_fractions ):
        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        activities=np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            activities[e]=1.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] > 1e-10:
                    activities[e]=activities[e]*np.power(site_occupancies[occ],self.endmember_occupancies[e][occ]*self.site_multiplicities[occ])
            normalisation_constant=np.exp(self.endmember_configurational_entropies[e]/R)
            activities[e]=normalisation_constant*activities[e]
        return activities

 
class AsymmetricRegularSolution ( IdealSolution ):
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
        IdealSolution.__init__(self, endmembers )
        
    def configurational_entropy( self, molar_fractions ):
        return IdealSolution.configurational_entropy( self, molar_fractions )

    def non_ideal_interactions( self, molar_fractions ):
        # -sum(sum(qi.qj.Wij*)
        # equation (2) of Holland and Powell 2003
        phi=np.array([self.alpha[i]*molar_fractions[i] for i in range(self.n_endmembers)])
        phi=np.divide(phi, np.sum(phi))


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


    def excess_gibbs_free_energy( self, pressure, temperature, molar_fractions ):

        ideal_gibbs = IdealSolution.ideal_gibbs_excess( self, temperature, molar_fractions )


        phi=np.array([self.alpha[i]*molar_fractions[i] for i in range(self.n_endmembers)])
        phi=np.divide(phi, np.sum(phi))

        H_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Wh,phi))
        S_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Ws,phi))
        V_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Wv,phi))

        non_ideal_gibbs = H_excess - temperature*S_excess + pressure*V_excess
     
        return ideal_gibbs + non_ideal_gibbs


    def excess_volume ( self, pressure, temperature, molar_fractions ):

        phi=np.array([self.alpha[i]*molar_fractions[i] for i in range(self.n_endmembers)])
        phi=np.divide(phi, np.sum(phi))

        V_excess=np.dot(self.alpha.T,molar_fractions)*np.dot(phi.T,np.dot(self.Wv,phi))
        return V_excess



class SymmetricRegularSolution ( AsymmetricRegularSolution ):
    """
    Solution model implementing the symmetric regular solution model
    """
    def __init__( self, endmembers, enthalpy_interaction, volume_interaction = None, entropy_interaction = None ):
        alphas = np.ones( len(endmembers) )
        AsymmetricRegularSolution.__init__(self, endmembers, alphas, enthalpy_interaction, volume_interaction, entropy_interaction )

