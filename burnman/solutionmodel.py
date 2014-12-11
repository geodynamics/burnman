# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import warnings
import burnman
from burnman.processchemistry import *

R = 8.3145 # J/K/mol

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

    def excess_gibbs_free_energy( self, pressure, temperature, molar_fractions ):
        return self.ideal_gibbs_excess( temperature, molar_fractions )

    def ideal_gibbs_excess ( self, temperature, molar_fractions ): 

        activities = self.ideal_activities( molar_fractions )
        gibbs_excess_ideal = 0.0

        tol = 1.e-10
        for i in range(self.n_endmembers):
            if molar_fractions[i] > tol:
                gibbs_excess_ideal = gibbs_excess_ideal +  \
                                     molar_fractions[i] * R * temperature * np.log(activities[i])
        
        return gibbs_excess_ideal

    def ideal_activities ( self, molar_fractions ):

        site_occupancies=np.dot(molar_fractions, self.endmember_occupancies)
        activities=np.empty(shape=(self.n_endmembers))

        for e in range(self.n_endmembers):
            activities[e]=1.0
            normalisation_constant=1.0
            for occ in range(self.n_occupancies):
                if self.endmember_occupancies[e][occ] != 0: #integer?  
                    activities[e]=activities[e]*np.power(site_occupancies[occ],self.site_multiplicities[occ])
                    normalisation_constant=normalisation_constant/np.power(self.endmember_occupancies[e][occ],self.site_multiplicities[occ])
            activities[e]=normalisation_constant*activities[e]

        return activities

 
class AsymmetricRegularSolution ( IdealSolution ):
    """
    Solution model implementing the Asymmetric Van Laar formulation
    (best reference?).
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
    Solution model implementing the Symmetric Van Laar solution
    """
    def __init__( self, endmembers, enthalpy_interaction, volume_interaction = None, entropy_ineraction = None ):
        #symmetric case if all the alphas are equal?
        alphas = np.ones( len(endmembers) )
        AsymmetricRegularSolution.__init__(self, endmembers, alphas, enthalpy_interaction, volume_ineteraction, entropy_interaction )

