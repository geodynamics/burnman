# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt

import burnman
import burnman.gibbsminimization as gm


class EquilibriumAssemblage(burnman.Material):
    """
    Class for taking a certain assemblage of elements,
    then calculating the equilibrium assemblage of minerals
    that minimizes the Gibbs free energy.  This should
    have similar functionality to :class:`burnman.Composite`,
    but instead of having a fixed set of minerals, it 
    dynamically updates the minerals when you call set_state().
    As such, it should be significantly slower to use.
    """ 

    def __init__(self, composition, phases):
        """
        Initialize the equilibrium assemblage with the elements which
        comprise it as well as a list of phases into which those 
        elements may go.
  
        Parameters
        ----------
        composition : dictionary
            Dictionary where the keys are strings corresponding to 
            elements (e.g. 'Mg') and the values are molar fractions
            of that element.
        phases : list of :class:`burnman.Mineral` or :class:`burnman.SolidSolution`
            List of phases over which the equilibrium assemblage will try
            to minimize the gibbs free energy. This class only understands
            how to use instances of :class:`burnman.Mineral` or 
            :class:`burnman.SolidSolution`. 
        """
        self.composition = composition
        self.phases = phases

        # Calculate the stoichiometric matrix
        stoich, stoich_el, formulae = gm.assemble_stoichiometric_matrix(phases)
        #The elements in the composition object should be a subset of the
        #set of elements in the various phases
        assert( set(composition.keys()).issubset(set(stoich_el) ))

        self.endmember_formulae = formulae
        self.stoichiometric_matrix = stoich
        self.elements = stoich_el

        #There may now be elements in the stoichiometric matrix that were not in the initial composition
        #dictionary because they appeared in one or more of the minerals.  We need to account for this
        #when we make the bulk composition vector
        self.bulk_composition_vector = np.array([ (composition[e] if e in composition.keys() else 0.0) \
                                                   for e in self.elements] )
 
        #Calculate the nullspace and the baseline assemblage vector
        self.__setup_subspaces()
        self.__compute_baseline_assemblage()


    def set_method(self, method):
        """
        set the same equation of state method for all the phases in the assemblage
        """
        for phase in self.phases:
            phase.set_method(method)


    def set_state( self, pressure, temperature):
        """
        Set the state of the assemblage, as a function of pressure and temperature.  
        The performs a minimization of the Gibbs free energy to determine the equilibrium
        assemblage, and as such, is not a cheap operation. 
  
        Parameters
        ----------
        pressure : float
            The pressure in [Pa] at which to set the state. 
        temperature : float
            The temperature in [K] at which to set the state. 
        """
        self._minimize_gibbs( pressure, temperature )
#        self._solve_equilibrium_equations (pressure, temperature )


    def _solve_equilibrium_equations( self, pressure, temperature ): 
        
        ref_gibbs = self.__compute_gibbs( pressure, temperature, self.__compute_species_vector(self.reduced_species_vector * 0.))
        print ref_gibbs
        equilibrium = lambda x : (self._equilibrium_equations( pressure, temperature, self.__compute_species_vector( x ) ))
        sol = opt.root(equilibrium, self.reduced_species_vector, method='hybr', options={'xtol':1.e-11})

        self.reduced_species_vector = sol.x
        self.species_vector = self.__compute_species_vector(self.reduced_species_vector)
        self.gibbs = self.__compute_gibbs(pressure, temperature, self.species_vector)

    def _equilibrium_equations ( self, pressure, temperature, species_vector):
        partial_gibbs = self._compute_partial_gibbs( pressure, temperature, species_vector)
        gibbs_inequalities = np.dot(partial_gibbs, self.right_nullspace)

        if np.any(species_vector < -1.e-6):
           gibbs_inequalitites = np.ones_like( gibbs_inequalities) * 10000000.

        return gibbs_inequalities

    def _minimize_gibbs( self, pressure, temperature):

        n = len(self.endmember_formulae)
        
        # Calculate a reference gibbs free energy at the baseline assemblage.  This is a kind of ill-conditioned
        # minimization problem, and it seems to work better if we subtract off a reference state and normalize.
        ref_gibbs = self.__compute_gibbs( pressure, temperature, self.__compute_species_vector(self.reduced_species_vector * 0.))
 
        #define the function to minimize and then do it.
        minimize_gibbs = lambda x : (self.__compute_gibbs( pressure, temperature, self.__compute_species_vector(x) ) - ref_gibbs)/np.abs(ref_gibbs)
#        minimize_gibbs = lambda x : (self.__compute_gibbs( pressure, temperature, self.__compute_species_vector(x) ))
        equilibrium = lambda x : (self._equilibrium_equations( pressure, temperature, self.__compute_species_vector( x ) )/np.abs(ref_gibbs))
#        equilibrium=None
    
        #Set the solution
        constraints={'type':'ineq', 'fun':self.__compute_species_vector}
        sol = opt.minimize( minimize_gibbs, self.reduced_species_vector, method='SLSQP', jac=equilibrium, constraints=constraints )
        self.reduced_species_vector = sol.x
        self.species_vector = self.__compute_species_vector(self.reduced_species_vector)
        self.gibbs = self.__compute_gibbs(pressure, temperature, self.species_vector)
        



    def print_assemblage(self):
        """
        Print the current abundance of each endmember in the assemblage, as a molar fraction.
        """
        tot = np.sum(self.species_vector)
        for f,s in zip(self.endmember_formulae, self.species_vector):
            print f, s/tot

    def _compute_partial_gibbs( self, pressure, temperature, species_vector):
        partial_gibbs = np.empty_like(species_vector)

        #Loop over the various phases and compute the gibbs free energy
        #of each one at P,T.  We have to treat solid solutions and
        #single phases somewhat differently.
        i = 0
        for phase in self.phases:
            if isinstance (phase, burnman.SolidSolution):

                n = len(phase.base_material)
                frac = np.sum(species_vector[i:(i+n)])
                molar_fractions = ( species_vector[i:(i+n)]/frac if (frac > 1.e-6)  else np.ones( n )/n )
                phase.set_composition( molar_fractions )
                phase.set_state( pressure, temperature )
                
                partial_gibbs[i:(i+n)] = phase.partial_gibbs
                i+=n
                   
            elif isinstance(phase, burnman.Mineral):
                phase.set_state( pressure, temperature )
                partial_gibbs[i] = phase.gibbs
                i+=1
            else:
                raise Exception('Unsupported mineral type, can only read burnman.Mineral or burnman.SolidSolution')

        return partial_gibbs

 
    def __compute_gibbs( self, pressure, temperature, species_vector ):
        """
        Given a pressure, temperature, and vector in the nullspace, 
        calculate the gibbs free energy of the assemblage.  This 
        is basically the function to minimize when taking the 
        assemblage to a P-T, subject to the bulk composition contraint
        (which is parameterized by vectors in the nullspace)
        """

        assert( len(species_vector) == len(self.endmember_formulae) )

        gibbs = np.dot( species_vector, self._compute_partial_gibbs(pressure, temperature, species_vector) )
        return gibbs


    def __compute_species_vector ( self, reduced_vector ):
        """
        Given a vector in the nullspace, return a full species vector
        in the endmember space.
        """
        species_vector = self.baseline_assemblage + np.dot( self.right_nullspace, np.transpose(reduced_vector) )
        return species_vector


    def __compute_bulk_composition (self, species_vector):
        """
        Given a vector in the endmember space, return the 
        bulk composition.
        """
        return np.dot(self.stoichiometric_matrix, species_vector)


    def __setup_subspaces (self):
        """
        Calculate the nullspace for the bulk composition constraint, and make 
        an attempt to sparsify it.
        """
        self.right_nullspace = gm.compute_nullspace( self.stoichiometric_matrix )
        self.right_nullspace = gm.sparsify_basis(self.right_nullspace)


    def __compute_baseline_assemblage(self):
        """
        The nullspace gives us vectors that do not change the bulk composition, 
        but we need a particular vector that satisfies our composition constraint,
        so that general vectors may be represented with the particular vector plus
        excursions in the nullspace.  This calculates our particular vector, called
        the baseline assemblage.
        """

        # Initially I calculated this using the SVD and projecting onto the nullspace,
        # but it was more complicated and less robust, and the non-negativity constraint
        # was kind of tricky.  Non-negative least squares does the same thing, is simpler,
        # and is more robust. [IR]
        eps = 1.e-10
        baseline_assemblage = opt.nnls( self.stoichiometric_matrix, self.bulk_composition_vector)

        # It is possible, even easy, to not be able to represent a given composition by
        # a given set of phases.  Raise an exception if this occurs.
        if  baseline_assemblage[1] > eps :
            raise Exception( "Composition cannot be represented by the given minerals." )

        assert( np.all(np.abs(np.dot(self.stoichiometric_matrix, baseline_assemblage[0])\
                                                     - self.bulk_composition_vector) < eps) )
          
        #Set the baseline assemblage
        self.baseline_assemblage = baseline_assemblage[0]

        #Our starting vector in the nullspace is just going to be the zero vector
        self.reduced_species_vector = np.zeros( self.right_nullspace.shape[1] )
        self.species_vector = self.__compute_species_vector( self.reduced_species_vector)
