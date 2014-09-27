# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import scipy.linalg as linalg
import warnings
import matplotlib.pyplot as plt

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

        stoich, stoich_el, formulae = gm.assemble_stoichiometric_matrix(phases)
        #The elements in the composition object should be a subset of the
        #set of elements in the various phases
        assert( set(composition.keys()).issubset(set(stoich_el) ))

        self.endmember_formulae = formulae
        self.stoichiometric_matrix = stoich
        self.elements = stoich_el

        self.bulk_composition_vector = np.array([ (composition[e] if e in composition.keys() else 0.0) \
                                                   for e in self.elements] )
 
        self.__setup_subspaces()
        self.__compute_baseline_assemblage()


    def set_method(self, method):
        for phase in self.phases:
            phase.set_method(method)


    def set_state( self, pressure, temperature):
        
        n = len(self.endmember_formulae)
        ref_gibbs = self.__compute_gibbs( pressure, temperature, self.reduced_species_vector *0.)
        minimize_gibbs = lambda x : (self.__compute_gibbs( pressure, temperature, x ) - ref_gibbs)/np.abs(ref_gibbs)
        sol = opt.fmin_slsqp( minimize_gibbs, self.reduced_species_vector, f_ieqcons = lambda x : self.__compute_species_vector(x), full_output=1, disp=False)
        self.reduced_species_vector = sol[0]
        self.species_vector = self.__compute_species_vector(self.reduced_species_vector)
        self.gibbs = sol[1]*np.abs(ref_gibbs) + ref_gibbs


    def print_assemblage(self):
        tot = np.sum(self.species_vector)
        for f,s in zip(self.endmember_formulae, self.species_vector):
            print f, s/tot

 
    def __compute_gibbs( self, pressure, temperature, reduced_vector ):

        species_vector = self.__compute_species_vector(reduced_vector)
        assert( len(species_vector) == len(self.endmember_formulae) )

        tmp_gibbs = 0.0
        i = 0

        for phase in self.phases:
            if isinstance (phase, burnman.SolidSolution):
                n = len(phase.base_material)
                phase.set_method('slb3')
                frac = np.sum(species_vector[i:(i+n)])
                molar_fractions = ( species_vector[i:(i+n)]/frac if frac > 0.0 else np.ones( n )/n )
                phase.set_composition( molar_fractions )
                phase.set_state( pressure, temperature )
                tmp_gibbs += phase.gibbs
                i+=n
            elif isinstance(phase, burnman.Mineral):
                phase.set_method('slb3')
                phase.set_state( pressure, temperature )
                tmp_gibbs += phase.gibbs * species_vector[i]
                i+=1
            else:
                raise Exception('Unsupported mineral type, can only read burnman.Mineral or burnman.SolidSolution')

        return tmp_gibbs


    def __compute_species_vector ( self, reduced_vector ):
        return self.baseline_assemblage + np.dot( self.right_nullspace, np.transpose(reduced_vector) )


    def __compute_bulk_composition (self, species_vector):
        return np.dot(self.stoichiometric_matrix, species_vector)


    def __setup_subspaces (self):
 
        self.right_nullspace = gm.compute_nullspace( self.stoichiometric_matrix )
        self.right_nullspace = gm.sparsify_basis(self.right_nullspace)


    def __compute_baseline_assemblage(self):

        eps = 1.e-10
        baseline_assemblage = opt.nnls( self.stoichiometric_matrix, self.bulk_composition_vector)
        if  baseline_assemblage[1] > eps :
            raise Exception( "Composition cannot be represented by the given minerals." )

        assert( np.all(np.abs(np.dot(self.stoichiometric_matrix, baseline_assemblage[0])\
                                                     - self.bulk_composition_vector) < eps) )
          
        self.baseline_assemblage = baseline_assemblage[0]
        self.reduced_species_vector = np.zeros( self.right_nullspace.shape[1] )
        self.species_vector = self.__compute_species_vector( self.reduced_species_vector)
