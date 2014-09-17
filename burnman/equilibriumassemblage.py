# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import warnings

import burnman
import burnman.gibbsminimization as gm


class EquilibriumAssemblage(burnman.Material):
    """
    Class for taking a certain assemblage of elements,
    then calculating the equilibrium assemblage of minerals
    that minimizes the Gibbs free energy.  This should
    have similar functionality to ``burnman.Composite'',
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
        elements : dictionary
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

    def set_method(self, method):
        for phase in self.phases:
            phase.set_method(method)

 
    def set_state( self, pressure, temperature):
        
        n = len(self.endmember_formulae)
        minimize_gibbs = lambda x : self.__compute_gibbs( pressure, temperature, x )
        sol = opt.fmin_slsqp( minimize_gibbs, np.ones(n)/n, f_eqcons = self.__composition_constraint, bounds=[(0.0, 1.0),]*n, full_output=0)
        self.species_vector = sol
        print zip(self.endmember_formulae, self.species_vector)

    def __compute_gibbs( self, pressure, temperature, species_vector ):

        assert( len(species_vector) == len(self.endmember_formulae) )

        tmp_gibbs = 0.0
        i = 0

        for phase in self.phases:
            if isinstance (phase, burnman.SolidSolution):
                n = len(phase.base_material)
                total_frac = np.sum( species_vector[i:(i+n)] )
                phase.set_method('slb3')
                phase.set_composition( np.array( species_vector[i:(i+n)]/total_frac) )
                phase.set_state( pressure, temperature )
                tmp_gibbs += phase.gibbs * total_frac
                i+=n
            elif isinstance(phase, burnman.Mineral):
                phase.set_method('slb3')
                phase.set_state( pressure, temperature )
                tmp_gibbs += phase.gibbs * species_vector[i]
                i+=1
            else:
                raise Exception('Unsupported mineral type, can only read burnman.Mineral or burnman.SolidSolution')

        return tmp_gibbs
                
        
    def __composition_constraint (self, species_vector ):
        assert( len(species_vector) == len(self.endmember_formulae) )
        return np.dot( self.stoichiometric_matrix, species_vector) - self.bulk_composition_vector


