# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import scipy.linalg as linalg
import scipy.optimize as opt
from .mineral import Mineral
from .solidsolution import SolidSolution
import warnings

def assemble_stoichiometric_matrix ( minerals):
    """
    This takes a list of minerals and assembles a matrix where 
    the rows are elements and the columns are species (or endmembers).
    If a solid solution is passed in, then the endmembers are extracted
    from it. 
    
    Parameters
    ----------
    minerals : list of minerals
        List of objects of type :class:`burnman.Mineral` or `burnman.SolidSolution`.
        Other types cannot be understood, and an exception will be thrown.
    
    Returns
    -------
    stoichiometric_matrix: 2D numpy array 
        2D numpy array that is n_elements by n_endmembers.  The matrix entries
        are the number of atoms of an element in that particular endmember.
        The elements are ordered alphabetically, and the endmembers are ordered
        in the order in which they are passed in.
    elements: list of strings
        List of elements that are in the stoichiometric matrix, in alphabetical
        order.
    formulae: list of strings
        List of endmember formulae constructed from the minerals passed in.
        They are ordered in the same order as they are passed in, but with 
        The endmembers in the solid solutions also included.
    """
 
    elements = set()
    formulae = []
    n_phase_endmembers = []
    
    # Make a list of the different formulae, as well as a set of 
    # the elements that we have present in the mineral list
    for m in minerals:
        # Add the endmembers if it is a solid solution
        if isinstance(m, SolidSolution):
            n_phase_endmembers.append(len(m.endmembers))
            for e in m.endmembers:
                f = e[0].params['formula']
                formulae.append(f)
                for k in f.keys():
                    elements.add(k)
        # Add formula if it is a simple mineral
        elif isinstance(m, Mineral):
            n_phase_endmembers.append(1)
            f = m.params['formula']
            formulae.append(f)
            for k in f.keys():
                elements.add(k)
        else:
            raise Exception('Unsupported mineral type, can only read burnman.Mineral or burnman.SolidSolution')

    #Listify the elements and sort them so they have a well-defined order.
    #This will be the ordering of the rows.  The ordering of the columns
    #will be the ordering of the endmembers as they are passed in.
    elements = list(elements)
    elements.sort()

    #Populate the stoichiometric matrix
    stoichiometric_matrix = np.empty( [ len(elements), len(formulae) ] )
    for i,e in enumerate(elements):
        for j,f in enumerate(formulae):
            stoichiometric_matrix[i,j] = ( f[e]  if e in f else 0.0 )

    return stoichiometric_matrix, elements, formulae, n_phase_endmembers


def compute_column_and_null_spaces ( stoichiometric_matrix ):
    """
    Given a stoichiometric matrix, compute a basis for the column and nullspaces.
    These bases corresponds to the subspaces that do and do not change the 
    bulk composition.  The vectors of the nullspace correspond to chemical
    reactions.  This merely identifies a set of linearly independent
    reactions, without saying anything about which ones are likely to
    occur.

    The calculation of the nullspace is done with SVD.
    
    Parameters
    ----------
    stoichiometric_matrix: numpy array
        The stoichiometric matrix, presumably calculated using :func:`compute_stoichiometric_matrix`.
    
    Returns
    -------
    col_space: numpy array
        A 2D array corresponding to the column space of the stoichiometric matrix, 
        with the columns corresponding to basis vectors.
    null_space: numpy array
        A 2D array corresponding to the nullspace of the stoichiometric matrix, 
        with the columns corresponding to basis vectors.
    """

    eps = 1.e-10
    # Do an SVD of the stoichiometric matrix.
    U, S, Vh = linalg.svd( stoichiometric_matrix)

    # The columns of V that correspond to large and small (or nonexistent)
    # singular values are the colume and left null spaces for the matrix.
    # select them, appropriately transpose them, and return them
    col_mask = ( np.append(S, np.zeros(len(Vh)-len(S))) > eps)
    col_space = np.compress(col_mask, Vh, axis=0)
    
    null_mask = ( np.append(S, np.zeros(len(Vh)-len(S))) <= eps)
    null_space = np.compress(null_mask, Vh, axis=0)
    
    return np.transpose(col_space), np.transpose(null_space)
    

def sparsify_basis ( basis ):
    """
    The vectors computed using :func:`compute_nullspace` are non-unique, and are
    not likely to correspond to the most physically meaningful set of reactions.  
    Instead we would often like to know a set of nullspace vectors that are the 
    "sparsest", meaning that they involve a minimal set of reactions.  This is known
    to be an NP-hard problem.  This function makes a fairly crude attempt to 
    take a set of vectors and rotate them to the sparsest coordinate system. It does 
    this by doing many L1 minimizations of the vectors, and as such might be 
    pretty slow, especially for large sets of bases.
    
    Parameters
    ----------
    basis: numpy array
        A 2D array correspnding to the basis, 
        with the columns corresponding to basis vectors.
   
    Returns
    -------
    new_basis: numpy array
        A 2D array correspnding to the attempt at a sparser basis, 
        with the columns corresponding to basis vectors.
    """
   
    eps = 1.e-6
    new_basis = basis.copy()
    n_cols = new_basis.shape[1]
    n_rows = new_basis.shape[0]

    # Okay, this is kind of ugly.  The idea is that we want to make a new basis by
    # making linear combinations of the old basis vectors, while attempting to 
    # minimize the L1 norm of the new basis vectors.  So we loop over each basis
    # vector and try to make a new one of all the vectors AFTER it in the list.
    # After this minimization is complete, we project that (now sparse) vector
    # out of the rest of them in a standard Gram-Schmidt way, and move on to the
    # next vector.  After all are done, we return the new basis.  

    #lambda function for computing L1 norm of a vector
    l1 = lambda x : np.sum( np.abs (x) )
 
    # Add a linear combination of all but the first column of B into
    # the first column of B, according to x
    combine = lambda B, x: np.dot( B, np.append( np.array([1.0,]), x) )
    
    #Loop over all the columns
    for i in range( n_cols ):

        #Find a linear combination of all the columns after the ith one that minimizes
        # the L1 norm of the ith column
        sp = opt.fmin( lambda x : l1(combine( new_basis[:, i:n_cols], x )), np.ones(n_cols-i-1), disp=0, xtol = eps)
        new_basis[:,i] = np.reshape(combine( new_basis[:, i:n_cols], sp), (n_rows,))
        new_basis[:,i] = new_basis[:,i]/linalg.norm(new_basis[:,i])

        #Now project that column out of all the others.
        for j in range (i+1, n_cols):
            new_basis[:,j] = new_basis[:,j] - np.dot(new_basis[:,i], new_basis[:,j])*new_basis[:,i]
            new_basis[:,j] = new_basis[:,j]/linalg.norm(new_basis[:,j])


    #Finally, there should now be a lot of near-zero entries in the new basis.
    #Explicitly zero them out.
    new_basis[ np.abs(new_basis) < eps ] = 0.
    return new_basis
 
def endmember_fractions_to_cvector(fvector, endmembers_per_phase):
    """
    Converts a list of molar fractions of endmembers into the
    solution variables (mineral fractions and compositions). 
    
    Parameters
    ----------
    fvector: list
        A list corresponding to molar fractions of all the endmembers
        in each of the phases.

    endmembers_per_phase : list
        A list containing the number of endmembers for each phase
        (equal to 1 for a pure phase)
   
    Returns
    -------
    composition : list
        A list of phase fractions and their compositions,
        ordered as phase fraction, followed by the molar fractions
        of the endmembers after the first endmember, i.e.
        [f(phase[i]), f(phase[i].mbr[1]), f(phase[i].mbr[2]), ...]
    """
    p=0
    composition = []
    for i, n_endmembers in enumerate(endmembers_per_phase):
        amount_phase = np.sum(fvector[p:p+n_endmembers])
        composition.append(amount_phase)
        composition.extend(fvector[p+1:p+n_endmembers]/amount_phase)
        p += n_endmembers

    return composition

def cvector_to_mineral_mbr_fractions(cvector, endmembers_per_phase):
    p=0
    composition = [[], []]
    for i, n_endmembers in enumerate(endmembers_per_phase):
        amount_phase = cvector[p]
        if amount_phase < 0.:
            amount_phase = 0.
        composition[0].append(amount_phase)
        molar_fractions = [1. - sum(cvector[p+1:p+n_endmembers])]
        molar_fractions.extend(cvector[p+1:p+n_endmembers])
        composition[1].append(molar_fractions)
        p += n_endmembers
            
    return composition
    
def mineral_mbr_fractions_to_endmember_fractions(cvectors):
    mvector = []
    for i, c in enumerate(cvectors[1]):
        mvector.extend([f*cvectors[0][i] for f in c])
    return np.array(mvector)

    
def set_eqns(PTX, assemblage, endmembers_per_phase, guessed_composition, col, null, constraints):
    """
    Set up the equations we need to solve a general equilibrium gibbs problem
    """
    
    P = PTX[0]
    T = PTX[1]
    X = PTX[2:]

    # Here are the two compositional (or P or T) constraints
    eqns = []
    for constraint in constraints:
        if constraint[0] == 'P':
            eqns.append(P - constraint[1])
        elif constraint[0] == 'T':
            eqns.append(T - constraint[1])
        elif constraint[0] == 'X':
            eqns.append(constraint[3] - np.dot(X, constraint[1])/np.dot(X, constraint[2]))
        else:
            print('Constraint type not recognised')
            exit()

    # Here we convert our guesses from a single vector
    # to a vector of phase_fractions and molar_fractions
    c = cvector_to_mineral_mbr_fractions(X, endmembers_per_phase)
    c_guess = cvector_to_mineral_mbr_fractions(guessed_composition, endmembers_per_phase)

    assemblage.set_fractions(c[0])
    for i, composition in enumerate(c[1]):
        if len(composition) > 1:
            assemblage.phases[i].set_composition(composition)

    assemblage.set_state(P, T)
    partial_gibbs = []
    for (phase, fraction) in zip(*assemblage.unroll()):
        if isinstance(phase, SolidSolution):
            partial_gibbs.extend(phase.partial_gibbs)
        else:
            partial_gibbs.append(phase.gibbs)

    # The equilibrium relation is 0 = sum(G + RT ln a)
    eqns.extend(np.dot(partial_gibbs, null))

    # On top of this, we should make sure that the bulk composition is correct
    # (i.e., that the reaction vector is in the reaction nullspace)
    guessed_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c_guess)
    new_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c)
    eqns.extend(np.dot((guessed_endmember_fractions - new_endmember_fractions), col))
    return eqns

def gibbs_minimizer(composition, assemblage, constraints, guesses=None):
    """
    Sets up a gibbs minimization problem at constant composition and attempts to solve it
    """
    
    # The next two lines set up a matrix of the endmember compositions
    # and then finds the nullspace, which corresponds to a set of independent reactions
    stoichiometric_matrix, elements, formulae, endmembers_per_phase = assemble_stoichiometric_matrix ( assemblage.phases )
    col, null = compute_column_and_null_spaces(stoichiometric_matrix)
    null = sparsify_basis(null)
    col = sparsify_basis(col)
    
    # Create a compositional vector with elements in the same order as the 
    # stoichiometric matrix
    comp_vector=np.empty( (len(elements)) )
    for idx, element in enumerate(elements):
        comp_vector[idx]=composition[element]
        
    # Check that the bulk composition can be described by a linear set of the
    # endmembers
    potential_endmember_amounts,resid,rank,s = linalg.lstsq(stoichiometric_matrix,comp_vector)
    
    resid = np.dot(stoichiometric_matrix,potential_endmember_amounts) - comp_vector
    
    ctol=1.e-3 # compositional tolerance. Users should express bulk compositions to about 3 d.p. to make sure that compositions can be recast
    if (any(abs(resid[i]) > ctol for i in range(len(resid)))):
        raise Exception("Bulk composition is well outside the compositional range of the chosen minerals (Maximum residual: "+str(max(resid))+"). Exiting.")
    else:
        # Recast composition
        comp_vector=np.dot(stoichiometric_matrix,potential_endmember_amounts)


    initial_composition = endmember_fractions_to_cvector(potential_endmember_amounts, endmembers_per_phase)


    # If an initial guess hasn't been given, make one
    if guesses == None:
        guesses = [10.e9, 1200.]
        guesses.extend(initial_composition)
        
        for constraint in constraints:
            if constraint[0] == 'P':
                guesses[0] = constraint[1]
            if constraint[1] == 'T':
                guesses[1] = constraint[1]

    
    # Set up the problem and attempt to solve it
    # Ignore warning due to changing the phase fractions such that they don't equal 1. They're renormalised anyway.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        soln=opt.fsolve(set_eqns, guesses,
                        args=(assemblage, endmembers_per_phase, initial_composition, col, null, constraints),
                        xtol=1e-10, full_output=True)
    if soln[2]==1:
        return soln[0]
    else:
        return soln[3]
        sys.exit()


def binary_composition(x, composition1, composition2):
    composition = {}
    for key, value in composition1.items():
        composition[key] = value*(1. - x)
    for key, value in composition2.items():
        if key in composition:
            composition[key] += value*x
        else:
            composition[key] = value*x

    return composition

        
def gibbs_bulk_minimizer(composition1, composition2, guessed_bulk, assemblage, constraints, guesses=None):

    def minimize_P(x, composition1, composition2, assemblage, constraints, fixed_P, guesses):
        composition = binary_composition(x, composition1, composition2)
        return fixed_P - gibbs_minimizer(composition, assemblage, constraints)[0]

    
    def minimize_T(x, composition1, composition2, assemblage, constraints, fixed_T, guesses):
        composition = binary_composition(x, composition1, composition2)
        return fixed_T - gibbs_minimizer(composition, assemblage, constraints)[1]


    c = 0
    new_constraints = []
    for constraint in constraints:
        if (constraint[0] == 'P' or constraint[0] == 'T') and c==0:
            master_constraint = constraint
            c=1
        else:
            new_constraints.append(constraint)

    if master_constraint[0] == 'P':
        func = minimize_P
    else:
        func = minimize_T
        
    soln = opt.fsolve(func, [guessed_bulk], args=(composition1, composition2,
                                                  assemblage, new_constraints, master_constraint[1], guesses),
                      full_output=True, xtol=1.e-10)
    if soln[2]==1:
        return soln[0]
    else:
        return soln[3]
        sys.exit()
