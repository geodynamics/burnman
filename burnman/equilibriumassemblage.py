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
from .composite import Composite
import warnings

def get_formulae_indices_endmembers(composition, minerals):
    #Listify the elements and sort them so they have a consistent order.
    #This will be the ordering of the rows.  The ordering of the columns
    #will be the ordering of the endmembers as they are passed in.
    elements = list(set(composition.keys()))
    
    formulae = []
    indices = []
    endmembers_per_phase = []
    for m_idx, m in enumerate(minerals):
        # Add the endmembers if it is a solid solution
        if isinstance(m, SolidSolution):            
            endmembers_per_phase.append(len(m.endmembers))
            for e_idx, e in enumerate(m.endmembers):
                f = e[0].params['formula']
                if all(keys in elements for keys in f.keys()):
                    formulae.append(f)
                    indices.append([m_idx, e_idx])
        # Add formula if it is a simple mineral
        elif isinstance(m, Mineral):
            endmembers_per_phase.append(1)
            f = m.params['formula']
            if all(keys in elements for keys in f.keys()):
                formulae.append(f)
                indices.append([m_idx, 0])
        else:
            raise Exception('Unsupported mineral type, can only read burnman.Mineral or burnman.SolidSolution')
    return formulae, indices, endmembers_per_phase

def potential_amounts(comp_vector, stoichiometric_matrix, new_constraints=None):
    if new_constraints is not None:
        total_atoms = sum(comp_vector)
        A_additions = np.array([[ 0. for f in stoichiometric_matrix[0] ] for new_constraint in new_constraints])
        for i, new_constraint in enumerate(new_constraints):
            for index in new_constraint[0]:
                A_additions[i][index] = sum(stoichiometric_matrix.T[i])
        A = np.concatenate((stoichiometric_matrix, A_additions), axis=0)
        b = np.concatenate((comp_vector, [proportion*total_atoms for indices, proportion in new_constraints]))
    else:
        A = stoichiometric_matrix
        b = comp_vector

    def xeval(x, A):
        return np.dot(A, x)

    def objective(x, b, A):
        return np.sum((b - xeval(x, A)) ** 2)

    guessed_amounts,resid,rank,s = linalg.lstsq(A, b)

    # Apply inequality constraints to make sure that the phase amounts are positive
    # At the moment, this doesn't take into account the possibility of negative endmember
    # amounts in solid solutions (to create dependent endmembers; we need to fix this)
    cons = [{'type': 'ineq', 'fun': lambda x, i=[i]: x[i[0]]} for i in xrange(len(stoichiometric_matrix[0]))]
    
    sol = opt.minimize(objective, guessed_amounts, args=(b, A), method='SLSQP', constraints=cons)
    potential_endmember_amounts = sol.x
    resid = np.dot(A,potential_endmember_amounts) - b
    
    ctol=1.e-12 # compositional tolerance.
    if (any(abs(residual) > ctol for residual in resid)):
        return potential_endmember_amounts, False
    else:
        return potential_endmember_amounts, True

def assemble_compositional_tensors ( composition, minerals, constraints ):
    """
    This function takes a list of minerals and assembles a matrix where 
    the rows are elements and the columns are species (or endmembers).
    If a solid solution is passed in, then the endmembers are extracted
    from it. 
    
    Parameters
    ----------
    minerals : list of minerals
        List of objects of type :class:`burnman.Mineral` or `burnman.SolidSolution`.
        Other types cannot be understood, and an exception will be thrown.
    composition : dict of floats
    
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
    indices: list of list of integers
    
    """
    
    formulae, indices, endmembers_per_phase = get_formulae_indices_endmembers(composition, minerals)
    elements = list(set(composition.keys()))
    
    #Populate the stoichiometric matrix
    stoichiometric_matrix = np.empty( [ len(elements), len(formulae) ] )
    for i,e in enumerate(elements):
        for j,f in enumerate(formulae):
            stoichiometric_matrix[i,j] = ( f[e]  if e in f else 0.0 )
    
    # Check that the bulk composition can be described by a linear set of the endmembers
    comp_vector = np.array([composition[element] for element in elements])

    compositional_constraints = []
    for constraint in constraints:
        if constraint[0] == 'X':
            c = cvector_to_mineral_mbr_fractions(constraint[1], indices, endmembers_per_phase)
            c1 = cvector_to_mineral_mbr_fractions(constraint[2], indices, endmembers_per_phase)
            if (any(c[0]) == 1. and sum(c[0]) == 1.) and all(c1[0]) == 1.:
                compositional_constraints.append([c[0].index(1.), constraint[3]])
            
    potential_endmember_amounts, check_composition = potential_amounts(comp_vector, stoichiometric_matrix)
    if check_composition == False:
        raise Exception("Bulk composition is well outside the compositional range of the chosen minerals. Exiting.")
    else:
        # Recast composition. We already checked that the composition is ok, so this should have a negligible effect on the bulk composition.
        comp_vector = np.dot(stoichiometric_matrix,potential_endmember_amounts)
        
    # Try to remove endmembers which cannot be constituents of the solution
    null_endmember_indices = [ i for i, amount in enumerate(potential_endmember_amounts) if amount < 1.e-10]
    for i in reversed(null_endmember_indices):
        # Add constraint to stoichiometric matrix that the amount of the phase is very small, check that constraints can still be satisfied.
        # If not, remove endmember
        # We do this by adding another "element" to the stoichiometric matrix, with zeros everywhere except for a 1. for the endmember of interest
        # The compositional vector gets a delta for this last element
        delta = 1.e-5
        potential_endmember_amounts, check_composition = potential_amounts(comp_vector, stoichiometric_matrix, [[[i], delta]])
        if check_composition == False:
            #print('Removing endmember', i, 'with formula: ', formulae[i])
            stoichiometric_matrix = np.delete(stoichiometric_matrix, i, 1)
            del indices[i]
            del formulae[i]

    
    # Now we can find a potential solution to our problem (one that satisfies the bulk composition constraints)
    new_constraints = []
    for compositional_constraint in compositional_constraints:
        new_constraints.append([[i for i, idx in enumerate(indices) if idx[0] == compositional_constraint[0] ], compositional_constraint[1]])
    if new_constraints != []:
        potential_endmember_amounts, check_composition = potential_amounts(comp_vector, stoichiometric_matrix, new_constraints)
    else:
        potential_endmember_amounts, check_composition = potential_amounts(comp_vector, stoichiometric_matrix)
    
    col, null = compute_column_and_null_spaces(stoichiometric_matrix)
    null = sparsify_basis(null)
    col = sparsify_basis(col)

    initial_composition = endmember_fractions_to_cvector(potential_endmember_amounts, indices, endmembers_per_phase)
    # Finally, we don't want the solutions with modal abundances = zero
    # to have endmember compositions. We fairly arbitrarily choose small amounts of the secondary endmembers.
    old_idx = -1
    amount = 1.
    for i, index in enumerate(indices):
        i_idx, j_idx = index
        if i_idx == old_idx and amount < 1.e-5:
            initial_composition[i] = 1./endmembers_per_phase[i_idx]/2.
        else:
            old_idx = i_idx
            amount = initial_composition[i]

    return col, null, initial_composition, indices, endmembers_per_phase


def compute_column_and_null_spaces ( stoichiometric_matrix ):
    """
    Given a stoichiometric matrix, compute a basis for the column and nullspaces.
    These bases corresponds to the subspaces that do and do not change the 
    bulk composition.  The vectors of the nullspace correspond to chemical
    reactions, which merely identify a set of linearly independent
    reactions, without saying anything about which ones are likely to
    occur.

    The calculation of the column and nullspaces is done with SVD.
    
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
 
def endmember_fractions_to_cvector(partial_mvector, indices, endmembers_per_phase):
    """
    Converts a list of molar fractions of endmembers along with their
    indices into the solution variables
    (mineral fractions and compositions). 
    
    Parameters
    ----------
    mvector: list
        A list corresponding to molar fractions of all the endmembers
        in each of the phases.

    indices: list of list of integers

    endmembers_per_phase : list
        A list containing the number of endmembers for each phase
        (equal to 1 for a pure phase).
   
    Returns
    -------
    composition : list
        A list of phase fractions and their compositions,
        ordered as phase fraction, followed by the molar fractions
        of the endmembers after the first endmember, i.e.
        [f(phase[i]), f(phase[i].mbr[1]), f(phase[i].mbr[2]), ...]
    """
    p=0

    mvector = [[0. for i in xrange(endmembers_per_phase[i])] for i in xrange(len(endmembers_per_phase))]
    for i_idx, i in enumerate(indices):
        mvector[i[0]][i[1]] = partial_mvector[i_idx]

    # Need to return phase quantities and compositions in the reduced endmember space
    composition = []
    for i, n_endmembers in enumerate(endmembers_per_phase):
        amount_phase = np.sum(mvector[i])
        composition.append(amount_phase)

        endmember_indices = [j_idx for (i_idx, j_idx) in indices if i == i_idx]
        for j in xrange(len(endmember_indices)-1):
            composition.append(mvector[i][endmember_indices[j+1]]/amount_phase)

    return composition

def cvector_to_mineral_mbr_fractions(cvector, indices, endmembers_per_phase):
    """
    Converts a compositional list such as that computed by
    :func:`endmember_fractions_to_cvector` into two lists, one
    containing the amounts of the different phases
    (endmembers or solutions) and the second containing the molar
    fractions of the endmembers in each solid solution. If the
    phase is an endmember, the molar fraction is equal to 1.

    Parameters
    ----------
    cvector: list
        A list corresponding to the amounts of each phase and
        the molar fractions of all the endmembers (apart from the first)
        in each of the phases.

    endmembers_per_phase : list
        A list containing the number of endmembers for each phase
        (equal to 1 for a pure phase).

    Returns
    -------
    composition : list of lists
        Two lists, one containing phase amounts, and the other
        containing molar fractions of the endmembers in each phase.
    """
    
    phase_amounts = [0. for n_endmembers in endmembers_per_phase]
    molar_fractions = [[0. for i_mbr in xrange(n_endmembers)] for n_endmembers in endmembers_per_phase]

    old_idx = -1
    first_j_idx = 0
    
    for i, c in enumerate(cvector):
        i_idx, j_idx = indices[i]
        if i_idx != old_idx:
            if c >= 0.:
                phase_amounts[i_idx] = c
            else:
                phase_amounts[i_idx] = 0.
                
            molar_fractions[i_idx][j_idx] = 1.
            old_idx = i_idx
            first_j_idx = j_idx
        else:
            molar_fractions[i_idx][j_idx] = c
            molar_fractions[i_idx][first_j_idx] -= c

    
    return [phase_amounts, molar_fractions]
    
def mineral_mbr_fractions_to_endmember_fractions(cvectors, indices):
    """
    Converts a list of lists containing phase amounts and
    the molar fractions of the endmembers in each phase into
    a numpy array of endmember fractions.

    Parameters
    ----------
    composition : list of lists
        Two lists, one containing phase amounts, and the other
        containing molar fractions of the endmembers in each phase.

    Returns
    -------
    mvector: list
        A list corresponding to molar fractions of all the endmembers
        in each of the phases.

    """
    mvector = []
    for (i, j) in indices:
        mvector.append(cvectors[0][i]*cvectors[1][i][j])
        
    return np.array(mvector)

    
def set_eqns(PTX, assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints):
    """
    Set up the equations we need to solve a general equilibrium gibbs problem.

    Parameters
    ----------
    PTX : list of floats
        A list of pressure, temperature and compositional variables.
        The compositional part of the list described the
        composition of the assemblage, ordered as
        phase fraction, followed by the molar fractions
        of the endmembers after the first endmember, i.e.
        [f(phase[i]), f(phase[i].mbr[1]), f(phase[i].mbr[2]), ...].
        
    assemblage : composite
        The assemblage of minerals for which we want to find
        amounts and compositions.
        
    endmembers_per_phase : list of floats
        A list of length len(composite.phases) containing
        the number of endmembers in each phase.
        
    initial_composition : list of floats
        A list with the same structure as the X part of PTX,
        describing any composition which satisfies the bulk
        composition to be investigated. Some of the returned
        eqns are only zero if the bulk composition of X
        and initial composition are the same. 
        
    col : numpy array
        The column space of the stoichiometric matrix.
        
    null : numpy array
        The null space of the stoichiometric matrix.

    indices :
        The indices of the endmembers corresponding to the columns
        of the stoichiometric matrix
    
    constraints : list of lists
        Two pressure, temperature or compositional constraints
        on the problem of interest. The pressure
        (or temperature) constraints have the form
        [['P'], P_value].
        The compositional constraints have the form
        [['X'], list l0, list l1, float f]
        where list[i] is of len(initial_composition),
        and the constraint is
        dot(l0, X)/dot(l1, X) = f.
    
    Returns
    -------
    eqns: list of floats
        A list where all the elements are zero iff the
        variables in PTX satisfy the equilibrium relation and
        the bulk composition constraints.
    """
    
    P = PTX[0]
    T = PTX[1]
    X = PTX[2:]

    # Here are the two P, T or compositional constraints
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
    c = cvector_to_mineral_mbr_fractions(X, indices, endmembers_per_phase)
    c_initial = cvector_to_mineral_mbr_fractions(initial_composition, indices, endmembers_per_phase)

    if any(f < 0. for f in c_initial[0]):
        raise Exception('The starting guess contains a phase amount < 0. If you did not provide a starting guess, this is a bug. Please contact the burnman team.')

    if any(f < 0. for f in c[0]):
        raise Exception('Gibbs minimization failed as a result of a phase amount being < 0. This might indicate that the starting guess is poor, or that there is no solution to this problem')

    assemblage.set_fractions(np.array(c[0])/sum(c[0]))
    for i, composition in enumerate(c[1]):
        if len(composition) > 1:
            assemblage.phases[i].set_composition(composition)

    assemblage.set_state(P, T)
    full_partial_gibbs = []
    for (phase, fraction) in zip(*assemblage.unroll()):
        if isinstance(phase, SolidSolution):
            full_partial_gibbs.append(phase.partial_gibbs)
        else:
            full_partial_gibbs.append([phase.gibbs])

    # Remove endmembers from the partial_gibbs vector if they do not fall within the
    # compositional space given by "indices"
    partial_gibbs = []
    for (i_idx, j_idx) in indices:
        partial_gibbs.append(full_partial_gibbs[i_idx][j_idx])
            
    # The equilibrium relation is 0 = sum(G + RT ln a)
    eqns.extend(np.dot(partial_gibbs, null))
    
    # On top of this, we should make sure that the bulk composition is correct
    # (i.e., that the reaction vector is in the reaction nullspace)
    initial_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c_initial, indices)
    new_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c, indices)
    eqns.extend(np.dot((initial_endmember_fractions - new_endmember_fractions), col))

    return eqns

def compositional_variables(assemblage, indices):
    """
    Takes an assemblage and outputs names for the
    compositional variables which describe the bulk composition
    and compositions of all the phases.

    Parameters
    ----------
    assemblage : composite
        The assemblage of minerals for which we want to find
        amounts and compositions.
    
    Returns
    -------
    var_names : list of strings
        Strings are provided in the same order as the X part
        of the PTX variable input to :func:`set_eqns`. Phase
        amount names are given as 'x(phase.name)', while
        the molar fractions of endmembers in each phase are
        given as 'p(endmember.name)'
    """
    old_i = -1
    var_names=[]
    for (i, j) in indices:
        if i != old_i:
            var_names.append('x('+assemblage.phases[i].name+')')
            old_i = i
        else:
            var_names.append('p('+assemblage.phases[i].endmembers[j][0].name+')')

    return var_names
    
def gibbs_minimizer(composition, assemblage, constraints, guesses=None):
    """
    Sets up a gibbs minimization problem at constant composition and attempts to solve it

    Parameters
    ----------
    composition : dictionary of floats
        Dictionary contains the number of atoms of each element.
        
    assemblage : composite
        The assemblage of minerals for which we want to find
        amounts and compositions.

    constraints : list of lists
        Two pressure, temperature or compositional constraints
        on the problem of interest. The pressure
        (or temperature) constraints have the form
        [['P'], P_value].
        The compositional constraints have the form
        [['X'], list l0, list l1, float f]
        where list[i] is of len(initial_composition),
        and the constraint is
        dot(l0, X)/dot(l1, X) = f.

    guesses : optional list of floats
        List has the same form as PTX in :func:`set_eqns`
        
    Returns
    -------
    sol_dict : dictionary of floats
        Dictionary contains the solution vector of the minimization.
        Dictionary keys are 'P', 'T' and the variable names output by
        :func:`compositional_variables`
    """

    # Redefine composition removing negligible elements (where n atoms < 0.0001 % total atoms)
    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}
        
    # The function assemble_compositional_tensors sets up a matrix of the endmember compositions
    # and finds the column and left nullspace which corresponds to a set of independent reactions
    col, null, initial_composition, indices, endmembers_per_phase = assemble_compositional_tensors ( composition, assemblage.phases, constraints )
    
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
    soln=opt.fsolve(set_eqns, guesses,
                    args=(assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints),
                    xtol=1e-10, full_output=True)

    
    if soln[2]==1:
        sol_dict = {'P': soln[0][0], 'T': soln[0][1], 'c': soln[0][2:2+len(indices)], 'variables': compositional_variables(assemblage, indices)}
        return sol_dict
    else:
        raise Exception('Solution could not be found, error: \n'+soln[3]+'\n Guesses:'
                        +str(guesses))


def binary_composition(composition1, composition2, x):
    """
    Returns the composition within a binary system that
    is defined by a fraction of the second composition.
    
    Parameters
    ----------
    composition1 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    composition2 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    x : float
        Composition as a fraction of composition2.

    Returns
    -------
    composition : dictionary of floats
        (1-x)*composition1 + x*composition2.
    """
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
    """
    Sets up a gibbs minimization problem with variable composition and attempts to solve it

    Parameters
    ----------
    composition1 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    composition2 : dictionary of floats
        Dictionary contains the number of atoms of each element
        at one end of a binary.

    guessed_bulk : float
        Guessed composition satisfying the constraints, as a
        fraction of composition2
        
    assemblage : composite
        The assemblage of minerals for which we want to find
        amounts and compositions.

    constraints : list of lists
        Two pressure, temperature or compositional constraints
        on the problem of interest. The pressure
        (or temperature) constraints have the form
        [['P'], P_value].
        The compositional constraints have the form
        [['X'], list l0, list l1, float f]
        where list[i] is of len(initial_composition),
        and the constraint is
        dot(l0, X)/dot(l1, X) = f.

    guesses : optional list of floats
        List has the same form as PTX in :func:`set_eqns`
        
    Returns
    -------
    sol_dict : dictionary of floats
        Dictionary contains the solution vector of the minimization.
        Dictionary keys are 'P', 'T', 'X' and the variable names
        output by :func:`compositional_variables`
    """
    
    def minimize(x, composition1, composition2, assemblage, constraints, master_constraint, guesses):
        composition = binary_composition(composition1, composition2, x[0])
        vecsol = gibbs_minimizer(composition, assemblage, constraints, guesses)
        return master_constraint[1] - vecsol[master_constraint[0]]

    c = 0
    new_constraints = []
    for constraint in constraints:
        if (constraint[0] == 'P' or constraint[0] == 'T') and c==0:
            master_constraint = constraint
            c=1
        else:
            new_constraints.append(constraint)
            
    soln = opt.fsolve(minimize, [guessed_bulk], args=(composition1, composition2,
                                                      assemblage, new_constraints,
                                                      master_constraint, guesses),
                                                      full_output=True, xtol=1.e-10)
        
    if soln[2]==1:
        # One more run to get the full solution vector
        composition = binary_composition(composition1, composition2, soln[0][0])
        sol_dict = gibbs_minimizer(composition, assemblage, new_constraints, guesses)
        sol_dict['X'] = soln[0][0]
        return sol_dict
    else:
        raise Exception('Solution could not be found, error: '+soln[3])
        

def find_invariant(composition, phases, zero_phases, guesses=None):

    all_phases = zero_phases
    all_phases.extend(phases)
    assemblage = Composite(all_phases)
    base_assemblage = Composite(phases)

    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}

    
    formulae, indices, endmembers_per_phase = get_formulae_indices_endmembers(composition, assemblage.phases)
    
    c0a = [0. for index in indices]
    c0b = [0. for index in indices]
    c1 = [1.]
    c1.extend([float(t - s) for s, t in zip(zip(*indices)[0], zip(*indices)[0][1:])])
    c0a[0] = 1.
    c0b[c1.index(1., 1)] = 1.

    constraints= [['X', c0a, c1, 0.], ['X', c0b, c1, 0.]]
    soln_array = gibbs_minimizer(composition, assemblage, constraints, guesses)

    return (soln_array['P'], soln_array['T'])    
    
def find_univariant(composition, phases, zero_phase, condition_variable, condition_array, guesses=None):
    all_phases = [zero_phase]
    all_phases.extend(phases)
    assemblage = Composite(all_phases)
    
    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}

    formulae, indices, endmembers_per_phase = get_formulae_indices_endmembers(composition, assemblage.phases)
    
    c0 = [0. for index in indices]
    c0[0] = 1.
    c1 = [1.]
    c1.extend([float(t - s) for s, t in zip(zip(*indices)[0], zip(*indices)[0][1:])])
    
    soln_array = np.empty_like(condition_array)
    if condition_variable=='P':
        for i, P in enumerate(condition_array):
            constraints= [['P', P], ['X', c0, c1, 0.]]
            soln_array[i] = gibbs_minimizer(composition, assemblage, constraints, guesses)['T']
    elif condition_variable=='T':
        for i, T in enumerate(condition_array):
            constraints= [['T', T], ['X', c0, c1, 0.]]
            soln_array[i] = gibbs_minimizer(composition, assemblage, constraints, guesses)['P']

    return soln_array
