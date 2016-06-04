# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import scipy.linalg as linalg
import scipy.optimize as opt
import warnings

from .mineral import Mineral
from .solidsolution import SolidSolution
from .composite import Composite
from .processchemistry import process_solution_chemistry, dependent_endmembers


def _get_formulae_indices_endmembers(composition, minerals):
    """
    Takes a bulk composition and a endmembers and a list of minerals.
    Creates a list of endmember formulae, a unique index according to the
    (mineral, endmember) pair, and a list of the number of endmembers per phase.
    Only includes endmembers in the lists if all of their elements are
    contained in the given composition.

    Parameters
    ----------
    composition : dictionary of floats
        Compositional dictionary, given in number of atoms per element. 
        Used to produce a consistent number and ordering of elements and 
        remove endmembers which lie outside the given composition.

    minerals : List of burnman.Mineral 
        List of minerals

    Returns
    -------
    formulae : list of dictionaries
        Formulae for each of the endmembers in the minerals

    indices : list of list of integers
        List of integer pairs corresponding to the indices of the 
        mineral in the list of minerals and endmember in that mineral. 

    endmembers_per_phase : list of integers
        Number of endmembers per mineral. Pure phases have 1 endmember.
    """
    
    #Listify the elements and sort them so they have a consistent order.
    #This will be the ordering of the rows.  The ordering of the columns
    #will be the ordering of the endmembers as they are passed in.
    elements = list(set(composition.keys()))

    # Create the desired lists
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
    """
    Takes a bulk composition and stoichiometric matrix of the endmembers, and 
    attempts to find phase proportions of the endmembers which produce the bulk 
    composition. This will not in general be a unique solution; additional 
    constraints will be required (energy minimization, partitioning data, for example).
    
    Parameters
    ----------
    comp_vector : array of floats
        Amount of each element in the bulk composition

    stoichiometric_matrix : 2d array of floats
        Number of atoms of each element for each endmember

    new_constraints : lists of floats
        Lists of bulk compositional constraints.
        For example, if the amount of particular phase(s) must be a given atom proportion, 
        of the total, the indices (in the stoichiometric matrix) and proportions 
        of that phase are provided as a constraint.

    Returns
    -------
    potential_endmember_amounts : array of floats
        numpy array providing one solution to the amounts of each endmember 
        which satisfy the bulk composition and additional compositional constraints

    check_composition : boolean
        A check that the residual between the bulk composition implied by 
        the endmember amounts and that provided by the input are equal 
        within a tolerance given in the function.
    """
    
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
    norm = 1./sum(guessed_amounts)
    cons = [{'type': 'ineq', 'fun': lambda x, i=[i]: x[i[0]]*norm} for i in xrange(len(stoichiometric_matrix[0]))]
    sol = opt.minimize(objective, guessed_amounts, args=(b*norm, A), method='SLSQP', constraints=cons, tol=1.e-15, options={'eps': 1.e-15, 'maxiter' :1000})
    check = sol.fun<1.e-15
    return sol.x/norm, check
    
def make_dpdnt_stoichiometric_matrix(stoichiometric_matrix,
                                     indices, minerals,
                                     endmembers_per_phase):

    dependent_stoichiometric_matrix = np.array([]).reshape(0, len(stoichiometric_matrix))
    dependency_matrix = np.identity(len(indices))
    
    n_start = 0
    for i, n_mbrs in enumerate(endmembers_per_phase):
        mbr_indices = [(index, mbr_idx) for index, (phase_idx, mbr_idx) in enumerate(indices) if phase_idx == i]
        if len(mbr_indices) > 1:
            dependent_vectors = dependent_endmembers([minerals[i].solution_model.formulas[idx[1]] for idx in mbr_indices])
            if dependent_vectors != []:
                vector_to_matrix = np.zeros((len(dependent_vectors), len(indices)))
                for idx, (i, j) in enumerate(mbr_indices):
                    for v_idx, v in enumerate(dependent_vectors):
                        vector_to_matrix[v_idx][i] = v[idx]

                dependency_matrix = np.concatenate((dependency_matrix, vector_to_matrix))
                
                independent_compositions = np.array([stoichiometric_matrix.T[idx[0]] for idx in mbr_indices])
                matrix=np.dot(independent_compositions.T, dependent_vectors.T)
                dependent_stoichiometric_matrix = np.concatenate((dependent_stoichiometric_matrix, matrix.T))
        n_start = n_start + n_mbrs

    dependent_stoichiometric_matrix = np.array(dependent_stoichiometric_matrix).T
    return dependent_stoichiometric_matrix, dependency_matrix
    
def assemble_compositional_tensors ( composition, minerals, constraints ):
    """
    This function takes a list of minerals and assembles a matrix where 
    the rows are elements and the columns are species (or endmembers).
    If a solid solution is passed in, then the endmembers are extracted
    from it. 
    
    Parameters
    ----------
    composition : dict of floats

    minerals : list of minerals
        List of objects of type :class:`burnman.Mineral` or `burnman.SolidSolution`.
        Other types cannot be understood, and an exception will be thrown.

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
    formulae, indices, endmembers_per_phase = _get_formulae_indices_endmembers(composition, minerals)
    elements = list(set(composition.keys()))
    
    #Populate the stoichiometric matrix
    stoichiometric_matrix = np.array( [[ ( f[e]  if e in f else 0.0 ) for f in formulae ]
                                       for e in elements])

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
        # Recast composition. We already checked that the composition is ok,
        # so this should have a negligible effect on the bulk composition.
        comp_vector = np.dot(stoichiometric_matrix,potential_endmember_amounts)
        
    # Try to remove endmembers which cannot be constituents of the solution
    # Do this by checking that the bulk composition can still be satisfied
    # if the amount of the endmember is some small value
    delta = 1.e-5
    null_endmember_indices = [ i for i, amount in enumerate(potential_endmember_amounts) if amount < 1.e-10]
    for i in reversed(null_endmember_indices):
        potential_endmember_amounts, composition_satisfied = potential_amounts(comp_vector,
                                                                               stoichiometric_matrix,
                                                                               [[[i], delta]])
        if not composition_satisfied:
            stoichiometric_matrix = np.delete(stoichiometric_matrix, i, 1)
            del indices[i]
            del formulae[i]

    # Now, let's find the dependent endmembers for each of the solid solutions:
    dependent_stoichiometric_matrix, dependency_matrix = make_dpdnt_stoichiometric_matrix(stoichiometric_matrix, indices, minerals, endmembers_per_phase)

    complete_stoichiometric_matrix = np.concatenate((stoichiometric_matrix,
                                                     dependent_stoichiometric_matrix), axis=1)
    
    # Now we can find a potential solution to our problem (one that satisfies the bulk composition constraints)
    new_constraints = []
    for compositional_constraint in compositional_constraints:
        new_constraints.append([[i for i, idx in enumerate(indices) if idx[0] == compositional_constraint[0] ], compositional_constraint[1]])
    if new_constraints != []:
        potential_endmember_amounts, composition_satisfied = potential_amounts(comp_vector, complete_stoichiometric_matrix, new_constraints)
    else:
        potential_endmember_amounts, composition_satisfied = potential_amounts(comp_vector, complete_stoichiometric_matrix)

    if not composition_satisfied:
        raise Exception('Composition not satisfied. This is a bug; please contact the BurnMan team.')

    # Now let's rework the dependent endmembers back into the independent set
    potential_endmember_amounts = np.dot(dependency_matrix.T, potential_endmember_amounts)
    
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
        if i_idx == old_idx and amount < 1.e-3:
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
        The stoichiometric matrix
    
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
    # Do a SVD of the stoichiometric matrix.
    U, S, Vh = linalg.svd( stoichiometric_matrix)

    # The columns of V that correspond to large and small (or nonexistent)
    # singular values are the column and left null spaces for the matrix.
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

    indices: list of lists of integers

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

    mvector = [[0. for j in xrange(endmembers_per_phase[i])] for i in xrange(len(endmembers_per_phase))]
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


