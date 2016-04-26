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

def assemble_stoichiometric_matrix ( minerals):
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
 
def endmember_fractions_to_cvector(mvector, endmembers_per_phase):
    """
    Converts a list of molar fractions of endmembers into the
    solution variables (mineral fractions and compositions). 
    
    Parameters
    ----------
    mvector: list
        A list corresponding to molar fractions of all the endmembers
        in each of the phases.

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
    composition = []
    for i, n_endmembers in enumerate(endmembers_per_phase):
        amount_phase = np.sum(mvector[p:p+n_endmembers])
        composition.append(amount_phase)
        composition.extend(mvector[p+1:p+n_endmembers]/amount_phase)
        p += n_endmembers

    return composition

def cvector_to_mineral_mbr_fractions(cvector, endmembers_per_phase):
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
    for i, c in enumerate(cvectors[1]):
        mvector.extend([f*cvectors[0][i] for f in c])
    return np.array(mvector)

    
def set_eqns(PTX, assemblage, endmembers_per_phase, initial_composition, col, null, constraints):
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
    c = cvector_to_mineral_mbr_fractions(X, endmembers_per_phase)
    c_initial = cvector_to_mineral_mbr_fractions(initial_composition, endmembers_per_phase)

    assemblage.set_fractions(c[0])
    for i, composition in enumerate(c[1]):
        if len(composition) > 1:
            assemblage.phases[i].set_composition(composition)

    assemblage.set_state(P, T)
    partial_gibbs = []
    for (phase, fraction) in zip(*assemblage.unroll()):
        if isinstance(phase, SolidSolution):
            partial_gibbs.extend(phase.partial_gibbs)
            print(phase.activities)
        else:
            partial_gibbs.append(phase.gibbs)

    # The equilibrium relation is 0 = sum(G + RT ln a)
    eqns.extend(np.dot(partial_gibbs, null))
    print('p', partial_gibbs)
    print(eqns)
    print(null)
    exit()
    # On top of this, we should make sure that the bulk composition is correct
    # (i.e., that the reaction vector is in the reaction nullspace)
    initial_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c_initial)
    new_endmember_fractions = mineral_mbr_fractions_to_endmember_fractions(c)
    eqns.extend(np.dot((initial_endmember_fractions - new_endmember_fractions), col))

    return eqns

def compositional_variables(assemblage):
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
    stoichiometric_matrix, elements, formulae, endmembers_per_phase = assemble_stoichiometric_matrix ( assemblage.phases )
    var_names=[]
    for i, phase in enumerate(assemblage.phases):
        var_names.append('x('+phase.name+')')
        for j in xrange(endmembers_per_phase[i] - 1):
            var_names.append('p('+phase.endmembers[j+1][0].name+')')
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
        sol_dict = {'P': soln[0][0], 'T': soln[0][1]}
        for i, variable in enumerate(compositional_variables(assemblage)):
            sol_dict[variable] = soln[0][2+i]
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
        composition = binary_composition(composition1, composition2, x)
        vecsol = gibbs_minimizer(composition, assemblage, constraints)
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

    stoichiometric_matrix, elements, formulae, endmembers_per_phase = assemble_stoichiometric_matrix ( assemblage.phases )
    c0a = []
    c0b = []
    c1 = []
    for n_endmembers in endmembers_per_phase:
        c0a.append(0.)
        c0b.append(0.)
        c1.append(1.)
        for i in xrange(n_endmembers - 1):
            c0.append(0.)
            c1.append(0.)
    c0a[0] = 1.
    c0b[endmembers_per_phase[0]] = 1.

    constraints= [['X', c0a, c1, 0.], ['X', c0b, c1, 0.]]
    soln_array = gibbs_minimizer(composition, assemblage, constraints, guesses)

    return (soln_array['P'], soln_array['T'])    
    
def find_univariant(composition, phases, zero_phase, condition_variable, condition_array, guesses=None):
    all_phases = [zero_phase]
    all_phases.extend(phases)
    assemblage = Composite(all_phases)

    
    stoichiometric_matrix, elements, formulae, endmembers_per_phase = assemble_stoichiometric_matrix ( assemblage.phases )
    c0 = []
    c1 = []
    for n_endmembers in endmembers_per_phase:
        c0.append(0.)
        c1.append(1.)
        for i in xrange(n_endmembers - 1):
            c0.append(0.)
            c1.append(0.)
    c0[0] = 1.
    
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
