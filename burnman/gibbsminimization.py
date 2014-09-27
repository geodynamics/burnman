# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.linalg as linalg
import scipy.optimize as opt
import burnman
from burnman.processchemistry import *


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

    # Make a list of the different formulae, as well as a set of 
    # the elements that we have present in the mineral list
    for m in minerals:
        # Add the endmembers if it is a solid solution
        if isinstance(m, burnman.SolidSolution):
            for e in m.base_material:
                f = e[0].params['formula']
                formulae.append(f)
                for k in f.keys():
                    elements.add(k)
        # Add formula if it is a simple mineral
        elif isinstance(m, burnman.Mineral):
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

    return stoichiometric_matrix, elements, formulae


def compute_nullspace ( stoichiometric_matrix ):
    """
    Given a stoichiometric matrix, compute a basis for the nullspace.
    This basis corresponds to the subspace that does not change the 
    bulk composition.  Therefore the vectors correspond to chemical
    reactions.  This merely identifies a set of linearly independent
    reactions, without saying anything about which ones are likely to
    occur.

    The calculation of the nullspace is done with and SVD.
    
    Parameters
    ----------
    stoichiometric_matrix: numpy array
        The stoichiometric matrix, presumably calculated using :func:`compute_stoichiometric_matrix`.
    
    Returns
    -------
    nullspace: numpy array
        A 2D array correspnding to the nullspace of the stoichiometric matrix, 
        with the columns corresponding to basis vectors.
    """

    eps = 1.e-10
    # Do an SVD of the stoichiometric matrix.
    U, S, Vh = linalg.svd( stoichiometric_matrix)

    # The columns of V that correspond to small singular values
    # (or do not exist) are the left nullspace for the matrix.
    # select them, appropriately transpose them, and return them
    S=np.append(S, [0. for i in range(len(Vh)-len(S))])
    null_mask = (S <= eps)
    null_space = np.compress(null_mask, Vh, axis=0)

    return np.transpose(null_space)
    

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
 
    
