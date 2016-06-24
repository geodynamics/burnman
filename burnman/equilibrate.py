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
from .processchemistry import process_solution_chemistry, compositional_variables, binary_composition
from .compositealgebra import _get_formulae_indices_endmembers, assemble_compositional_tensors, cvector_to_mineral_mbr_fractions, mineral_mbr_fractions_to_endmember_fractions

def _set_eqns(PTX, assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints):
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
            raise Exception('Constraint type not recognised. Equilibration calculations accept pressure (P), temperature (T) or compositional (X) constraints.')

    # Here we convert our guesses from a single vector
    # to a vector of phase_fractions and molar_fractions
    c = cvector_to_mineral_mbr_fractions(X, indices, endmembers_per_phase)
    c_initial = cvector_to_mineral_mbr_fractions(initial_composition, indices, endmembers_per_phase)

    if any(f < 0. for f in c_initial[0]):
        raise Exception('The starting guess contains a phase amount < 0. If you did not provide a starting guess, this is a bug. Please contact the burnman team.')

    if any(f < 0. for f in c[0]):
        raise Exception('Equilibration failed as a result of a phase amount being < 0. This might indicate that the starting guess is poor, or that there is no solution to this problem')

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

    # Add endmembers to the partial_gibbs vector only if they fall within the
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

def equilibrate(composition, assemblage, constraints, guesses=None):
    """
    Sets up a equilibration problem at constant composition and attempts to solve it

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
        List has the same form as PTX in :func:`_set_eqns`
        
    Returns
    -------
    sol_dict : dictionary of floats
        Dictionary contains the solution vector of the minimization.
        Dictionary keys are 'P', 'T' and the variable names output by
        :func:`compositional_variables`
    """

    # Redefine composition removing negligible elements (where n atoms < 0.0001 % total atoms)
    # Set up the inputs to the gibbs minimizer 
    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}
    col, null, initial_composition, indices, endmembers_per_phase = assemble_compositional_tensors ( composition, assemblage.phases, constraints )
        
    return _equilibrate(assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints, guesses)

def _equilibrate(assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints, guesses=None):

    if guesses == None:
        # Make an initial guess
        new_guesses = [10.e9, 1200.]
        new_guesses.extend(initial_composition)
    else:    
        # Remove unstable endmembers from guesses
        new_guesses=guesses[0:2] # keep P and T
        i_guess=2
        for i_phase, n_mbrs in enumerate(endmembers_per_phase):
            new_guesses.append(guesses[i_guess])
            i_mbrs=[j for i, j in indices if i==i_phase][1:]
            new_guesses.extend([guesses[i_guess + i_mbr] for i_mbr in i_mbrs])
            i_guess += n_mbrs

    # If P or T are fixed constraints, make sure that the initial guess is equal to that constraint!
    for constraint in constraints:
        if constraint[0] == 'P':
            new_guesses[0] = constraint[1]
        if constraint[1] == 'T':
            new_guesses[1] = constraint[1]

    # Set up the problem and attempt to solve it
    # Ignore warning due to changing the phase fractions such that they don't equal 1. They're renormalised anyway.
    soln=opt.fsolve(_set_eqns, new_guesses,
                    args=(assemblage, endmembers_per_phase,
                          initial_composition, col, null,
                          indices, constraints),
                    xtol=1e-10, full_output=True)
    
    if soln[2]==1:
        sol_dict = {'P': soln[0][0], 'T': soln[0][1], 'c': soln[0][2:2+len(indices)], 'variables': compositional_variables(assemblage, indices)}
        return sol_dict
    else:
        raise Exception('Solution could not be found, error: \n'+soln[3]+'\n Guesses:'
                        +str(guesses))

def find_equilibrium_composition(composition1, composition2, guessed_bulk, assemblage, constraints, guesses=None):
    """
    Sets up a equilibration problem with variable composition and attempts to solve it

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
        List has the same form as PTX in :func:`_set_eqns`
        
    Returns
    -------
    sol_dict : dictionary of floats
        Dictionary contains the solution vector of the minimization.
        Dictionary keys are 'P', 'T', 'X' and the variable names
        output by :func:`compositional_variables`
    """
    
    def minimize(x, composition1, composition2, assemblage, constraints, master_constraint, guesses):
        composition = binary_composition(composition1, composition2, x[0])
        vecsol = equilibrate(composition, assemblage, constraints, guesses)
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
        sol_dict = equilibrate(composition, assemblage, new_constraints, guesses)
        sol_dict['X'] = soln[0][0]
        return sol_dict
    else:
        raise Exception('Solution could not be found, error: '+soln[3])
        

def equilibrate_invariant(composition, phases, zero_phases, guesses=None):
    """
    Sets up a invariant equilibration problem and attempts to solve it

    Parameters
    ----------
    composition : dictionary of floats
        Dictionary contains the number of atoms of each element 
        in the bulk composition.
        
    phases : list of Minerals
        The assemblage of minerals for which we want to find
        amounts and compositions.

    zero_phases : list of Mineral
        The two minerals which becomes unstable at the invariant.

    guesses : optional list of floats
        List has the same form as PTX in :func:`_set_eqns`.
        
    Returns
    -------
    sol_array : List of floats
        List of P, T, compositional variables for the problem.

    """
    
    all_phases = zero_phases
    all_phases.extend(phases)
    assemblage = Composite(all_phases)
    base_assemblage = Composite(phases)

    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}

    
    formulae, indices, endmembers_per_phase = _get_formulae_indices_endmembers(composition, assemblage.phases)
    
    c0a = [0. for index in indices]
    c0b = [0. for index in indices]
    c1 = [1.]
    c1.extend([float(t - s) for s, t in zip(zip(*indices)[0], zip(*indices)[0][1:])])
    c0a[0] = 1.
    c0b[c1.index(1., 1)] = 1.

    constraints= [['X', c0a, c1, 0.], ['X', c0b, c1, 0.]]
    
    col, null, initial_composition, indices, endmembers_per_phase = assemble_compositional_tensors ( composition, assemblage.phases, constraints )
    
    sol = _equilibrate(assemblage, endmembers_per_phase, initial_composition, col, null, indices, constraints, guesses)
    sol_array = [sol['P'], sol['T']]
    sol_array.extend(sol['c'])
    sol_array.extend(sol['variables'])
    return sol_array
    
def equilibrate_univariant(composition, phases, zero_phase, condition_variable, condition_array, guesses=None):
    """
    Sets up a univariant equilibration problem and attempts to solve it

    Parameters
    ----------
    composition : dictionary of floats
        Dictionary contains the number of atoms of each element 
        in the bulk composition.
        
    phases : list of Minerals
        The assemblage of minerals for which we want to find
        amounts and compositions.

    zero_phase : Mineral
        The mineral which becomes unstable along the univariant.

    condition_variable : string
        Either 'P' or 'T', depending on the extensive quantity to loop over.

    condition_array : list (or array) of floats 
        List of the condition variable to loop over.

    guesses : optional list of floats
        List has the same form as PTX in :func:`_set_eqns`.
        
    Returns
    -------
    sol_array : List of list of floats
        List of lists of P, T, compositional variables for the problem. 
        If a solution is not found for one of the values in condition_array
        nothing is appended to sol_array.

    """

    all_phases = [zero_phase]
    all_phases.extend(phases)
    assemblage = Composite(all_phases)
    
    s = sum(composition.values())
    composition = {element:amount for element, amount in composition.items() if amount/s > 1.e-6}

    formulae, indices, endmembers_per_phase = _get_formulae_indices_endmembers(composition, assemblage.phases)
    
    c0 = [0. for index in indices]
    c0[0] = 1.
    c1 = [1.]
    c1.extend([float(t - s) for s, t in zip(zip(*indices)[0], zip(*indices)[0][1:])])

    X_constraint = ['X', c0, c1, 0.]
    col, null, initial_composition, indices, endmembers_per_phase = assemble_compositional_tensors ( composition, assemblage.phases, [X_constraint] )

    sol_array = []
    if condition_variable=='P':
        for i, P in enumerate(condition_array):
            constraints= [['P', P], X_constraint]
            try:
                sol = _equilibrate(assemblage, endmembers_per_phase, initial_composition,
                                       col, null, indices, constraints, guesses)
                guesses = [sol['P'], sol['T']]
                guesses.extend(sol['c'])
                sol_array.append(guesses)
            except:
                print('No solution found at {0:.5f} GPa'.format(P/1.e9))
    elif condition_variable=='T':
        for i, T in enumerate(condition_array):
            constraints= [['T', T], X_constraint]
            try:
                sol = _equilibrate(assemblage, endmembers_per_phase, initial_composition,
                                       col, null, indices, constraints, guesses)
                guesses = [sol['P'], sol['T']]
                guesses.extend(sol['c'])
                sol_array.append(guesses)
                sol_array[-1].extend(sol['variables'])
            except:
                print('No solution found at {0:.1f} K'.format(T))

    return sol_array
