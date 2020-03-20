# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2019 by the BurnMan team, released under the GNU
# GPL v2 or later.
from __future__ import absolute_import
from __future__ import print_function
import numpy as np

from scipy.optimize import linprog
from scipy.linalg import lu_factor, lu_solve

from .nonlinear_solvers import damped_newton_solve
from collections import namedtuple

from .mineral import Mineral
from .solidsolution import SolidSolution
# from .solutionbases import feasible_solution_basis_in_component_space
# from .solutionbases import transform_solution_to_new_basis

from sympy import Matrix, nsimplify


def columnspace_with_pivots(matrix):
    """
    Returns the columnspace including the column indices of the original matrix
    """
    reduced, pivots = matrix.rref()

    basis = []
    # create a set of vectors for the basis
    for i in range(matrix.cols):
        if i in pivots:
            basis.append(matrix.col(i))
    return [matrix._new(b) for b in basis], pivots


def remove_redundant_inequalities(inequalities_matrix, inequalities_vector):
    # Set of existing inequalities
    # Ax <= -b
    # We want to test if there is any solution to the additional constraint
    # Sx >= -t
    # If not, that constraint is redundant

    # Set up bounds and null vector for the linear programming problem
    n_dim = len(inequalities_matrix[0])
    bounds = [[None, None] for j in range(n_dim)]
    c = np.zeros(n_dim)

    # First, test that there is a value of x
    # where all constraints are satisfied!
    lp = linprog(c, A_ub=inequalities_matrix,
                 b_ub=-inequalities_vector, bounds=bounds)
    assert(lp.success)

    A = np.copy(inequalities_matrix)
    b = np.copy(inequalities_vector)

    # Sweep forward and backward to remove all redundant constraints
    for dirn in ['forward', 'backward']:
        A_build = np.copy(A[0, np.newaxis])
        b_build = np.ones(1)*b[0]
        for i in range(1, len(A)):
            A_ub = np.concatenate((A_build, -A[i, np.newaxis]))
            b_ub = np.concatenate((b_build, -b[i, np.newaxis]+1.e-10))

            lp = linprog(c, A_ub=A_ub, b_ub=-b_ub, bounds=bounds)

            if lp.success:
                A_build = np.concatenate((A_build, A[i, np.newaxis]))
                b_build = np.concatenate((b_build, b[i, np.newaxis]))

        A = np.copy(A_build[::-1])
        b = np.copy(b_build[::-1])

    return (A, b)


def simplify_matrix(arr):
    def f(i, j):
        return nsimplify(arr[i][j])
    return Matrix(len(arr), len(arr[0]), f)


def get_formulae_and_nendmembers(assemblage, elements):
    # Listify the elements and sort them so they have a consistent order.
    # This will be the ordering of the rows.  The ordering of the columns
    # will be the ordering of the endmembers as they are passed in.

    formulae = []
    endmembers_per_phase = []
    for ph_idx, ph in enumerate(assemblage.phases):

        if isinstance(ph, SolidSolution):
            formulae.extend(ph.endmember_formulae)
            endmembers_per_phase.append(ph.n_endmembers)

        elif isinstance(ph, Mineral):
            formulae.append(ph.formula)
            endmembers_per_phase.append(1)

        else:
            raise Exception('Unsupported mineral type, '
                            'can only read burnman.Mineral '
                            'or burnman.SolidSolution')

    return formulae, endmembers_per_phase


def calculate_constraints(assemblage, n_phase_mbrs):
    """

    """
    n_mbrs = sum(n_phase_mbrs)

    bounds = []
    n_constraints = 0
    n_raw_constraints = 0
    for i, n in enumerate(n_phase_mbrs):
        n_constraints += 1
        if n == 1:
            bounds.append(np.array([[]]))
            n_raw_constraints += 1
        else:
            bounds.append(assemblage.phases[i].solution_model.endmember_occupancies)
            n_constraints += len(bounds[-1][0])
            n_raw_constraints += len(bounds[-1][0])

    c_vector = np.zeros((n_constraints+2))
    c_matrix = np.zeros((n_constraints+2, n_mbrs+2))  # includes P, T
    raw_c_matrix = np.zeros((n_raw_constraints, n_mbrs))

    c_matrix[0, 0] = -1  # P>0
    c_matrix[1, 1] = -1  # T>0

    rcidx = 0  # index of current raw compositional constraint
    cidx = 2  # index of current compositional constraint
    pidx = 0  # starting index of current phase
    for i, n in enumerate(n_phase_mbrs):
        m = len(bounds[i][0])
        # The first endmember proportion is not a free variable
        # (all endmembers in any solution must sum to one)
        # Re-express the constraints without the first endmember
        c_matrix[cidx, pidx+2] = -1.  # need phase proportions > 0
        cidx += 1
        if m == 0:
            raw_c_matrix[rcidx, pidx] = 1.
            rcidx += 1
        else:
            c_vector[cidx:cidx+m] = -bounds[i][0]
            c_matrix[cidx:cidx+m,
                     pidx+1+2:pidx+n+2] = (np.einsum('i, j',
                                                     bounds[i][0],
                                                     np.ones_like(bounds[i][1:, 0]))
                                           - bounds[i].T[:, 1:])
            raw_c_matrix[rcidx:rcidx+m, pidx:pidx+n] = bounds[i].T
            cidx += m
        rcidx += m
        pidx += n

    c_matrix, c_vector = remove_redundant_inequalities(c_matrix, c_vector)

    return c_vector, c_matrix, raw_c_matrix


def get_parameters(assemblage):
    n_moles_phase = assemblage.n_moles * np.array(assemblage.molar_fractions)

    try:
        params = [assemblage.pressure, assemblage.temperature]
    except AttributeError:
        raise Exception('You need to set_state before getting parameters')

    for i, ph in enumerate(assemblage.phases):
        params.append(n_moles_phase[i])
        if isinstance(ph, SolidSolution):
            params.extend([f for f
                           in assemblage.phases[i].molar_fractions[1:]])
    return np.array(params)


def get_endmember_amounts(assemblage):
    phase_amounts = assemblage.n_moles*np.array(assemblage.molar_fractions)
    amounts = []
    for i, ph in enumerate(assemblage.phases):
        if isinstance(ph, SolidSolution):
            amounts.extend(phase_amounts[i]
                           * assemblage.phases[i].molar_fractions)
        else:
            amounts.append(phase_amounts[i])
    return np.array(amounts)


def set_compositions_and_state_from_parameters(assemblage, parameters):
    p = parameters
    assemblage.set_state(p[0], p[1])
    i = 2
    phase_amounts = np.zeros(len(assemblage.phases))
    for phase_idx, ph in enumerate(assemblage.phases):
        phase_amounts[phase_idx] = p[i]
        if isinstance(ph, SolidSolution):
            n_mbrs = len(ph.endmembers)
            f = [0.]*n_mbrs
            f[1:] = p[i+1:i+n_mbrs]
            f[0] = 1. - sum(f)
            ph.set_composition(f)
            i += n_mbrs
        else:
            i += 1

    assert(np.all(phase_amounts > -1.e-8))
    phase_amounts = np.abs(phase_amounts)
    assemblage.n_moles = sum(phase_amounts)
    assemblage.set_fractions(phase_amounts/assemblage.n_moles)
    return None


def F(x, assemblage, equality_constraints, prm):
    """
    x is a list of the pressure, temperature and compositional parameters
    (in that order).

    If equality_constraints[i][0] = 'P', F[i] = P - equality_constraints[i][1]
    If equality_constraints[i][0] = 'T', F[i] = T - equality_constraints[i][1]

    If one of equality_constraints[i][0] = 'PT',
    then the constraint is equivalent to
    the P-T point lying on an ellipse centered on a given point.

    For each equality_constraints = X, then an extra
    """

    set_compositions_and_state_from_parameters(assemblage, x)

    new_endmember_amounts = get_endmember_amounts(assemblage)

    n_reactions = len(prm.stoic_nullspace[:, 0])
    partial_gibbs_vector = np.zeros((prm.n_endmembers))
    j = 0
    for i, n_mbr in enumerate(prm.endmembers_per_phase):
        ph = assemblage.phases[i]
        if n_mbr == 1:  # for endmembers
            partial_gibbs_vector[j] = ph.gibbs
        else:  # for solid solutions
            partial_gibbs_vector[j:j+n_mbr] = ph.partial_gibbs
        j += n_mbr

    # We want to find the root of the following equations
    eqns = np.zeros((2 + prm.n_endmembers))
    for i, (type_c, eq_c) in enumerate(equality_constraints):
        if type_c == 'P':
            eqns[i] = x[0] - eq_c
        elif type_c == 'T':
            eqns[i] = x[1] - eq_c
        elif type_c == 'S':
            eqns[i] = assemblage.molar_entropy*assemblage.n_moles - eq_c
        elif type_c == 'V':
            eqns[i] = assemblage.molar_volume*assemblage.n_moles - eq_c
        elif type_c == 'PT_ellipse':
            v_scaled = (x[0:2] - eq_c[0])/eq_c[1]
            eqns[i] = np.linalg.norm(v_scaled) - 1.
        elif type_c == 'X':
            eqns[i] = np.dot(eq_c[0], x) - eq_c[1]  # i.e. Ax = b
        else:
            raise Exception('constraint type not recognised')

    eqns[2:2+n_reactions] = np.dot(prm.stoic_nullspace, partial_gibbs_vector)
    eqns[2+n_reactions:] = (np.dot(prm.stoic_colspace, new_endmember_amounts)
                            - prm.reduced_composition_vector)
    return eqns


def jacobian(x, assemblage, equality_constraints, prm):
    # The solver always calls the Jacobian with the same
    # x parameters as used previously for the root functions
    # Therefore we don't need to set compositions or state again here
    # If the next two lines are uncommented, do it anyway.

    # the jacobian J = dfi/dxj

    # First, we find out the effect of the two constraint parameters on the
    # pressure and temperature functions:
    # i.e. dF(i, constraints)/dx[0, 1]
    full_hessian = np.zeros((prm.n_endmembers+2, prm.n_endmembers+2))
    for i, (type_c, eq_c) in enumerate(equality_constraints):
        if type_c == 'P':  # dP/dx
            full_hessian[i, 0] = 1.  # full_hessian[i, j!=0] = 0
        elif type_c == 'T':  # dT/dx
            full_hessian[i, 1] = 1.  # full_hessian[i, j!=1] = 0
        elif type_c == 'S':  # dS/dx
            # dS/dP = -aV, dS/dT = Cp/T
            full_hessian[i, 0:2] = [-assemblage.n_moles
                                    * assemblage.alpha
                                    * assemblage.molar_volume,
                                    assemblage.n_moles
                                    * assemblage.molar_heat_capacity_p / x[1]]
            j = 2
            for k, n in enumerate(prm.endmembers_per_phase):
                full_hessian[i, j] = assemblage.phases[k].molar_entropy
                if n > 1:  # for solutions with >1 endmember
                    full_hessian[i, j+1:j+n] = (assemblage.n_moles
                                                * assemblage.molar_fractions[k]
                                                * (assemblage.phases[k].partial_entropies[1:]
                                                   - assemblage.phases[k].partial_entropies[0]))
                j += n
        elif type_c == 'V':  # dV/dx
            # dV/dP = -V/K_T, dV/dT = aV
            full_hessian[i, 0:2] = [-assemblage.n_moles
                                    * assemblage.molar_volume / assemblage.K_T,
                                    assemblage.n_moles*assemblage.molar_volume]
            j = 2
            for k, n in enumerate(prm.endmembers_per_phase):
                full_hessian[i, j] = assemblage.phases[k].molar_volume
                if n > 1:  # for solutions with >1 stable endmember
                    full_hessian[i, j+1:j+n] = (assemblage.n_moles
                                                * assemblage.molar_fractions[k]
                                                * (assemblage.phases[k].partial_volumes[1:]
                                                   - assemblage.phases[k].partial_volumes[0]))
                j += n
        elif type_c == 'PT_ellipse':
            v_scaled = (x[0:2] - eq_c[0])/eq_c[1]
            full_hessian[i, 0:2] = v_scaled/(np.linalg.norm(v_scaled)*eq_c[1])
        elif type_c == 'X':
            full_hessian[i, :] = eq_c[0]
        else:
            raise Exception('constraint type not recognised')

    # Next, let's get the effect of pressure and temperature
    # on each of the independent reactions
    # i.e. dF(i, reactions)/dx[0] and dF(i, reactions)/dx[1]
    partial_volumes_vector = np.zeros((prm.n_endmembers))
    partial_entropies_vector = np.zeros((prm.n_endmembers))
    j = 0
    for i, n in enumerate(prm.endmembers_per_phase):
        if n == 1:  # for endmembers
            partial_volumes_vector[j] = assemblage.phases[i].molar_volume
            partial_entropies_vector[j] = assemblage.phases[i].molar_entropy
        else:  # for solid solutions
            partial_volumes_vector[j:j+n] = assemblage.phases[i].partial_volumes
            partial_entropies_vector[j:j+n] = assemblage.phases[i].partial_entropies
        j += n
    reaction_volumes = np.dot(prm.stoic_nullspace, partial_volumes_vector)
    reaction_entropies = np.dot(prm.stoic_nullspace, partial_entropies_vector)

    # dGi/dP = deltaVi; dGi/dT = -deltaSi
    full_hessian[2:2+len(reaction_volumes), 0] = reaction_volumes
    full_hessian[2:2+len(reaction_volumes), 1] = -reaction_entropies

    # Pressure and temperature have no effect on the bulk
    # compositional constraints
    # i.e. dF(i, bulk)/dx[0] and dF(i, bulk)/dx[1] = 0

    # Finally, let's build the compositional Hessian d2G/dfidfj = dmui/dfj
    # where fj is the fraction of endmember j in a phase
    phase_amounts = np.array(assemblage.molar_fractions) * assemblage.n_moles
    comp_hessian = np.zeros((prm.n_endmembers, prm.n_endmembers))
    dfi_dxj = np.zeros((prm.n_endmembers, prm.n_endmembers))
    dpi_dxj = np.zeros((prm.n_endmembers, prm.n_endmembers))
    j = 0
    for i, n in enumerate(prm.endmembers_per_phase):
        # ignore pure phases, as changing proportions of a phase
        # does not change its molar gibbs free energy
        if n == 1:
            dpi_dxj[j, j] = 1.
        else:
            comp_hessian[j:j+n, j:j+n] = assemblage.phases[i].gibbs_hessian
            dfi_dxj[j:j+n, j:j+n] = np.eye(n)
            dfi_dxj[j, j:j+n] -= 1.  # x[0] = p(phase) and x[1:] = f[1:] - f[0]

            dpi_dxj[j:j+n, j:j+n] = dfi_dxj[j:j+n, j:j+n] * phase_amounts[i]
            dpi_dxj[j:j+n, j] = np.array(assemblage.phases[i].molar_fractions)
        j += n

    # dfi_dxj converts the endmember hessian to the parameter hessian.
    reaction_hessian = prm.stoic_nullspace.dot(comp_hessian).dot(dfi_dxj)
    bulk_hessian = prm.stoic_colspace.dot(dpi_dxj)
    full_hessian[2:, 2:] = np.concatenate((reaction_hessian, bulk_hessian))
    return full_hessian


def lambda_bounds(dx, x, endmembers_per_phase):
    """
    x is a list of the pressure, temperature and compositional parameters
    (in that order).

    dx is the proposed newton step
    """

    max_steps = np.ones((len(x)))*100000.

    # first two constraints are P and T
    max_steps[0:2] = [20.e9, 500.]  # biggest reasonable P and T steps

    j = 2
    for i, n in enumerate(endmembers_per_phase):
        # if the phase proportion constraint would be broken,
        # set a step that is marginally smaller
        if x[j] + dx[j] < 0.:
            max_steps[j] = max(x[j]*0.999, 0.001)
            #max_steps[j] = max(x[j]*0.5, 0.001) # allowed to cut amount by half...
        max_steps[j+1:j+n] = [max(xi*0.99, 0.01) for xi in x[j+1:j+n]]  # maximum compositional step
        j += n

    max_lmda = min([1. if step <= max_steps[i] else max_steps[i]/step
                    for i, step in enumerate(np.abs(dx))])

    # return (1.e-8, 1.)
    return (1.e-8, max_lmda)


# The following two functions construct two common examples of
# compositional constraints. Note that constraints in the form
# sum(ax + b)/sum(cx + d) - f = 0
# can be recast as:
# (a-fc)*x = fd - b
# which is less easy for a human to understand
# (in terms of chemical constraints), but much easier to solve.
def phase_fraction_constraint(phase, assemblage, proportion,
                                prm):
    phase_idx = assemblage.phases.index(phase)

    constraints = []
    for p in proportion:
        constraints.append(['X', [np.zeros((prm.n_parameters)), 0.]])
        constraints[-1][-1][0][prm.phase_fraction_indices] = -p
        constraints[-1][-1][0][prm.phase_fraction_indices[phase_idx]] += 1.

    return constraints


def phase_composition_constraint(phase, assemblage, indices, constraint):
    phase_idx = assemblage.phases.index(phase)
    start_idx = int(sum([len(i) for i in indices[:phase_idx]])) + 3 # +3 comes from P, T and the proportion of the phase of interest
    n_indices = sum([len(i) for i in indices])
    mbr_indices = indices[phase_idx]

    site_names, numerator, denominator, value = constraint
    site_indices = [phase.solution_model.site_names.index(name) for name in site_names]
    atoms = np.dot(phase.solution_model.endmember_noccupancies[mbr_indices][:, site_indices], np.array([numerator, denominator]).T)

    atoms0 = atoms[0]
    atoms -= atoms[0]
    numer, denom = atoms.T[:, 1:]

    constraints = []
    for v in value:
        f = v*atoms0[1] - atoms0[0]
        constraints.append(['X', [np.zeros((n_indices+2)), f]])
        constraints[-1][1][0][start_idx:start_idx+len(mbr_indices)-1] = numer - v*denom

    print(constraints)
    return constraints


def get_equilibration_parameters(assemblage, composition):
    # Initialize a named tuple for the equilibration parameters
    prm = namedtuple('assemblage_parameters', [])

    # Process elements
    prm.elements = sorted(set(composition.keys()))
    fn = get_formulae_and_nendmembers(assemblage, prm.elements)
    prm.formulae, prm.endmembers_per_phase = fn
    prm.n_endmembers = sum(prm.endmembers_per_phase)

    # Process parameter names
    prm.parameter_names = ['Pressure (Pa)', 'Temperature (K)']
    for i, n in enumerate(prm.endmembers_per_phase):
        prm.parameter_names.append('x({0})'.format(assemblage.phases[i].name))
        if n > 1:
            p_names = ['p({0} in {1})'.format(n, assemblage.phases[i].name)
                       for n in assemblage.phases[i].endmember_names[1:]]
            prm.parameter_names.extend(p_names)

    prm.n_parameters = len(prm.parameter_names)
    prm.phase_fraction_indices = [i for i in range(len(prm.parameter_names))
                                  if 'x(' in prm.parameter_names[i]]

    # Find the bulk composition vector
    prm.bulk_composition_vector = np.array([composition[e]
                                            for e in prm.elements])

    # Populate the stoichiometric matrix
    def f(i, j):
        e = prm.elements[i]
        if e in prm.formulae[j]:
            return nsimplify(prm.formulae[j][e])
        else:
            return 0
    S = Matrix(len(prm.elements), len(prm.formulae), f)
    prm.stoichiometric_matrix = S

    cspace, pivots = columnspace_with_pivots(S.T)
    prm.stoic_colspace = np.array(cspace)
    prm.stoic_nullspace = np.array(S.nullspace())
    right_nullspace = np.array(S.T.nullspace())
    if right_nullspace.shape[0] != 0:
        if (np.abs(right_nullspace.dot(prm.bulk_composition_vector)[0]) > 1.e-12):
            raise Exception('The bulk composition is not within the '
                            'compositional space of the assemblage')

    prm.reduced_composition_vector = [prm.bulk_composition_vector[i]
                                      for i in pivots]

    prm.constraint_vector, prm.constraint_matrix, prm.raw_constraint_matrix = calculate_constraints(assemblage, prm.endmembers_per_phase)
    return prm


def process_eq_constraints(equality_constraints, assemblage, prm):
    eq_constraint_lists = []
    for i in range(2):
        if equality_constraints[i][0] == 'phase_proportion':
            phase = equality_constraints[i][1][0]
            proportion = equality_constraints[i][1][1]
            if isinstance(proportion, float):
                proportion = np.array([proportion])
            if not isinstance(proportion, np.ndarray):
                raise Exception('The constraint proportion in equality {0} '
                                'should be a float or numpy array'.format(i+1))
            eq_constraint_lists.append(phase_fraction_constraint(phase,
                                                                 assemblage,
                                                                 proportion,
                                                                 prm))
        elif equality_constraints[i][0] == 'phase_composition':
            phase = equality_constraints[i][1][0]
            constraint = equality_constraints[i][1][1]
            if isinstance(constraint[3], float):
                constraint = (constraint[0], constraint[1], constraint[2],
                              np.array([constraint[3]]))
            if not isinstance(constraint[3], np.ndarray):
                raise Exception('The last constraint parameter in equality '
                                '{0} should be a float '
                                'or numpy array'.format(i+1))

            eq_constraint_lists.append(phase_composition_constraint(phase,
                                                                    assemblage,
                                                                    constraint))
        elif equality_constraints[i][0] == 'X':
            constraint = equality_constraints[i][1]
            if isinstance(constraint[-1], float):
                constraint = (constraint[0], np.array([constraint[-1]]))
            if not isinstance(constraint[-1], np.ndarray):
                raise Exception('The last constraint parameter in '
                                'equality {0} should be '
                                'a float or numpy array'.format(i+1))

            eq_constraint_lists.append([['X', [constraint[0], p]]
                                        for p in constraint[1]])

        elif (equality_constraints[i][0] == 'P'
              or equality_constraints[i][0] == 'T'
              or equality_constraints[i][0] == 'PT_ellipse'
              or equality_constraints[i][0] == 'S'
              or equality_constraints[i][0] == 'V'):
            if isinstance(equality_constraints[i][1], float):
                equality_constraints[i] = (equality_constraints[i][0],
                                           np.array([equality_constraints[i][1]]))
            if not isinstance(equality_constraints[i][1], np.ndarray):
                raise Exception('The last parameter in equality {0} should be '
                                'a float or numpy array'.format(i+1))
            eq_constraint_lists.append([[equality_constraints[i][0], p]
                                        for p in equality_constraints[i][1]])
        else:
            raise Exception('The type of equality_constraint is '
                            'not recognised for constraint {0}.\n'
                            'Should be one of P, T, S, V, X,\n'
                            'PT_ellipse, phase_proportion, '
                            'or phase_composition.'.format(i+1))
    return eq_constraint_lists


def equilibrate(composition, assemblage, equality_constraints,
                tol=1.e-3,
                store_iterates=False, store_assemblage=False,
                max_iterations=100., verbose=True):

    for ph in assemblage.phases:
        if isinstance(ph, SolidSolution):
            if not hasattr(ph, 'molar_fractions'):
                ph.set_composition(ph.guess)

    if assemblage.molar_fractions is None:
        n_phases = len(assemblage.phases)
        f = 1./float(n_phases)
        assemblage.set_fractions([f for i in range(n_phases)])

    assemblage.n_moles = (sum(composition.values())
                          / sum(assemblage.formula.values()))

    # Make parameter tuple
    prm = get_equilibration_parameters(assemblage, composition)

    # Check equality constraints have the correct structure
    # Convert into the format readable by the function and jacobian functions
    eq_constraint_lists = process_eq_constraints(equality_constraints,
                                                 assemblage, prm)

    # Set up solves
    sol_list = []
    n_c0 = len(eq_constraint_lists[0])
    n_c1 = len(eq_constraint_lists[1])

    # Find the initial state (could be none here)
    initial_state = [assemblage.pressure, assemblage.temperature]

    # Reset initial state if  equality constraints
    # are related to pressure or temperature
    for i in range(2):
        if eq_constraint_lists[i][0][0] == 'P':
            initial_state[0] = eq_constraint_lists[i][0][1]
        elif eq_constraint_lists[i][0][0] == 'T':
            initial_state[1] = eq_constraint_lists[i][0][1]
        elif eq_constraint_lists[i][0][0] == 'PT_ellipse':
            initial_state = eq_constraint_lists[i][0][1][1]

    if initial_state[0] is None:
        initial_state[0] = 5.e9
    if initial_state[1] is None:
        initial_state[1] = 1200.

    assemblage.set_state(*initial_state)
    prm.initial_parameters = get_parameters(assemblage)
    parameters = get_parameters(assemblage)

    # Solve the system of equations, loop over input parameters
    sol_list = np.empty(shape=(n_c0, n_c1)+(0,)).tolist()
    for i_c0 in range(n_c0):
        new_c0 = True
        for i_c1 in range(n_c1):
            if verbose:
                string = 'Processing solution'
                if n_c0 > 1:
                    string += ' {0}/{1}'.format(i_c0+1, n_c0)
                if n_c1 > 1:
                    string += ' {0}/{1}'.format(i_c1+1, n_c1)
                print(string+':')

            equality_constraints = [eq_constraint_lists[0][i_c0],
                                    eq_constraint_lists[1][i_c1]]

            # Set the initial proportions and compositions
            # of the phases in the assemblage:
            # set_compositions_and_state_from_parameters(assemblage,
            #                                           prm.initial_parameters)
            try:
                sol = damped_newton_solve(F=lambda x: F(x, assemblage,
                                                        equality_constraints, prm),
                                          J=lambda x: jacobian(x, assemblage,
                                                               equality_constraints,
                                                               prm),
                                          lambda_bounds=lambda dx, x: lambda_bounds(dx, x, prm.endmembers_per_phase),
                                          guess=parameters,
                                          linear_constraints=(prm.constraint_matrix,
                                                              prm.constraint_vector),
                                          tol=tol,
                                          store_iterates=store_iterates,
                                          max_iterations=max_iterations)
            except:
                pass

            if store_assemblage:
                sol.assemblage = assemblage.copy()

            if verbose:
                print(sol.text)

            sol_list[i_c0][i_c1] = sol
            new_c0 = False

            prev = []
            if i_c1 < n_c1 - 1:
                next_ecs = [i_c0, i_c1 + 1]
            elif i_c1 == n_c1 - 1 and i_c0 < n_c0 - 1:
                next_ecs = [i_c0+1, 0]
            else:  # last value
                next_ecs = None

            if next_ecs is not None:
                cs = [eq_constraint_lists[0][next_ecs[0]],
                      eq_constraint_lists[1][next_ecs[1]]]
                prev_sol = []
                if next_ecs[0] != 0:
                    prev_sol.append(sol_list[next_ecs[0] - 1][next_ecs[1]])
                if next_ecs[1] != 0:
                    prev_sol.append(sol_list[next_ecs[0]][next_ecs[1] - 1])

                updated_params = False
                for s in prev_sol:
                    if s.success and not updated_params:
                        # next guess based on a Newton step
                        dF = F(s.x, assemblage, cs, prm)
                        luJ = lu_factor(s.J)
                        new_parameters = s.x + lu_solve(luJ, -dF)
                        c = (prm.constraint_matrix.dot(new_parameters)
                             + prm.constraint_vector)
                        if all(c <= 0.):  # accept new guess
                            parameters = new_parameters
                        else:  # use the parameters from this step
                            parameters = s.x
                            exhausted_phases = [assemblage.phases[phase_idx].name
                                                for phase_idx, v in
                                                enumerate(new_parameters[prm.phase_fraction_indices]) if v<0.]
                            if len(exhausted_phases) > 0 and verbose:
                                print('A phase might be exhausted before the next step: {0}'.format(exhausted_phases))

                        updated_params = True
                if not updated_params:
                    parameters = prm.initial_parameters

    # Finally, make dimensions of sol_list equal the input dimensions
    if len(sol_list[0]) == 1:
        sol_list = list(zip(*sol_list))[0]
    if len(sol_list) == 1:
        sol_list = sol_list[0]
    return sol_list, prm
