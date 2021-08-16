# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np

from scipy.linalg import lu_factor, lu_solve
from collections import namedtuple

from .nonlinear_solvers import damped_newton_solve
from .solidsolution import SolidSolution


def calculate_constraints(assemblage):
    """
    This function calculates the linear inequality constraints bounding
    the valid parameter space for a given assemblage.

    The constraints are as follows:

      - Pressure and temperature must be positive
      - All phase fractions must be positive
      - All site-species occupancies must be positive

    The constraints are stored in a vector (b) and matrix (A).
    The sign convention is chosen such that the constraint is satisfied
    if A.x + b < eps.

    Parameters
    ----------
    assemblage : burnman.Composite object
        The assemblage for which the constraints are calculated.

    Returns
    -------
    c_vector :
        The constraints vector.

    c_matrix :
        The constraints matrix.
    """
    bounds = []
    n_constraints = 0
    for i, n in enumerate(assemblage.endmembers_per_phase):
        n_constraints += 1
        if n == 1:
            bounds.append(np.array([[]]))
        else:
            bounds.append(assemblage.phases[i].solution_model.endmember_occupancies)
            n_constraints += len(bounds[-1][0])

    c_vector = np.zeros((n_constraints+2))
    c_matrix = np.zeros((n_constraints+2, assemblage.n_endmembers+2))  # includes P, T

    c_matrix[0, 0] = -1  # P>0
    c_matrix[1, 1] = -1  # T>0

    cidx = 2  # index of current compositional constraint
    pidx = 0  # starting index of current phase
    for i, n in enumerate(assemblage.endmembers_per_phase):
        m = len(bounds[i][0])
        # The first endmember proportion is not a free variable
        # (all endmembers in any solution must sum to one)
        # Re-express the constraints without the first endmember
        c_matrix[cidx, pidx+2] = -1.  # need phase proportions > 0
        cidx += 1
        if m != 0:
            c_vector[cidx:cidx+m] = -bounds[i][0]
            c_matrix[cidx:cidx+m,
                     pidx+1+2:pidx+n+2] = (np.einsum('i, j', bounds[i][0],
                                                     np.ones_like(bounds[i][1:, 0]))
                                           - bounds[i].T[:, 1:])
            cidx += m
        pidx += n

    return c_vector, c_matrix


def get_parameters(assemblage):
    """
    Gets the parameters vector  for the current equilibrium problem.
    These are:

      - pressure
      - temperature
      - absolute amount of each phase
      - if that phase is a solution with >1 endmember, the mole fractions of
      the second to the last independent endmember in the solution.

    Parameters
    ----------
    assemblage : burnman.Composite object
        The assemblage for which equilibrium is to be calculated.

    Returns
    -------
    params : numpy array
        An array containing the current parameter values.
    """
    params = np.empty(assemblage.n_endmembers + 2)
    n_moles_phase = assemblage.n_moles * np.array(assemblage.molar_fractions)

    try:
        params[:2] = [assemblage.pressure, assemblage.temperature]
    except AttributeError:
        raise Exception('You need to set_state before getting parameters')

    j=2
    for i, ph in enumerate(assemblage.phases):
        params[j] = n_moles_phase[i]
        if isinstance(ph, SolidSolution):
            params[j+1:j+assemblage.endmembers_per_phase[i]] = assemblage.phases[i].molar_fractions[1:]
        j += assemblage.endmembers_per_phase[i]

    return params


def get_endmember_amounts(assemblage):
    """
    Gets the absolute amounts of all the endmembers in the solution.

    Parameters
    ----------
    assemblage : burnman.Composite object
        The assemblage for which equilibrium is to be calculated.

    Returns
    -------
    amounts : numpy array
        An array containing the current amounts of all the endmembers.
    """
    phase_amounts = assemblage.n_moles * assemblage.molar_fractions
    amounts = np.empty(assemblage.n_endmembers)
    j=0
    for i, ph in enumerate(assemblage.phases):
        if isinstance(ph, SolidSolution):
            amounts[j:j+assemblage.endmembers_per_phase[i]] = phase_amounts[i] * assemblage.phases[i].molar_fractions
        else:
            amounts[j] = phase_amounts[i]
        j += assemblage.endmembers_per_phase[i]

    return amounts


def set_compositions_and_state_from_parameters(assemblage, parameters):
    """
    Sets the phase compositions, amounts and state of the assemblage
    from a list of parameter values.

    Parameters
    ----------
    assemblage : burnman.Composite object
        The assemblage for which equilibrium is to be calculated.

    parameters : numpy array
        An array containing the current parameter values.
    """
    assemblage.set_state(parameters[0], parameters[1])
    i = 2
    phase_amounts = np.zeros(len(assemblage.phases))
    for phase_idx, ph in enumerate(assemblage.phases):
        phase_amounts[phase_idx] = parameters[i]
        if isinstance(ph, SolidSolution):
            n_mbrs = len(ph.endmembers)
            f = [0.]*n_mbrs
            f[1:] = parameters[i+1:i+n_mbrs]
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


def F(x, assemblage, equality_constraints, reduced_composition_vector):
    """
    Returns a vector of values which are zero at equilibrium.
    The equations corresponding to the first two vector values depend on the
    equality_constraints chosen. For example, if
      - eq[i][0] = 'P', F[i] = P - eq[i][1]
      - eq[i][0] = 'T', F[i] = T - eq[i][1]
      - eq[i][0] = 'S', F[i] = entropy - eq[i][1]
      - eq[i][0] = 'V', F[i] = volume - eq[i][1]
      - eq[i][0] = 'PT_ellipse', F[i] = norm(([P, T] - eq[i][1][0])/eq[i][1][1]) - 1
      - eq[i][0] = 'X', np.dot(eq[i][1][0], x) - eq[i][1][1]

    The next set of vector values correspond to the reaction affinities.
    The final set of vector values correspond to the bulk
    composition constraints.

    Parameters
    ----------
    x : numpy array
        Parameter values for the equilibrium problem to be solved.

    assemblage : burnman.Composite object
        The assemblage for which equilibrium is to be calculated.

    equality_constraints : list of lists
        A list of the equality constraints (see above).

    reduced_composition_vector : numpy array
        The vector corresponding to the amounts of the independent
        elements.

    Returns
    -------
    eqns : numpy array
        An array containing the vector values which
        are equal to zero at equilibrium.
    """

    set_compositions_and_state_from_parameters(assemblage, x)
    new_endmember_amounts = get_endmember_amounts(assemblage)

    # We want to find the root of the following equations
    eqns = np.zeros((2 + assemblage.n_endmembers))
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

    eqns[2:2+assemblage.n_reactions] = assemblage.reaction_affinities
    eqns[2+assemblage.n_reactions:] = (np.dot(assemblage.reduced_stoichiometric_array.T,
                                              new_endmember_amounts)
                                       - reduced_composition_vector)
    return eqns


def jacobian(x, assemblage, equality_constraints):
    """
    The Jacobian of the equilibrium problem (dF/dx).
    See documentation for F for more details.

    Parameters
    ----------
    x : numpy array
        Parameter values for the equilibrium problem to be solved.

    assemblage : burnman.Composite object
        The assemblage for which equilibrium is to be calculated.

    equality_constraints : list of lists
        A list of the equality constraints (see documentation for F).

    Returns
    -------
    jacobian : 2D numpy array
        An array containing the Jacobian for the equilibrium problem
    """
    # The solver always calls the Jacobian with the same
    # x parameters as used previously for the root functions
    # Therefore we don't need to set compositions or state again here
    # If the next two lines are uncommented, do it anyway.

    # the jacobian J = dfi/dxj

    # First, we find out the effect of the two constraint parameters on the
    # pressure and temperature functions:
    # i.e. dF(i, constraints)/dx[0, 1]
    jacobian = np.zeros((assemblage.n_endmembers+2, assemblage.n_endmembers+2))
    for i, (type_c, eq_c) in enumerate(equality_constraints):
        if type_c == 'P':  # dP/dx
            jacobian[i, 0] = 1.  # jacobian[i, j!=0] = 0
        elif type_c == 'T':  # dT/dx
            jacobian[i, 1] = 1.  # jacobian[i, j!=1] = 0
        elif type_c == 'S':  # dS/dx
            # dS/dP = -aV, dS/dT = Cp/T
            jacobian[i, 0:2] = [-assemblage.n_moles
                                * assemblage.alpha
                                * assemblage.molar_volume,
                                assemblage.n_moles
                                * assemblage.molar_heat_capacity_p / x[1]]
            j = 2
            for k, n in enumerate(assemblage.endmembers_per_phase):
                jacobian[i, j] = assemblage.phases[k].molar_entropy
                if n > 1:  # for solutions with >1 endmember
                    jacobian[i, j+1:j+n] = (assemblage.n_moles
                                            * assemblage.molar_fractions[k]
                                            * (assemblage.phases[k].partial_entropies[1:]
                                               - assemblage.phases[k].partial_entropies[0]))
                j += n
        elif type_c == 'V':  # dV/dx
            # dV/dP = -V/K_T, dV/dT = aV
            jacobian[i, 0:2] = [-assemblage.n_moles
                                * assemblage.molar_volume / assemblage.K_T,
                                assemblage.n_moles*assemblage.molar_volume]
            j = 2
            for k, n in enumerate(assemblage.endmembers_per_phase):
                jacobian[i, j] = assemblage.phases[k].molar_volume
                if n > 1:  # for solutions with >1 stable endmember
                    jacobian[i, j+1:j+n] = (assemblage.n_moles
                                            * assemblage.molar_fractions[k]
                                            * (assemblage.phases[k].partial_volumes[1:]
                                               - assemblage.phases[k].partial_volumes[0]))
                j += n
        elif type_c == 'PT_ellipse':
            v_scaled = (x[0:2] - eq_c[0])/eq_c[1]
            jacobian[i, 0:2] = v_scaled/(np.linalg.norm(v_scaled)*eq_c[1])
        elif type_c == 'X':
            jacobian[i, :] = eq_c[0]
        else:
            raise Exception('constraint type not recognised')

    # Next, let's get the effect of pressure and temperature
    # on each of the independent reactions
    # i.e. dF(i, reactions)/dx[0] and dF(i, reactions)/dx[1]
    partial_volumes_vector = np.zeros((assemblage.n_endmembers))
    partial_entropies_vector = np.zeros((assemblage.n_endmembers))
    j = 0
    for i, n in enumerate(assemblage.endmembers_per_phase):
        if n == 1:  # for endmembers
            partial_volumes_vector[j] = assemblage.phases[i].molar_volume
            partial_entropies_vector[j] = assemblage.phases[i].molar_entropy
        else:  # for solid solutions
            partial_volumes_vector[j:j+n] = assemblage.phases[i].partial_volumes
            partial_entropies_vector[j:j+n] = assemblage.phases[i].partial_entropies
        j += n
    reaction_volumes = np.dot(assemblage.reaction_basis, partial_volumes_vector)
    reaction_entropies = np.dot(assemblage.reaction_basis, partial_entropies_vector)

    # dGi/dP = deltaVi; dGi/dT = -deltaSi
    jacobian[2:2+len(reaction_volumes), 0] = reaction_volumes
    jacobian[2:2+len(reaction_volumes), 1] = -reaction_entropies

    # Pressure and temperature have no effect on the bulk
    # compositional constraints
    # i.e. dF(i, bulk)/dx[0] and dF(i, bulk)/dx[1] = 0

    # Finally, let's build the compositional Hessian d2G/dfidfj = dmui/dfj
    # where fj is the fraction of endmember j in a phase
    phase_amounts = np.array(assemblage.molar_fractions) * assemblage.n_moles
    comp_hessian = np.zeros((assemblage.n_endmembers, assemblage.n_endmembers))
    dfi_dxj = np.zeros((assemblage.n_endmembers, assemblage.n_endmembers))
    dpi_dxj = np.zeros((assemblage.n_endmembers, assemblage.n_endmembers))
    j = 0
    for i, n in enumerate(assemblage.endmembers_per_phase):
        # ignore pure phases, as changing the amount of a phase
        # does not change its molar gibbs free energy
        if n == 1:
            dpi_dxj[j, j] = 1.
        else:
            comp_hessian[j:j+n, j:j+n] = assemblage.phases[i].gibbs_hessian
            dfi_dxj[j:j+n, j:j+n] = np.eye(n)
            dfi_dxj[j, j:j+n] -= 1.  # x[0] = p(phase) and x[1:] = f[1:] - f[0]

            dpi_dxj[j:j+n, j:j+n] = dfi_dxj[j:j+n, j:j+n] * phase_amounts[i]
            dpi_dxj[j:j+n, j] = assemblage.phases[i].molar_fractions
        j += n

    # dfi_dxj converts the endmember hessian to the parameter hessian.
    reaction_hessian = assemblage.reaction_basis.dot(comp_hessian).dot(dfi_dxj)
    bulk_hessian = assemblage.reduced_stoichiometric_array.T.dot(dpi_dxj)
    jacobian[2:, 2:] = np.concatenate((reaction_hessian, bulk_hessian))
    return jacobian


def lambda_bounds(dx, x, endmembers_per_phase):
    """
    Returns the lambda bounds for the damped affine invariant modification
    to Newton's method for nonlinear problems (Deuflhard, 1974;1975;2004).

    Parameters
    ----------
    dx : numpy array
        The proposed newton step.

    x : numpy array
        Parameter values for the equilibrium problem to be solved.

    endmembers_per_phase : list of integers
        A list of the number of endmembers in each phase.

    Returns
    -------
    lmda_bounds : tuple of floats
        minimum and maximum allowed fractions of the full newton step (dx).
    """

    max_steps = np.ones((len(x)))*100000.

    # first two constraints are P and T
    max_steps[0:2] = [20.e9, 500.]  # biggest reasonable P and T steps

    j = 2
    for i, n in enumerate(endmembers_per_phase):
        # if the phase fraction constraint would be broken,
        # set a step that is marginally smaller
        if x[j] + dx[j] < 0.:
            max_steps[j] = max(x[j]*0.999, 0.001)
            # max_steps[j] = max(x[j]*0.5, 0.001) # allowed to cut amount by half...
        max_steps[j+1:j+n] = [max(xi*0.99, 0.01) for xi in x[j+1:j+n]]  # maximum compositional step
        j += n

    max_lmda = min([1. if step <= max_steps[i] else max_steps[i]/step
                    for i, step in enumerate(np.abs(dx))])

    return (1.e-8, max_lmda)


def phase_fraction_constraints(phase, assemblage, fractions, prm):
    """
    Converts a phase fraction constraint into standard linear form
    that can be processed by the root finding problem.

    Parameters
    ----------
    phase : burnman.SolidSolution or burnman.Mineral
        The phase for which the fraction is to be constrained

    assemblage : burnman.Composite
        The assemblage for which equilibrium is to be calculated.

    fractions : numpy array
        The phase fractions to be satified at equilibrium.

    prm : namedtuple
        A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).

    Returns
    -------
    constraints : list
        An list of the phase fraction constraints.
    """
    phase_idx = assemblage.phases.index(phase)

    constraints = []
    for fraction in fractions:
        constraints.append(['X', [np.zeros((prm.n_parameters)), 0.]])
        constraints[-1][-1][0][prm.phase_amount_indices] = -fraction
        constraints[-1][-1][0][prm.phase_amount_indices[phase_idx]] += 1.

    return constraints


def phase_composition_constraints(phase, assemblage, constraints):
    """
    Converts a phase composition constraint into standard linear form
    that can be processed by the root finding problem.

    We start with constraints in the form (site_names, n, d, v), where
    n*x/d*x = v and n and d are fixed vectors of site coefficients.
    So, one could for example choose a constraint
    ([Mg_A, Fe_A], [1., 0.], [1., 1.], [0.5]) which would
    correspond to equal amounts Mg and Fe on the A site.

    These are then converted by this function into endmember proportions
    (n'*p/d'*p = v). Because the proportions must add up to zero,
    we can reexpress this ratio into a linear constraint:
    [(n'[1:] - n'[0]) - v*(d'[1:] - d'[0])]*xi = v*d0 - n0
    which is less easy for a human to understand
    (in terms of chemical constraints), but easier to use as a constraint
    in a nonlinear solve.

    Parameters
    ----------
    phase : burnman.SolidSolution
        The phase for which the composition is to be constrained

    assemblage : burnman.Composite
        The assemblage for which equilibrium is to be calculated.

    constraints : a 4-tuple (list of strings, numpy array * 3)
        A tuple corresponding to the desired constraints, in the form
        (site_names, numerator, denominator, values).

    Returns
    -------
    x_constraints : list
        An list of the phase composition constraints in standard form.
    """
    phase_idx = assemblage.phases.index(phase)

    site_names, numerator, denominator, values = constraints
    site_indices = [phase.solution_model.site_names.index(name)
                    for name in site_names]
    noccs = phase.solution_model.endmember_noccupancies

    # Converts site constraints into endmember constraints
    # Ends up with shape (n_endmembers, 2)
    endmembers = np.dot(noccs[:, site_indices],
                        np.array([numerator, denominator]).T)

    numer0, denom0 = endmembers[0]
    endmembers -= endmembers[0]
    numer, denom = endmembers.T[:, 1:]

    # We start from the correct index
    start_idx = sum(assemblage.endmembers_per_phase[:phase_idx]) + 3
    n_indices = assemblage.endmembers_per_phase[phase_idx] - 1

    x_constraints = []
    for v in values:
        f = v*denom0 - numer0
        x_constraints.append(['X', [np.zeros((prm.n_parameters)), f]])
        x_constraints[-1][1][0][start_idx:start_idx+n_indices] = numer - v*denom

    return x_constraints


def get_equilibration_parameters(assemblage, composition):
    """
    Builds a named tuple containing the parameter names and
    various other parameters needed by the equilibrium solve.

    Parameters
    ----------
    assemblage : burnman.Composite
        The assemblage for which equilibrium is to be calculated.

    composition : dictionary
        The bulk composition for the equilibrium problem.

    Returns
    -------
    prm : namedtuple
        A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).
    """
    # Initialize a named tuple for the equilibration parameters
    prm = namedtuple('assemblage_parameters', [])

    # Process parameter names
    prm.parameter_names = ['Pressure (Pa)', 'Temperature (K)']
    for i, n in enumerate(assemblage.endmembers_per_phase):
        prm.parameter_names.append('x({0})'.format(assemblage.phases[i].name))
        if n > 1:
            p_names = ['p({0} in {1})'.format(n, assemblage.phases[i].name)
                       for n in assemblage.phases[i].endmember_names[1:]]
            prm.parameter_names.extend(p_names)

    prm.n_parameters = len(prm.parameter_names)
    prm.phase_amount_indices = [i for i in range(len(prm.parameter_names))
                                if 'x(' in prm.parameter_names[i]]

    # Find the bulk composition vector
    prm.bulk_composition_vector = np.array([composition[e]
                                            for e in assemblage.elements])

    if assemblage.compositional_null_basis.shape[0] != 0:
        if (np.abs(assemblage.compositional_null_basis.dot(prm.bulk_composition_vector)[0]) > 1.e-12):
            raise Exception('The bulk composition is not within the '
                            'compositional space of the assemblage')

    prm.reduced_composition_vector = [prm.bulk_composition_vector[i]
                                      for i in assemblage.independent_element_indices]

    prm.constraint_vector, prm.constraint_matrix = calculate_constraints(assemblage)
    return prm


def process_eq_constraints(equality_constraints, assemblage, prm):
    """
    A function that processes the equality constraints
    into a form that can be processed by the F and jacobian functions.

    This function has two main tasks: it turns phase_fraction and
    phase_composition constraints into standard linear constraints in the
    solution parameters. It also turns vector-valued constraints into a
    list of scalar-valued constraints.

    Parameters
    ----------
    equality_constraints : list
        A list of equality constraints as provided by the user.

    assemblage : burnman.Composite
        The assemblage for which equilibrium is to be calculated.

    prm : namedtuple
        A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).

    Returns
    -------
    eq_constraint_lists : list of lists
        A list of lists of equality constraints in a form that can be processed
        by the F and jacobian functions.
    """
    eq_constraint_lists = []
    for i in range(2):
        if equality_constraints[i][0] == 'phase_fraction':
            phase = equality_constraints[i][1][0]
            fraction = equality_constraints[i][1][1]
            if isinstance(fraction, float):
                fraction = np.array([fraction])
            if not isinstance(fraction, np.ndarray):
                raise Exception('The constraint fraction in equality {0} '
                                'should be a float or numpy array'.format(i+1))
            eq_constraint_lists.append(phase_fraction_constraints(phase,
                                                                  assemblage,
                                                                  fraction,
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

            eq_constraint_lists.append(phase_composition_constraints(phase,
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
                            'PT_ellipse, phase_fraction, '
                            'or phase_composition.'.format(i+1))
    return eq_constraint_lists


def equilibrate(composition, assemblage, equality_constraints,
                tol=1.e-3,
                store_iterates=False, store_assemblage=False,
                max_iterations=100., verbose=True):
    """
    A function that equilibrates an assemblage subject to given
    bulk composition and equality constraints by
    solving the equilibrium relations
    (chemical affinities for feasible reactions in the system
    should be equal to zero).

    Parameters
    ----------
    composition : dictionary
        The bulk composition that the assemblage must satisfy

    assemblage : burnman.Composite object
        The assemblage to be equilibrated

    equality_constraints : list
        A list of equality constraints.

    tol : float
        The tolerance for the nonlinear solver.

    store_iterates : boolean
        Whether to store the parameter values for each iteration in
        each solution object.

    store_assemblage : boolean
        Whether to store a copy of the assemblage object in each
        solution object.

    max_iterations : integer
        The maximum number of iterations for the nonlinear solver.

    verbose : boolean
        Whether to print output updating the user on the status of
        equilibration.

    Returns
    -------
    sol_list : single, list, or 2D list of solver solution objects

    prm : namedtuple object
        A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).
    """
    for ph in assemblage.phases:
        if isinstance(ph, SolidSolution) and not hasattr(ph, 'molar_fractions'):
            raise Exception(f'set_composition for solution {ph} before running equilibrate.')

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

            # Set the initial fractions and compositions
            # of the phases in the assemblage:
            # set_compositions_and_state_from_parameters(assemblage,
            #                                           prm.initial_parameters)
            try:
                sol = damped_newton_solve(F=lambda x: F(x, assemblage,
                                                        equality_constraints,
                                                        prm.reduced_composition_vector),
                                          J=lambda x: jacobian(x, assemblage,
                                                               equality_constraints),
                                          lambda_bounds=lambda dx, x: lambda_bounds(dx, x, assemblage.endmembers_per_phase),
                                          guess=parameters,
                                          linear_constraints=(prm.constraint_matrix,
                                                              prm.constraint_vector),
                                          tol=tol,
                                          store_iterates=store_iterates,
                                          max_iterations=max_iterations)
            except AttributeError:
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
                        dF = F(s.x, assemblage, cs,
                               prm.reduced_composition_vector)
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
                                                enumerate(new_parameters[prm.phase_amount_indices]) if v < 0.]
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
