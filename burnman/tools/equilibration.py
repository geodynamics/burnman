# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from itertools import product
from scipy.linalg import lu_factor, lu_solve
from collections import namedtuple

from ..optimize.nonlinear_solvers import damped_newton_solve
from ..classes.solution import Solution


def calculate_constraints(assemblage, n_free_compositional_vectors):
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

    :param assemblage: The assemblage for which the constraints are calculated.
    :type assemblage: :class:`burnman.Composite`

    :returns: The constraints vector and matrix.
    :rtype: tuple
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

    c_vector = np.zeros((n_constraints + 2))
    c_matrix = np.zeros(
        (n_constraints + 2, assemblage.n_endmembers + 2 + n_free_compositional_vectors)
    )  # includes P, T

    c_matrix[0, 0] = -1  # P>0
    c_matrix[1, 1] = -1  # T>0

    cidx = 2  # index of current compositional constraint
    pidx = 0  # starting index of current phase
    for i, n in enumerate(assemblage.endmembers_per_phase):
        m = len(bounds[i][0])
        # The first endmember proportion is not a free variable
        # (all endmembers in any solution must sum to one)
        # Re-express the constraints without the first endmember
        c_matrix[cidx, pidx + 2] = -1.0  # need phase proportions > 0
        cidx += 1
        if m != 0:
            c_vector[cidx : cidx + m] = -bounds[i][0]
            c_matrix[cidx : cidx + m, pidx + 1 + 2 : pidx + n + 2] = (
                np.einsum("i, j", bounds[i][0], np.ones_like(bounds[i][1:, 0]))
                - bounds[i].T[:, 1:]
            )
            cidx += m
        pidx += n

    return c_vector, c_matrix


def get_parameters(assemblage, n_free_compositional_vectors=0):
    """
    Gets the starting parameters vector (x) for the current equilibrium problem.
    These are:

      - pressure
      - temperature
      - absolute amount of each phase. if a phase is a solution
        with >1 endmember, the following parameters are the mole fractions
        of the independent endmembers in the solution, except for the first
        endmember (as the mole fractions must sum to one).

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :returns: The current values of all the parameters.
    :rtype: numpy.array
    """
    params = np.zeros(assemblage.n_endmembers + 2 + n_free_compositional_vectors)
    n_moles_phase = assemblage.n_moles * np.array(assemblage.molar_fractions)

    try:
        params[:2] = [assemblage.pressure, assemblage.temperature]
    except AttributeError:
        raise Exception("You need to set_state before getting parameters")

    j = 2
    for i, ph in enumerate(assemblage.phases):
        params[j] = n_moles_phase[i]
        if isinstance(ph, Solution):
            params[j + 1 : j + assemblage.endmembers_per_phase[i]] = assemblage.phases[
                i
            ].molar_fractions[1:]
        j += assemblage.endmembers_per_phase[i]

    return params


def get_endmember_amounts(assemblage):
    """
    Gets the absolute amounts of all the endmembers in the solution.

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :returns: The current amounts of all the endmembers.
    :rtype: numpy.array
    """
    phase_amounts = assemblage.n_moles * assemblage.molar_fractions
    amounts = np.empty(assemblage.n_endmembers)
    j = 0
    for i, ph in enumerate(assemblage.phases):
        if isinstance(ph, Solution):
            amounts[j : j + assemblage.endmembers_per_phase[i]] = (
                phase_amounts[i] * assemblage.phases[i].molar_fractions
            )
        else:
            amounts[j] = phase_amounts[i]
        j += assemblage.endmembers_per_phase[i]

    return amounts


def set_compositions_and_state_from_parameters(assemblage, parameters):
    """
    Sets the phase compositions, amounts and state of the assemblage
    from a list of parameter values.

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param parameters: The current parameter values.
    :type parameters: numpy.array
    """
    assemblage.set_state(parameters[0], parameters[1])
    i = 2
    phase_amounts = np.zeros(len(assemblage.phases))
    for phase_idx, ph in enumerate(assemblage.phases):
        phase_amounts[phase_idx] = parameters[i]
        if isinstance(ph, Solution):
            n_mbrs = len(ph.endmembers)
            f = [0.0] * n_mbrs
            f[1:] = parameters[i + 1 : i + n_mbrs]
            f[0] = 1.0 - sum(f)
            ph.set_composition(f)
            i += n_mbrs
        else:
            i += 1

    assert np.all(phase_amounts > -1.0e-8)
    phase_amounts = np.abs(phase_amounts)
    assemblage.n_moles = sum(phase_amounts)
    assemblage.set_fractions(phase_amounts / assemblage.n_moles)
    return None


def F(
    x,
    assemblage,
    equality_constraints,
    reduced_composition_vector,
    reduced_free_composition_vectors,
):
    """
    The vector-valued function for which the root is sought.
    The first two vector values depend on the
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

    :param x: Parameter values for the equilibrium problem to be solved.
    :type x: numpy array

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param equality_constraints: A list of the equality constraints
        (see above for valid formats).
    :type equality_constraints: list of lists

    :param reduced_composition_vector: The amounts of the independent
        elements.
    :type reduced_composition_vector: numpy.array

    :param reduced_free_composition_vectors: The amounts of the
        independent elements in each of the free_compositional_vectors.
    :type reduced_free_composition_vectors: 2D numpy.array

    :returns: The vector corresponding to F(x).
    :rtype: numpy.array
    """

    set_compositions_and_state_from_parameters(assemblage, x)
    new_endmember_amounts = get_endmember_amounts(assemblage)

    # We want to find the root of the following equations
    n_equality_constraints = len(equality_constraints)
    eqns = np.zeros((assemblage.n_endmembers + n_equality_constraints))
    i = 0
    for i, (type_c, eq_c) in enumerate(equality_constraints):
        if type_c == "P":
            eqns[i] = x[0] - eq_c
        elif type_c == "T":
            eqns[i] = x[1] - eq_c
        elif type_c == "S":
            eqns[i] = assemblage.molar_entropy * assemblage.n_moles - eq_c
        elif type_c == "V":
            eqns[i] = assemblage.molar_volume * assemblage.n_moles - eq_c
        elif type_c == "PT_ellipse":
            v_scaled = (x[0:2] - eq_c[0]) / eq_c[1]
            eqns[i] = np.linalg.norm(v_scaled) - 1.0
        elif type_c == "X":
            eqns[i] = np.dot(eq_c[0], x) - eq_c[1]  # i.e. Ax = b
        else:
            raise Exception("constraint type not recognised")
    i += 1
    if n_equality_constraints > 2:
        new_reduced_composition_vector = reduced_composition_vector + x[
            2 - n_equality_constraints :
        ].dot(reduced_free_composition_vectors)
    else:
        new_reduced_composition_vector = reduced_composition_vector
    eqns[i : i + assemblage.n_reactions] = assemblage.reaction_affinities
    eqns[i + assemblage.n_reactions :] = (
        np.dot(assemblage.reduced_stoichiometric_array.T, new_endmember_amounts)
        - new_reduced_composition_vector
    )
    return eqns


def jacobian(x, assemblage, equality_constraints, reduced_free_composition_vectors):
    """
    The Jacobian of the vector-valued function :math:`F` for which the
    root is sought (:math:`\\partial F / \\partial x`).
    See documentation for :func:`F` and :func:`get_parameters`
    (which return :math:`F` and :math:`x` respectively) for more details.

    :param x: Parameter values for the equilibrium problem to be solved.
    :type x: numpy.array

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param equality_constraints: A list of the equality constraints
        (see documentation for :func:`burnman.tools.equilbration.F`).
    :type equality_constraints: list of lists

    :param reduced_free_composition_vectors: The amounts of the
        independent elements in each of the free_compositional_vectors.
    :type reduced_free_composition_vectors: 2D numpy array

    :returns: The Jacobian for the equilibrium problem.
    :rtype: 2D numpy.array

    """
    # The solver always calls the Jacobian with the same
    # x parameters as used previously for the root functions
    # Therefore we don't need to set compositions or state again here.

    # First, we find out the effect of the two constraint parameters F[:2]
    # on the pressure (x[0]) and temperature (x[1]):
    n_equality_constraints = len(equality_constraints)
    jacobian = np.zeros(
        (
            assemblage.n_endmembers + n_equality_constraints,
            assemblage.n_endmembers + n_equality_constraints,
        )
    )
    ic = 0
    for ic, (type_c, eq_c) in enumerate(equality_constraints):
        if type_c == "P":  # dP/dx
            jacobian[ic, 0] = 1.0  # jacobian[i, j!=0] = 0
        elif type_c == "T":  # dT/dx
            jacobian[ic, 1] = 1.0  # jacobian[i, j!=1] = 0
        elif type_c == "S":  # dS/dx
            # dS/dP = -aV, dS/dT = Cp/T
            jacobian[ic, 0:2] = [
                -assemblage.n_moles * assemblage.alpha * assemblage.molar_volume,
                assemblage.n_moles * assemblage.molar_heat_capacity_p / x[1],
            ]
            j = 2
            for k, n in enumerate(assemblage.endmembers_per_phase):
                jacobian[ic, j] = assemblage.phases[k].molar_entropy
                if n > 1:  # for solutions with >1 endmember
                    jacobian[ic, j + 1 : j + n] = (
                        assemblage.n_moles
                        * assemblage.molar_fractions[k]
                        * (
                            assemblage.phases[k].partial_entropies[1:]
                            - assemblage.phases[k].partial_entropies[0]
                        )
                    )
                j += n
        elif type_c == "V":  # dV/dx
            # dV/dP = -V/K_T, dV/dT = aV
            jacobian[ic, 0:2] = [
                -assemblage.n_moles
                * assemblage.molar_volume
                / assemblage.isothermal_bulk_modulus_reuss,
                assemblage.n_moles * assemblage.molar_volume,
            ]
            j = 2
            for k, n in enumerate(assemblage.endmembers_per_phase):
                jacobian[ic, j] = assemblage.phases[k].molar_volume
                if n > 1:  # for solutions with >1 stable endmember
                    jacobian[ic, j + 1 : j + n] = (
                        assemblage.n_moles
                        * assemblage.molar_fractions[k]
                        * (
                            assemblage.phases[k].partial_volumes[1:]
                            - assemblage.phases[k].partial_volumes[0]
                        )
                    )
                j += n
        elif type_c == "PT_ellipse":
            v_scaled = (x[0:2] - eq_c[0]) / eq_c[1]
            jacobian[ic, 0:2] = v_scaled / (np.linalg.norm(v_scaled) * eq_c[1])
        elif type_c == "X":
            jacobian[ic, :] = eq_c[0]
        else:
            raise Exception("constraint type not recognised")
    ic += 1

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
        else:  # for solutions
            partial_volumes_vector[j : j + n] = assemblage.phases[i].partial_volumes
            partial_entropies_vector[j : j + n] = assemblage.phases[i].partial_entropies
        j += n
    reaction_volumes = np.dot(assemblage.reaction_basis, partial_volumes_vector)
    reaction_entropies = np.dot(assemblage.reaction_basis, partial_entropies_vector)

    # dGi/dP = deltaVi; dGi/dT = -deltaSi
    jacobian[ic : ic + len(reaction_volumes), 0] = reaction_volumes
    jacobian[ic : ic + len(reaction_volumes), 1] = -reaction_entropies

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
        if n == 1:
            # changing the amount (p) of a pure phase
            # does not change its fraction in that phase,
            # so dfi_dxj remains unchanged
            dpi_dxj[j, j] = 1.0
        else:
            comp_hessian[j : j + n, j : j + n] = assemblage.phases[i].gibbs_hessian

            # x[0] = p(phase) and x[1:] = f[1:] - f[0]
            # Therefore
            # df[0]/dx[0] = 0
            # df[0]/dx[1:] = -1
            # (because changing the fraction of any endmember
            # depletes the fraction of the first endmember)
            # df[1:]/dx[1:] = 1 on diagonal, 0 otherwise
            # (because all other fractions are independent of each other)
            dfi_dxj[j : j + n, j : j + n] = np.eye(n)
            dfi_dxj[j, j : j + n] -= 1.0
            # Total amounts of endmembers (p) are the fractions
            # multiplied by the amounts of their representative phases
            dpi_dxj[j : j + n, j : j + n] = (
                dfi_dxj[j : j + n, j : j + n] * phase_amounts[i]
            )
            # The derivative of the amount of each endmember with respect
            # to the amount of each phase is equal to the molar fractions
            # of the endmembers.
            dpi_dxj[j : j + n, j] = assemblage.phases[i].molar_fractions
        j += n

    # dfi_dxj converts the endmember hessian to the parameter hessian.
    reaction_hessian = assemblage.reaction_basis.dot(comp_hessian).dot(dfi_dxj)
    bulk_hessian = assemblage.reduced_stoichiometric_array.T.dot(dpi_dxj)

    if reaction_hessian.shape[0] > 0:
        jacobian[ic:, 2 : 2 + len(reaction_hessian[0])] = np.concatenate(
            (reaction_hessian, bulk_hessian)
        )
    else:
        jacobian[ic:, 2 : 2 + len(bulk_hessian[0])] = bulk_hessian

    if len(reduced_free_composition_vectors) > 0:
        jacobian[
            -reduced_free_composition_vectors.shape[1] :, 2 + len(reaction_hessian[0]) :
        ] = -reduced_free_composition_vectors.T

    return jacobian


def lambda_bounds(dx, x, endmembers_per_phase):
    """
    Returns the lambda bounds for the damped affine invariant modification
    to Newton's method for nonlinear problems (Deuflhard, 1974;1975;2004).

    :param dx: The proposed newton step.
    :type dx: numpy.array

    :param x: Parameter values for the equilibrium problem to be solved.
    :type x: numpy.array

    :param endmembers_per_phase: A list of the number of endmembers in each phase.
    :type endmembers_per_phase: list of int

    :returns: Minimum and maximum allowed fractions of the full newton step (dx).
    :rtype: tuple of floats
    """

    max_steps = np.ones((len(x))) * 100000.0

    # first two constraints are P and T
    max_steps[0:2] = [20.0e9, 500.0]  # biggest reasonable P and T steps

    j = 2
    for i, n in enumerate(endmembers_per_phase):
        # if the phase fraction constraint would be broken,
        # set a step that is marginally smaller
        if x[j] + dx[j] < 0.0:
            max_steps[j] = max(x[j] * 0.999, 0.001)
        max_steps[j + 1 : j + n] = [
            max(xi * 0.99, 0.01) for xi in x[j + 1 : j + n]
        ]  # maximum compositional step
        j += n

    max_lmda = min(
        [
            1.0 if step <= max_steps[i] else max_steps[i] / step
            for i, step in enumerate(np.abs(dx))
        ]
    )

    return (1.0e-8, max_lmda)


def phase_fraction_constraints(phase, assemblage, fractions, prm):
    """
    Converts a phase fraction constraint into standard linear form
    that can be processed by the root finding problem.

    We start with a single fraction or an array of fractions
    for a particular phase (:math:`n_p / \\sum n = f`).
    These are then converted into the "X" form of constraint
    by multiplying by :math:`\\sum n` and moving all terms to
    the LHS of the equation:

    :math:`-f n_0 - f n_1 - \\ldots + (1-f) n_p - \\ldots = 0`

    This form is less readable, but easier to use as a constraint
    in a nonlinear solve.

    :param phase: The phase for which the fraction is to be constrained
    :type phase: :class:`burnman.Solution` or :class:`burnman.Mineral`

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param fractions: The phase fractions to be satified at equilibrium.
    :type fractions: numpy.array

    :param prm: A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).
    :type prm: namedtuple

    :returns: The phase fraction constraints.
    :rtype: list
    """
    phase_idx = assemblage.phases.index(phase)

    constraints = []
    for fraction in fractions:
        constraints.append(["X", [np.zeros((prm.n_parameters)), 0.0]])
        constraints[-1][-1][0][prm.phase_amount_indices] = -fraction
        constraints[-1][-1][0][prm.phase_amount_indices[phase_idx]] += 1.0

    return constraints


def phase_composition_constraints(phase, assemblage, constraints, prm):
    """
    Converts a phase composition constraint into standard linear form
    that can be processed by the root finding problem.

    We start with constraints in the form (site_names, n, d, v), where
    :math:`(n x)/ (d x) = v` and :math:`n` and :math:`d` are fixed vectors of
    site coefficients. So, one could for example choose a constraint
    ([Mg_A, Fe_A], [1., 0.], [1., 1.], [0.5]) which would
    correspond to equal amounts Mg and Fe on the A site.

    This function converts the user-defined vectors of site constraints
    :math:`n` and :math:`d` into vectors of endmember proportion
    constraints :math:`n'` and :math:`d'`, such that
    :math:`(n' x)/ (d' x) = v`. This is done via linear transformation
    using the site occupancy matrix provided by :class:`burnman.Solution`.
    By multiplying by the denominator, we have the following scalar
    comparison: :math:`(n' x) = v (d' x)`

    The equilibration function does not use the proportion of
    the first endmember (as the endmember proportions must sum to one),
    and so we split :math:`x`, :math:`n'` and :math:`d'` into the first
    element and following elements:
    :math:`(n'_0 x_0 + n'_i x_i) = v (d'_0 x_0 + d'_i x_i)`
    where :math:`i` is taken over all elements apart from the first.

    With some more rearranging we can express the constraint in standard
    linear form:
    :math:`(n'_0 (1 - \\sum_j x_j) + n'_i x_i) = v (d'_0 (1 - \\sum_j x_j) + d'_i x_i)`

    :math:`(n'_0 + (n'_i - 1_i n'_0) x_i) = v (d'_0 + (d'_i - 1_i d'_0) x_i)`

    :math:`(((n'_i - 1_i n'_0) - v(d'_i - 1_i d'_0)) x_i) = (v d'_0 - n'_0)`

    This form is less readable, but easier to use as a constraint
    in a nonlinear solve.

    :param phase: The phase for which the composition is to be constrained.
    :type phase: :class:`burnman.Solution`

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param constraints: The desired constraints in the form:
        site_names (list of strings), numerator (numpy.array),
        denominator (numpy.array), values (numpy.array).
    :type constraints: tuple

    :returns: The phase composition constraints in standard form.
    :rtype: list
    """
    phase_idx = assemblage.phases.index(phase)

    site_names, numerator, denominator, values = constraints
    site_indices = [phase.solution_model.site_names.index(name) for name in site_names]
    noccs = phase.solution_model.endmember_noccupancies

    # Converts site constraints into endmember constraints
    # Ends up with shape (n_endmembers, 2)
    endmembers = np.dot(noccs[:, site_indices], np.array([numerator, denominator]).T)

    numer0, denom0 = endmembers[0]
    endmembers -= endmembers[0]
    numer, denom = endmembers.T[:, 1:]

    # We start from the correct index
    start_idx = sum(assemblage.endmembers_per_phase[:phase_idx]) + 3
    n_indices = assemblage.endmembers_per_phase[phase_idx] - 1

    x_constraints = []
    for v in values:
        f = v * denom0 - numer0
        x_constraints.append(["X", [np.zeros((prm.n_parameters)), f]])
        x_constraints[-1][1][0][start_idx : start_idx + n_indices] = numer - v * denom

    return x_constraints


def get_equilibration_parameters(assemblage, composition, free_compositional_vectors):
    """
    Builds a named tuple containing the parameter names and
    various other parameters needed by the equilibrium solve.

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param composition: The bulk composition for the equilibrium problem.
    :type composition: dict

    :param free_compositional_vectors: The bulk compositional
        degrees of freedom for the equilibrium problem.
    :type free_compositional_vectors: list of dictionaries

    :returns: A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).
    :rtype: namedtuple
    """
    # Initialize a named tuple for the equilibration parameters
    prm = namedtuple("assemblage_parameters", [])

    # Process parameter names
    prm.parameter_names = ["Pressure (Pa)", "Temperature (K)"]
    for i, n in enumerate(assemblage.endmembers_per_phase):
        prm.parameter_names.append("x({0})".format(assemblage.phases[i].name))
        if n > 1:
            p_names = [
                "p({0} in {1})".format(n, assemblage.phases[i].name)
                for n in assemblage.phases[i].endmember_names[1:]
            ]
            prm.parameter_names.extend(p_names)

    n_free_compositional_vectors = len(free_compositional_vectors)
    for i in range(n_free_compositional_vectors):
        prm.parameter_names.append(f"v_{i}")

    prm.n_parameters = len(prm.parameter_names)
    prm.phase_amount_indices = [
        i for i in range(len(prm.parameter_names)) if "x(" in prm.parameter_names[i]
    ]

    # Find the bulk composition vector
    prm.bulk_composition_vector = np.array(
        [composition[e] for e in assemblage.elements]
    )

    if n_free_compositional_vectors > 0:
        prm.free_compositional_vectors = np.array(
            [
                [
                    (
                        free_compositional_vectors[i][e]
                        if e in free_compositional_vectors[i]
                        else 0.0
                    )
                    for e in assemblage.elements
                ]
                for i in range(n_free_compositional_vectors)
            ]
        )
    else:
        prm.free_compositional_vectors = np.empty((0, len(assemblage.elements)))

    if assemblage.compositional_null_basis.shape[0] != 0:
        if (
            np.abs(
                assemblage.compositional_null_basis.dot(prm.bulk_composition_vector)[0]
            )
            > 1.0e-12
        ):
            raise Exception(
                "The bulk composition is not within the "
                "compositional space of the assemblage"
            )

    prm.reduced_composition_vector = prm.bulk_composition_vector[
        assemblage.independent_element_indices
    ]
    prm.reduced_free_composition_vectors = prm.free_compositional_vectors[
        :, assemblage.independent_element_indices
    ]
    prm.constraint_vector, prm.constraint_matrix = calculate_constraints(
        assemblage, n_free_compositional_vectors
    )
    return prm


def process_eq_constraints(equality_constraints, assemblage, prm):
    """
    A function that processes the equality constraints
    into a form that can be processed by the F and jacobian functions.

    This function has two main tasks: it turns phase_fraction and
    phase_composition constraints into standard linear constraints in the
    solution parameters. It also turns vector-valued constraints into a
    list of scalar-valued constraints.

    :param equality_constraints: A list of equality constraints.
        For valid types of constraints, please see the documentation for
        :func:`burnman.equilibrate`.
    :type equality_constraints: list

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param prm: A tuple with attributes n_parameters
        (the number of parameters for the current equilibrium problem)
        and phase_amount_indices (the indices of the parameters that
        correspond to phase amounts).
    :type prm: namedtuple

    :returns: Equality constraints in a form that can be processed
        by the F and jacobian functions.
    :rtype: list of lists
    """
    eq_constraint_lists = []
    for i in range(len(equality_constraints)):
        if equality_constraints[i][0] == "phase_fraction":
            phase = equality_constraints[i][1][0]
            fraction = equality_constraints[i][1][1]
            if isinstance(fraction, float):
                fraction = np.array([fraction])
            if not isinstance(fraction, np.ndarray):
                raise Exception(
                    "The constraint fraction in equality {0} "
                    "should be a float or numpy array".format(i + 1)
                )
            eq_constraint_lists.append(
                phase_fraction_constraints(phase, assemblage, fraction, prm)
            )
        elif equality_constraints[i][0] == "phase_composition":
            phase = equality_constraints[i][1][0]
            constraint = equality_constraints[i][1][1]
            if isinstance(constraint[3], float):
                constraint = (
                    constraint[0],
                    constraint[1],
                    constraint[2],
                    np.array([constraint[3]]),
                )
            if not isinstance(constraint[3], np.ndarray):
                raise Exception(
                    "The last constraint parameter in equality "
                    "{0} should be a float "
                    "or numpy array".format(i + 1)
                )

            eq_constraint_lists.append(
                phase_composition_constraints(phase, assemblage, constraint, prm)
            )
        elif equality_constraints[i][0] == "X":
            constraint = equality_constraints[i][1]
            if isinstance(constraint[-1], float):
                constraint = (constraint[0], np.array([constraint[-1]]))
            if not isinstance(constraint[-1], np.ndarray):
                raise Exception(
                    "The last constraint parameter in "
                    "equality {0} should be "
                    "a float or numpy array".format(i + 1)
                )

            eq_constraint_lists.append(
                [["X", [constraint[0], p]] for p in constraint[1]]
            )

        elif (
            equality_constraints[i][0] == "P"
            or equality_constraints[i][0] == "T"
            or equality_constraints[i][0] == "PT_ellipse"
            or equality_constraints[i][0] == "S"
            or equality_constraints[i][0] == "V"
        ):
            if isinstance(equality_constraints[i][1], float):
                equality_constraints[i] = (
                    equality_constraints[i][0],
                    np.array([equality_constraints[i][1]]),
                )
            if not isinstance(equality_constraints[i][1], np.ndarray):
                raise Exception(
                    "The last parameter in "
                    f"equality_constraint[{i+1}] should be a "
                    "float or numpy array"
                )
            eq_constraint_lists.append(
                [[equality_constraints[i][0], p] for p in equality_constraints[i][1]]
            )
        else:
            raise Exception(
                "The type of equality_constraint is "
                "not recognised for constraint {0}.\n"
                "Should be one of P, T, S, V, X,\n"
                "PT_ellipse, phase_fraction, "
                "or phase_composition.".format(i + 1)
            )
    return eq_constraint_lists


def equilibrate(
    composition,
    assemblage,
    equality_constraints,
    free_compositional_vectors=[],
    tol=1.0e-3,
    store_iterates=False,
    store_assemblage=True,
    max_iterations=100.0,
    verbose=False,
):
    """
    A function that finds the thermodynamic equilibrium state of an
    assemblage subject to given equality constraints by solving a set of
    nonlinear equations related to the chemical potentials and
    other state variables of the system.

    The user chooses an assemblage (e.g. olivine, garnet and orthopyroxene)
    and :math:`2+n_c` equality constraints, where :math:`n_c` is the number of
    bulk compositional degrees of freedom. The equilibrate function attempts to
    find the remaining unknowns that satisfy those constraints.

    There are a number of equality constraints implemented in burnman. These are
    given as a list of lists. Each constraint should have the form:
    [<constraint type>, <constraint>], where
    <constraint type> is one of 'P', 'T', 'S', 'V', 'X', 'PT_ellipse',
    'phase_fraction', or 'phase_composition'. The format of the
    <constraint> object depends on the constraint type:
        - P: float or numpy.array of
            pressures [Pa]
        - T: float or numpy.array of
            temperatures [K]
        - S: float or numpy.array of
            entropies [J/K]
        - V: float or numpy.array of
            volumes [m:math:`^3`]
        - PT_ellipse: list of two floats or numpy.arrays, where the equality
            satifies the equation norm(([P, T] - arr[0])/arr[1]) = 1
        - phase_fraction: tuple with the form (<phase> <fraction(s)>),
            where <phase> is one of the phase objects in the assemblage
            and <fraction(s)> is a float or numpy.array corresponding
            to the desired phase fractions.
        - phase_composition: tuple with the form (<phase> <constraint>),
            where <phase> is one of the phase objects in the assemblage
            and <constraint> has the form (site_names, n, d, v),
            where :math:`(nx)/(dx) = v`, n and d are constant vectors of
            site coefficients, and v is a float or numpy.array. For example,
            a constraint of the form ([Mg_A, Fe_A], [1., 0.], [1., 1.], [0.5])
            would correspond to equal amounts Mg and Fe on the A site
            (Mg_A / (Mg_A + Fe_A) = 0.5).
        - X: list of a numpy.array and a float or numpy.array,
            where the equality satisfies the linear equation arr[0] x = arr[1].
            The first numpy.array is fixed, and the second is to be looped over
            by the equilibrate function.
            This is a generic compositional equality constraint.

    :param composition: The bulk composition that the assemblage must satisfy.
    :type composition: dict

    :param assemblage: The assemblage to be equilibrated.
    :type assemblage: :class:`burnman.Composite`

    :param equality_constraints: The list of equality constraints. See above
        for valid formats.
    :type equality_constraints: list of list

    :param free_compositional_vectors: A list of dictionaries containing
        the compositional freedom of the solution. For example, if the list
        contains the vector {'Mg': 1., 'Fe': -1}, that implies that the bulk
        composition is equal to composition + :math:`a` (n_Mg - n_Fe),
        where the value of :math:`a` is to be determined by the solve.
        Vector given in atomic (molar) units of elements.
    :type free_compositional_vectors: list of dict

    :param tol: The tolerance for the nonlinear solver.
    :type tol: float

    :param store_iterates: Whether to store the parameter values for
        each iteration in each solution object.
    :type store_iterates: bool

    :param store_assemblage: Whether to store a copy of the assemblage
        object in each solution object.
    :type store_assemblage: bool

    :param max_iterations: The maximum number of iterations for the
        nonlinear solver.
    :type max_iterations: int

    :param verbose: Whether to print output updating the user on the status of
        equilibration.
    :type verbose: bool

    :returns: Solver solution object (or a list, or a 2D list of solution objects)
        created by :func:`burnman.optimize.nonlinear_solvers.damped_newton_solve`,
        and a namedtuple object created by
        :func:`burnman.tools.equilibration.get_equilibration_parameters`. See
        documentation of these functions for the return types. If store_assemblage
        is True, each solution object also has an attribute called `assemblage`,
        which contains a copy of the input assemblage with the equilibrated
        properties. So, for a 2D grid of solution objects, one could call either
        sols[0][1].x[0] or sols[0][1].assemblage.pressure to get the pressure.
    :rtype: tuple
    """
    for ph in assemblage.phases:
        if isinstance(ph, Solution) and not hasattr(ph, "molar_fractions"):
            raise Exception(
                f"set_composition for solution {ph} before running equilibrate."
            )

    if assemblage.molar_fractions is None:
        n_phases = len(assemblage.phases)
        f = 1.0 / float(n_phases)
        assemblage.set_fractions([f for i in range(n_phases)])

    assemblage.n_moles = sum(composition.values()) / sum(assemblage.formula.values())

    n_equality_constraints = len(equality_constraints)
    n_free_compositional_vectors = len(free_compositional_vectors)

    if n_equality_constraints != n_free_compositional_vectors + 2:
        raise Exception(
            "The number of equality constraints "
            f"(currently {n_equality_constraints}) "
            "must be two more than the number of "
            "free_compositional vectors "
            f"(currently {n_free_compositional_vectors})."
        )

    for v in free_compositional_vectors:
        if np.abs(sum(v.values())) > 1.0e-12:
            raise Exception(
                "The amounts of each free_compositional_vector" "must sum to zero"
            )

    # Make parameter tuple
    prm = get_equilibration_parameters(
        assemblage, composition, free_compositional_vectors
    )

    # Check equality constraints have the correct structure
    # Convert into the format readable by the function and jacobian functions
    eq_constraint_lists = process_eq_constraints(equality_constraints, assemblage, prm)

    # Set up solves
    nc = [len(eq_constraint_list) for eq_constraint_list in eq_constraint_lists]

    # Find the initial state (could be none here)
    initial_state = [assemblage.pressure, assemblage.temperature]

    # Reset initial state if equality constraints
    # are related to pressure or temperature
    for i in range(n_equality_constraints):
        if eq_constraint_lists[i][0][0] == "P":
            initial_state[0] = eq_constraint_lists[i][0][1]
        elif eq_constraint_lists[i][0][0] == "T":
            initial_state[1] = eq_constraint_lists[i][0][1]
        elif eq_constraint_lists[i][0][0] == "PT_ellipse":
            initial_state = eq_constraint_lists[i][0][1][1]

    if initial_state[0] is None:
        initial_state[0] = 5.0e9
    if initial_state[1] is None:
        initial_state[1] = 1200.0

    assemblage.set_state(*initial_state)
    parameters = get_parameters(assemblage, n_free_compositional_vectors)

    # Solve the system of equations, loop over input parameters
    sol_array = np.empty(shape=tuple(nc), dtype="object")

    # Loop over problems
    problems = list(product(*[list(range(nc[i])) for i in range(len(nc))]))
    n_problems = len(problems)
    for i_problem, i_c in enumerate(problems):
        if verbose:
            string = "Processing solution"
            for i in range(len(i_c)):
                string += " {0}/{1}".format(i_c[i] + 1, nc[i])

            print(string + ":")

        equality_constraints = [eq_constraint_lists[i][i_c[i]] for i in range(len(nc))]

        # Set the initial fractions and compositions
        # of the phases in the assemblage:
        sol = damped_newton_solve(
            F=lambda x: F(
                x,
                assemblage,
                equality_constraints,
                prm.reduced_composition_vector,
                prm.reduced_free_composition_vectors,
            ),
            J=lambda x: jacobian(
                x,
                assemblage,
                equality_constraints,
                prm.reduced_free_composition_vectors,
            ),
            lambda_bounds=lambda dx, x: lambda_bounds(
                dx, x, assemblage.endmembers_per_phase
            ),
            guess=parameters,
            linear_constraints=(prm.constraint_matrix, prm.constraint_vector),
            tol=tol,
            store_iterates=store_iterates,
            max_iterations=max_iterations,
        )

        if sol.success and len(assemblage.reaction_affinities) > 0.0:
            maxres = np.max(np.abs(assemblage.reaction_affinities)) + 1.0e-5
            assemblage.equilibrium_tolerance = maxres

        if store_assemblage:
            sol.assemblage = assemblage.copy()
            if sol.success and len(assemblage.reaction_affinities) > 0.0:
                sol.assemblage.equilibrium_tolerance = maxres

        if verbose:
            print(sol.text)

        sol_array[i_c] = sol

        # Next, we use the solution values and Jacobian
        # to provide a starting guess for the next problem.
        # First, we find the equality constraints for the next problem
        if i_problem < n_problems - 1:
            next_i_c = problems[i_problem + 1]

            next_equality_constraints = [
                eq_constraint_lists[i][next_i_c[i]] for i in range(len(nc))
            ]

            # We use the nearest solutions as potential starting points
            # to make the next guess
            prev_sols = []
            for i in range(len(nc)):
                if next_i_c[i] != 0:
                    prev_i_c = np.copy(next_i_c)
                    prev_i_c[i] -= 1
                    prev_sols.append(sol_array[tuple(prev_i_c)])

            updated_params = False
            for s in prev_sols:
                if s.success and not updated_params:
                    # next guess based on a Newton step
                    # using the old solution vector and Jacobian
                    # with the new constraints.
                    dF = F(
                        s.x,
                        assemblage,
                        next_equality_constraints,
                        prm.reduced_composition_vector,
                        prm.reduced_free_composition_vectors,
                    )
                    luJ = lu_factor(s.J)
                    new_parameters = s.x + lu_solve(luJ, -dF)
                    c = (
                        prm.constraint_matrix.dot(new_parameters)
                        + prm.constraint_vector
                    )
                    if all(c <= 0.0):  # accept new guess
                        parameters = new_parameters
                    else:  # use the parameters from this step
                        parameters = s.x
                        exhausted_phases = [
                            assemblage.phases[phase_idx].name
                            for phase_idx, v in enumerate(
                                new_parameters[prm.phase_amount_indices]
                            )
                            if v < 0.0
                        ]
                        if len(exhausted_phases) > 0 and verbose:
                            print(
                                "A phase might be exhausted before the "
                                f"next step: {exhausted_phases}"
                            )

                    updated_params = True

    # Finally, make dimensions of sol_array equal the input dimensions
    if np.prod(sol_array.shape) > 1:
        sol_array = np.squeeze(sol_array)
    else:
        sol_array = sol_array.flatten()[0]

    return sol_array, prm
