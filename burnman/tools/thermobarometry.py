import warnings
import numpy as np
import sympy as sp

from scipy.linalg import sqrtm, block_diag
from scipy.optimize import minimize

from ..classes.solution import Solution
from ..classes.combinedmineral import CombinedMineral

from .polytope import (
    greedy_independent_endmember_selection,
    solution_polytope_from_endmember_occupancies,
)


def get_reaction_matrix(assemblage, small_fraction_tol=0.0):
    """
    Set up a matrix of independent reactions for an assemblage,
    possibly reducing the number of endmembers by excluding components
    with small molar fractions.

    :param assemblage: The mineral assemblage for which to build the
        reaction matrix.
    :type assemblage: Assemblage

    :param small_fraction_tol: If > 0.0, reduces the number of endmembers
        in solution phases by transforming to a smaller set of independent
        endmembers using a greedy algorithm and excluding those with molar
        fractions smaller than this value during the inversion.
    :type small_fraction_tol: float, optional, default 0.0

    :return: The reaction matrix for the assemblage,
        returned in the original endmember basis.
    :rtype: 2D np.array
    """
    transformation_matrix = []
    if small_fraction_tol > 0.0:
        transformation_matrices = []
        for phase in assemblage.phases:
            if isinstance(phase, Solution):
                occs = phase.solution_model.endmember_occupancies
                poly = solution_polytope_from_endmember_occupancies(occs)
                S = poly.endmember_occupancies
                T = poly.endmembers_as_independent_endmember_amounts
                f = phase.molar_fractions
                v = f.dot(occs)
                mbrs = greedy_independent_endmember_selection(
                    S, v, small_fraction_tol=small_fraction_tol, norm_tol=1e-12
                )
                idx, _, _ = mbrs
                transformation_matrices.append(T[idx, :])
            else:
                transformation_matrices.append(np.array([[1.0]]))

        transformation_matrix = block_diag(*transformation_matrices)
    else:
        transformation_matrix = np.identity(len(assemblage.endmember_names))

    # Here we build the reaction matrix R in the reduced endmember basis
    # i.e. after applying the transformation to the stoichiometric matrix
    S = sp.Matrix(transformation_matrix.dot(assemblage.stoichiometric_array))
    S = S.applyfunc(lambda x: sp.nsimplify(x))
    R = np.array([v[:] for v in S.T.nullspace()], dtype=float)
    if len(R) == 0:
        raise Exception("No independent reactions found for the given assemblage!")
    R = R.dot(transformation_matrix)  # transform back to original basis
    return R


def _build_endmember_covariance_matrix(assemblage, dataset_covariances):
    """
    Build the Gibbs free energy covariance matrix for the given assemblage.

    :param assemblage: The mineral assemblage for which to build the covariance matrix.
    :type assemblage: Assemblage

    :param dataset_covariances: The covariance data from the thermodynamic dataset.
    :type dataset_covariances: dict, with keys 'endmember_names' and 'covariance_matrix'
    """
    # For endmembers that are made up of other endmembers (the "made" endmembers in HP-speak)
    # we need to build a transformation matrix from the "raw" endmembers (those in the covariance matrix)
    # to the "transformed" endmembers (those in our assemblage)
    cov = dataset_covariances

    A = []  # transformation matrix to go from raw to transformed endmembers

    # Loop over phases in the assemblage to build the matrix A
    for phase in assemblage.phases:
        if isinstance(phase, Solution):
            for endmember, _ in phase.endmembers:
                if isinstance(endmember, CombinedMineral):
                    indices = [
                        cov["endmember_names"].index(name)
                        for name in endmember.mixture.endmember_names
                    ]
                    b = endmember.mixture.molar_fractions
                    A.append([0] * len(cov["endmember_names"]))
                    for i, idx in enumerate(indices):
                        A[-1][idx] = b[i]
                else:
                    A.append([0] * len(cov["endmember_names"]))
                    A[-1][cov["endmember_names"].index(endmember.name)] = 1.0
        else:  # single endmember phase
            A.append([0] * len(cov["endmember_names"]))
            A[-1][cov["endmember_names"].index(phase.name)] = 1.0

    # Convert A to a numpy array
    A = np.array(A)

    # Transform the covariance matrix using matrix A
    cov_transformed = A.dot(cov["covariance_matrix"]).dot(A.T)
    eigenvalues, _ = np.linalg.eig(cov_transformed)
    assert np.all(
        eigenvalues > -1.0e-5
    ), f"Transformed covariance matrix has large negative eigenvalues: {min(eigenvalues)}!"

    return cov_transformed


def _build_solution_covariance_matrix(assemblage):
    """
    Build the Gibbs free energy covariance matrix contribution from
    the compositional uncertainties of solution phases in the assemblage.

    :param assemblage: The mineral assemblage for which to build the
        covariance matrix.
    :type assemblage: Assemblage

    :return: The Gibbs free energy covariance matrix contribution from
        the compositional uncertainties of solution phases in the assemblage.
    :rtype: 2D np.array
    """
    n_mbrs_all = len(assemblage.endmember_names)
    cov_RTlna = np.zeros((n_mbrs_all, n_mbrs_all))

    j = 0
    for phase in assemblage.phases:
        if isinstance(phase, Solution):
            H_G = phase.gibbs_hessian
            n_mbrs_phase = H_G.shape[0]
            cov_X_phase = np.array(phase.compositional_covariances)

            if cov_X_phase.ndim != 2:
                raise Exception(
                    "Compositional uncertainties for each phase must be "
                    "provided as a 2D covariance matrix."
                )
            if cov_X_phase.shape[0] != cov_X_phase.shape[1]:
                raise Exception(
                    "Compositional uncertainties covariance matrix must be square."
                )
            if cov_X_phase.shape[0] != n_mbrs_phase:
                raise Exception(
                    "Compositional uncertainties covariance matrix "
                    "must have the same number of rows as phases."
                )

            cov_RTlna_phase = H_G.dot(cov_X_phase).dot(H_G.T)
            cov_RTlna[j : j + n_mbrs_phase, j : j + n_mbrs_phase] = cov_RTlna_phase
            j += n_mbrs_phase
        else:  # Add formula if it is an endmember
            cov_RTlna[j, j] = 1.0e-12  # small uncertainty for pure phases
            j += 1
    return cov_RTlna


def assemblage_set_state_from_params(assemblage, params):
    """
    Set the state of the assemblage (P, T, compositions) from the given
    list of parameters.

    :param assemblage: The mineral assemblage for which to set the state.
        If there are free compositional vectors for any solution phases
        (e.g., to account for an unknown state of order at fixed composition),
        the corresponding phases must have the attributes
        `baseline_composition` and `free_compositional_vectors` set,
        and the length of params must match the number of
        free compositional vectors plus two (for P and T).
    :type assemblage: Assemblage

    :param params: List of parameters, where the first two are pressure (Pa)
        and temperature (K), and any additional parameters correspond to
        compositional degrees of freedom. Each compositional degree of freedom
        corresponds to a free compositional vector in one of the
        solution phases in the assemblage.
    :type params: list or np.array

    :return: None
    :rtype: None
    """
    # If there are compositional degrees of freedom, set the compositions
    # according to the values in params
    if len(params) > 2:
        vi = 2
        for phase in assemblage.phases:
            if hasattr(phase, "free_compositional_vectors"):
                nv_in_phase = phase.free_compositional_vectors.shape[0]
                v = np.array(params[vi : vi + nv_in_phase])
                dc = phase.free_compositional_vectors.T.dot(v)
                new_composition = phase.baseline_composition + dc
                phase.set_composition(new_composition)
                vi += nv_in_phase

    # Set pressure and temperature
    assemblage.set_state(params[0], params[1])
    return None


def assemblage_affinity_covariance_matrix(
    assemblage,
    reaction_matrix,
    dataset_covariances=None,
    include_state_uncertainties=False,
):
    """
    Compute the covariance matrix of the affinities of the independent
    reactions for the given assemblage. This can include contributions from
    compositional uncertainties, dataset uncertainties, and from uncertainties
    in pressure and temperature.

    :param assemblage: The mineral assemblage for which to compute
        the covariance matrix. Should have its state (P, T, compositions)
        set prior to calling this function.
    :type assemblage: Assemblage

    :param reaction_matrix: The reaction matrix for the assemblage.
    :type reaction_matrix: 2D np.array

    :param dataset_covariances: The covariance data from the thermodynamic dataset.
    :type dataset_covariances: dict, with keys 'endmember_names' and 'covariance_matrix'.
        Default is None, in which case only compositional uncertainties are considered.

    :param include_state_uncertainties: If True, includes the contribution
        from uncertainties in pressure and temperature to the covariance matrix.
        If True, the assemblage must have the attribute `state_covariances`.
    :type include_state_uncertainties: bool

    :return: Covariance matrix of the affinities of the independent reactions.
    :rtype: 2D np.array
    """
    R = reaction_matrix
    cov_mu = _build_solution_covariance_matrix(assemblage)
    if dataset_covariances is not None:
        cov_mu += _build_endmember_covariance_matrix(assemblage, dataset_covariances)
    if include_state_uncertainties:
        jacobian = np.array(
            [
                assemblage.endmember_partial_volumes,
                -assemblage.endmember_partial_entropies,
            ]
        ).T
        cov_mu += jacobian.dot(assemblage.state_covariances).dot(jacobian.T)
    acov = R.dot(cov_mu).dot(R.T)
    return acov


def assemblage_affinity_misfit(
    assemblage,
    reaction_matrix,
    dataset_covariances=None,
    include_state_uncertainties=False,
):
    """
    Compute the objective misfit function (chi-squared) for given P and T.

    :param assemblage: The mineral assemblage for which to compute the misfit.
        Should have its state (P, T, compositions) set prior to calling this function.
    :type assemblage: Assemblage

    :param reaction_matrix: The reaction matrix for the assemblage.
    :type reaction_matrix: 2D np.array

    :param dataset_covariances: The covariance data from the thermodynamic dataset.
    :type dataset_covariances: dict, with keys 'endmember_names' and 'covariance_matrix'.
        Default is None, in which case only compositional uncertainties are considered.

    :param include_state_uncertainties: If True, includes the contribution
        from uncertainties in pressure and temperature to the covariance matrix.
        If True, the assemblage must have the attribute `state_covariances`.
    :type include_state_uncertainties: bool

    :return: Chi-squared misfit value.
    :rtype: float
    """
    R = reaction_matrix
    cov_a = assemblage_affinity_covariance_matrix(
        assemblage, R, dataset_covariances, include_state_uncertainties
    )
    a = R.dot(assemblage.endmember_partial_gibbs)
    return a.T.dot(np.linalg.pinv(cov_a)).dot(a)


def assemblage_state_misfit(assemblage):
    """
    Compute the objective misfit function (chi-squared) for given P and T
    based on prior expectations for P and T and their covariance.
    The assemblage must have attributes `state_priors` and
    `state_inverse_covariances`.

    :param assemblage: The mineral assemblage for which to compute the misfit.
    :type assemblage: Assemblage

    :return: Chi-squared misfit value.
    :rtype: float
    """
    delta_conditions = np.array(
        [
            assemblage.pressure - assemblage.state_priors[0],
            assemblage.temperature - assemblage.state_priors[1],
        ]
    )
    invcov = assemblage.state_inverse_covariances
    return delta_conditions.T.dot(invcov).dot(delta_conditions)


def estimate_conditions(
    assemblage,
    dataset_covariances=None,
    include_state_misfit=False,
    guessed_conditions=[1.0e9, 1000.0],
    pressure_bounds=[1.0e5, 400.0e9],
    temperature_bounds=[300.0, 4000.0],
    P_scaling=1.0e6,
    small_fraction_tol=0.0,
    max_it=100,
):
    """
    Perform a least-squares inversion to find the optimal pressure and
    temperature for a given mineral assemblage.
    Algorithm modified from Powell and Holland (1994).

    :param assemblage: The mineral assemblage for which to perform the
        inversion. Each solution phase in the assemblage must have its
        composition set along with its compositional covariance matrix
        (called `compositional_covariances`).
        If there are compositional degrees of freedom, they can be added by
        setting the attribute `free_compositional_vectors` on the relevant
        solution phases, where each row of the array corresponds to a
        free compositional vector, and the columns correspond to the
        amounts of endmembers of that phase in each vector.
    :type assemblage: Assemblage

    :param dataset_covariances: The covariance data from the thermodynamic
        dataset.
    :type dataset_covariances: dict, with keys 'endmember_names' and
        'covariance_matrix'. Default is None, in which case only
        compositional uncertainties are considered.

    :param include_state_misfit: If True, includes the misfit from
        prior expectations on P and T. The assemblage must also have
        attributes `state_priors` and `state_inverse_covariances`.
    :type include_state_misfit: bool

    :param guessed_conditions: Initial guess for pressure (Pa) and
        temperature (K). If not provided, the initial guess will be taken
        from the current state of the assemblage.
    :type guessed_conditions: np.array of shape (2,), optional, default None

    :param pressure_bounds: Bounds for pressure (Pa) during optimization.
    :type pressure_bounds: list of two floats

    :param temperature_bounds: Bounds for temperature (K) during optimization.
    :type temperature_bounds: list of two floats

    :param P_scaling: Scaling factor for pressure to improve numerical stability.
    :type P_scaling: float

    :param small_fraction_tol: If > 0.0, reduces the number of endmembers in
        solution phases by transforming to a smaller set of independent
        endmembers using a greedy algorithm and excluding those with
        molar fractions smaller than this value during the inversion.
    :type small_fraction_tol: float, optional, default 0.0

    :param max_it: Maximum number of iterations for the optimization algorithm.
    :type max_it: int, optional, default 100

    :return: Result object from the optimization containing the optimal
        conditions (x, which includes P, T, and any free compositional vectors)
        and other properties of the solution. These include
        the covariance matrix of the estimated parameters (xcov),
        the standard deviations of the estimated parameters (var),
        the correlation matrix (xcorr),
        the affinities of the independent reactions at the optimal conditions (affinities),
        the covariance matrix of the affinities (acov),
        the affinities weighted by the inverse square root of their covariance matrix (weighted_affinities),
        the partial derivatives of the affinities with respect to P and T (dadx),
        the number of independent reactions (n_reactions), the number of fitted parameters (n_params),
        the degrees of freedom (degrees_of_freedom), the chi-squared value (fun),
        the reduced chi-squared value (reduced_chisqr,
        given by the reduced chi-squared divided by the number of degrees of freedom),
        and the fit quality (fit, given by the square root of the reduced chi-squared value).
    :rtype: OptimizeResult
    """

    if include_state_misfit:
        assert hasattr(assemblage, "state_priors") and hasattr(
            assemblage, "state_inverse_covariances"
        ), (
            "To include state misfit, the assemblage must have "
            "attributes 'state_priors' and 'state_inverse_covariances'."
        )

    # Set initial guess for P and T if guessed_conditions is not provided
    if guessed_conditions is None:
        try:
            guessed_conditions = np.array([assemblage.pressure, assemblage.temperature])
        except Exception as e:
            raise ValueError(
                "guessed_conditions was not passed as an "
                "argument to the function, and could not "
                "be inferred from the current state of "
                "the assemblage."
            ) from e

    # Count the number of free compositional vectors across all phases
    n_free_vectors = 0
    for phase in assemblage.phases:
        if hasattr(phase, "free_compositional_vectors"):
            n_free_vectors += phase.free_compositional_vectors.shape[0]
            phase.baseline_composition = phase.molar_fractions.copy()

    # Check if P and/or T are fixed
    P_fixed = np.isclose(pressure_bounds[0], pressure_bounds[1])
    T_fixed = np.isclose(temperature_bounds[0], temperature_bounds[1])
    if P_fixed and T_fixed and n_free_vectors == 0:
        raise Exception(
            "Both pressure and temperature cannot be fixed if there "
            "are no free compositional vectors!"
        )

    # Set up the reaction matrix for the assemblage and assign it to an
    # assemblage attribute called reaction_matrix_for_optimization
    R = get_reaction_matrix(assemblage, small_fraction_tol)

    # Check that there are enough constraints on the problem
    n_params = 2 + n_free_vectors - int(P_fixed) - int(T_fixed)
    n_constraints = 2 if include_state_misfit else 0
    n_reactions = R.shape[0]
    if n_reactions + n_constraints < n_params:
        raise Exception(
            f"Not enough independent reactions ({n_reactions}) to constrain "
            "the inversion! You may need to relax small_fraction_tol, "
            "fix P or T through the pressure_bounds or temperature_bounds, "
            "or add priors on P and T via condition_covariances "
            "and condition_priors."
        )

    # Define the chi-squared function to minimize
    def chisqr(args):
        assemblage_set_state_from_params(assemblage, [args[0] * P_scaling, *args[1:]])
        # Compute the misfit. We do not include the state uncertainties
        # in the affinity misfit because we are optimizing over P and T directly
        # i.e., we are finding the misfit assuming that P and T are known exactly
        chisqr = assemblage_affinity_misfit(
            assemblage, R, dataset_covariances, include_state_uncertainties=False
        )

        # We do, however, include the state misfit if priors are provided
        # because this adds additional constraints to the problem
        if include_state_misfit:
            chisqr += assemblage_state_misfit(assemblage)
        return chisqr

    # Set up initial guess, bounds, and options for the optimizer
    x0 = list(guessed_conditions)
    x0[0] /= P_scaling
    x0.extend([0.0] * n_free_vectors)
    bounds = [
        (pressure_bounds[0] / P_scaling, pressure_bounds[1] / P_scaling),
        (temperature_bounds[0], temperature_bounds[1]),
        *[(None, None)] * n_free_vectors,
    ]
    options = {"maxiter": max_it}

    # Solve first with Nelder-Mead if there are free compositional vectors
    # Nelder-Mead is more robust for problems where the feasible region
    # may be complex due to compositional degrees of freedom
    if n_free_vectors > 0:
        res = minimize(chisqr, x0, method="Nelder-Mead", bounds=bounds, options=options)
        x0 = res.x

    # Solve with SLSQP, which returns the Jacobian needed for uncertainty estimation
    res = minimize(chisqr, x0, method="SLSQP", bounds=bounds, options=options)

    # Rescale solution and Jacobian back to SI units (Pa for pressure)
    res.x[0] *= P_scaling
    res.jac[0] /= P_scaling

    # Post-process results
    # First we reset the assemblage to the optimal P and T
    assemblage_set_state_from_params(assemblage, res.x)

    # Compute the covariance matrix of the estimated parameters
    # and other statistics
    res.acov = assemblage_affinity_covariance_matrix(assemblage, R, dataset_covariances)
    dmudP = assemblage.endmember_partial_volumes
    dmudT = -assemblage.endmember_partial_entropies
    res.dadx = R.dot(np.array([dmudP, dmudT]).T)
    res.xcov = np.linalg.pinv(
        res.dadx.T.dot(np.linalg.pinv(res.acov)).dot(res.dadx)
        + (assemblage.state_inverse_covariances if include_state_misfit else 0)
    )
    res.var = np.sqrt(np.diag(res.xcov))
    res.xcorr = res.xcov / np.outer(res.var, res.var)

    res.affinities = R.dot(assemblage.endmember_partial_gibbs)
    res.weighted_affinities = np.linalg.pinv(sqrtm(res.acov)).dot(res.affinities)

    res.n_reactions = n_reactions
    res.n_constraints = n_constraints
    res.n_params = n_params
    res.degrees_of_freedom = res.n_reactions + res.n_constraints - res.n_params
    if res.degrees_of_freedom > 0:
        res.reduced_chisqr = res.fun / res.degrees_of_freedom
        res.fit = np.sqrt(res.reduced_chisqr)
    else:
        warnings.warn(
            "Degrees of freedom <= 0, cannot compute reduced chi-squared and fit quality.",
            RuntimeWarning,
        )

    return res
