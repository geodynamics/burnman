import numpy as np
from scipy.linalg import sqrtm, block_diag
from scipy.optimize import minimize
import sympy as sp

from ..classes.solution import Solution
from ..classes.combinedmineral import CombinedMineral

from .polytope import (
    greedy_independent_endmember_selection,
    solution_polytope_from_endmember_occupancies,
)


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
    n_mbrs_all = len(assemblage.endmember_names)
    cov_RTlna = np.zeros((n_mbrs_all, n_mbrs_all))

    j = 0
    for phase in assemblage.phases:
        if isinstance(phase, Solution):
            H_G = phase.gibbs_hessian
            n_mbrs_phase = H_G.shape[0]
            cov_X_phase = np.array(phase.compositional_covariances)

            assert (
                cov_X_phase.ndim == 2
            ), "Compositional uncertainties must be provided as a covariance matrix."
            assert (
                cov_X_phase.shape[0] == n_mbrs_phase
            ), "Compositional uncertainties covariance matrix must have the same number of rows as phases."
            assert (
                cov_X_phase.shape[1] == n_mbrs_phase
            ), "Compositional uncertainties covariance matrix must be square."

            cov_RTlna_phase = H_G.dot(cov_X_phase).dot(H_G.T)
            cov_RTlna[j : j + n_mbrs_phase, j : j + n_mbrs_phase] = cov_RTlna_phase
            j += n_mbrs_phase
        else:  # Add formula if it is an endmember
            cov_RTlna[j, j] = 1.0e-12  # small uncertainty for pure phases
            j += 1
    return cov_RTlna


def estimate_conditions(
    assemblage,
    dataset_covariances,
    guessed_conditions=np.array([1.0e9, 873.0]),
    pressure_bounds=[1.0e5, 400.0e9],
    temperature_bounds=[300.0, 4000.0],
    P_scaling=1.0e6,
    small_fraction_tol=0.0,
    max_it=100,
):
    """
    Perform the avPT thermobarometric inversion to find the optimal pressure and temperature
    for a given mineral assemblage. Algorithm based on Powell and Holland (1994).

    :param assemblage: The mineral assemblage for which to perform the inversion.
    :type assemblage: Assemblage

    :param dataset_covariances: The covariance data from the thermodynamic dataset.
    :type dataset_covariances: dict, with keys 'endmember_names' and 'covariance_matrix'

    :param guessed_conditions: Initial guess for pressure (Pa) and temperature (K).
    :type guessed_conditions: np.array of shape (2,)

    :param pressure_bounds: Bounds for pressure (Pa) during optimization.
    :type pressure_bounds: list of two floats

    :param temperature_bounds: Bounds for temperature (K) during optimization.
    :type temperature_bounds: list of two floats

    :param P_scaling: Scaling factor for pressure to improve numerical stability.
    :type P_scaling: float

    :param small_fraction_tol: If > 0.0, reduces the number of endmembers in solution phases by
        transforming to a smaller set of independent endmembers using a greedy algorithm
        and excluding those with molar fractions smaller than this value during the inversion.
    :type small_fraction_tol: float, optional, default 0.0

    :param max_it: Maximum number of iterations for the optimization algorithm.
    :type max_it: int, optional, default 100

    :return: Result object from the optimization containing optimal P and T (x) and other properties
        of the solution. These include the covariance matrix of the estimated parameters (xcov),
        the correlation coefficient between P and T (xcorr),
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
    # TODO: Implement compositional unknown parameters like redox state
    # to join P and T in the optimization

    P_fixed = np.isclose(pressure_bounds[0], pressure_bounds[1])
    T_fixed = np.isclose(temperature_bounds[0], temperature_bounds[1])
    if P_fixed and T_fixed:
        raise Exception("Both pressure and temperature cannot be fixed!")

    # Build the reaction matrix R, possibly reducing the number of endmembers
    # by excluding components with small molar fractions
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
    # i.e. after applying the transformation matrix to the stoichiometric matrix
    S = sp.Matrix(transformation_matrix.dot(assemblage.stoichiometric_array))
    S = S.applyfunc(lambda x: sp.nsimplify(x))
    R = np.array([v[:] for v in S.T.nullspace()], dtype=float)
    R = R.dot(transformation_matrix)  # transform back to original endmember basis
    if len(R) == 0:
        raise Exception("No independent reactions found for the given assemblage!")

    def calculate_cov_a(assemblage, theta, dataset_covariances):
        assemblage.set_state(*theta)
        cov1 = _build_endmember_covariance_matrix(assemblage, dataset_covariances)
        cov2 = _build_solution_covariance_matrix(assemblage)
        cov_mu = cov1 + cov2
        acov = R.dot(cov_mu).dot(R.T)
        return acov

    def chisqr(args):
        """
        Compute the objective misfit function (chi-squared) for given P and T.
        :param theta: Array of pressure (Pa) and temperature (K).
        :type theta: np.array of shape (2,)

        :return: Chi-squared misfit value.
        :rtype: float
        """
        # These are the chemical potential covariances for all endmembers
        cov_a = calculate_cov_a(
            assemblage, [args[0] * P_scaling, args[1]], dataset_covariances
        )
        a = R.dot(assemblage.endmember_partial_gibbs)
        chisqr = a.T.dot(np.linalg.inv(cov_a)).dot(a)
        return chisqr

    # Minimize the chi-squared function
    x0 = [guessed_conditions[0] / P_scaling, guessed_conditions[1]]
    bounds = [
        (pressure_bounds[0] / P_scaling, pressure_bounds[1] / P_scaling),
        (temperature_bounds[0], temperature_bounds[1]),
    ]
    res = minimize(
        chisqr,
        x0,
        method="SLSQP",
        bounds=bounds,
        options={"maxiter": max_it},
    )

    # Rescale pressure back to Pa
    res.x[0] *= P_scaling
    res.jac[0] /= P_scaling

    # Post-process results
    # First we reset the assemblage to the optimal P and T
    # and recalculate the affinity covariance matrix
    res.acov = calculate_cov_a(assemblage, res.x, dataset_covariances)

    # Compute additional diagnostics
    res.dadx = R.dot(
        np.array(
            [
                assemblage.endmember_partial_volumes,
                -assemblage.endmember_partial_entropies,
            ]
        ).T
    )

    res.xcov = np.linalg.pinv(res.dadx.T.dot(np.linalg.pinv(res.acov)).dot(res.dadx))
    res.xcorr = res.xcov[0, 1] / np.sqrt(res.xcov[0, 0] * res.xcov[1, 1])

    res.affinities = R.dot(assemblage.endmember_partial_gibbs)
    res.weighted_affinities = np.linalg.pinv(sqrtm(res.acov)).dot(res.affinities)

    res.n_reactions = len(res.affinities)
    res.n_params = len(res.x)
    res.degrees_of_freedom = res.n_reactions - res.n_params
    res.reduced_chisqr = res.fun / res.degrees_of_freedom
    res.fit = np.sqrt(res.reduced_chisqr)

    return res
