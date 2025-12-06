# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


from abc import ABC, abstractmethod
import numpy as np
from scipy.stats import t, norm, genextreme
import copy

from ..utils.math import unit_normalize
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse


class NonLinearModel(ABC):
    """
    Abstract base class that defines the required interface
    for models used in nonlinear least squares fitting. Models must support
    evaluation, normal vector computation, and parameter handling.
    """

    @abstractmethod
    def get_params(self):
        """
        :returns: Current model parameters as a NumPy array.
        :rtype: numpy.ndarray
        """
        raise NotImplementedError("get_params() must be implemented by subclass")

    @abstractmethod
    def set_params(self, param_values):
        """
        :param param_values: Parameter values to set.
        :type param_values: numpy.ndarray
        """
        raise NotImplementedError("set_params() must be implemented by subclass")

    @abstractmethod
    def function(self, x, flag=None):
        """
        :param x: Input coordinates.
        :type x: numpy.ndarray

        :param flag: Optional argument to control evaluation mode.
        :type flag: any

        :returns: Model output.
        :rtype: numpy.ndarray
        """
        raise NotImplementedError("function() must be implemented by subclass")

    @abstractmethod
    def normal(self, x, flag=None):
        """
        :param x: Input coordinates.
        :type x: numpy.ndarray

        :param flag: Optional argument to control evaluation mode.
        :type flag: any

        :returns: Normal vector.
        :rtype: numpy.ndarray
        """
        raise NotImplementedError("normal() must be implemented by subclass")

    def validate(self):
        """
        Ensure that all required model attributes are present before fitting.
        These include observed data, data covariance matrices, convergence
        tolerances for MLE projections, and finite difference step sizes.

        :raises AssertionError: If required attributes are missing.
        """
        assert hasattr(self, "data")
        assert hasattr(self, "data_covariances")
        assert hasattr(self, "mle_tolerances")
        assert hasattr(self, "delta_params")


def abs_line_project(M, n):
    """
    Project a covariance matrix M onto a unit direction vector n.
    This gives the variance of the distribution in the direction n.

    :param M: Covariance matrix.
    :type M: numpy.ndarray

    :param n: Direction vector.
    :type n: numpy.ndarray

    :returns: Projected variance.
    :rtype: float
    """
    n = unit_normalize(n)
    return n.dot(M).dot(n.T)


def mle_estimate(model, x, x_m, cov, flag):
    """
    Find an approximation to the maximum likelihood estimate
    of a point on the model surface corresponding to
    noisy observation x, under Gaussian noise with covariance
    cov, given a starting estimate on the model surface x_m.

    This approximation is the orthogonal projection of x
    onto the local tangent plane of the model surface,
    as determined using the Mahalanobis metric defined by cov.

    :param model: Model object.
    :type model: NonLinearModel

    :param x: Observed data point.
    :type x: numpy.ndarray

    :param x_m: Initial guess on the model surface.
    :type x_m: numpy.ndarray

    :param cov: Covariance matrix.
    :type cov: numpy.ndarray

    :param flag: Optional flag for evaluation.
    :type flag: any

    :returns: Tuple containing MLE position, residual, and projected variance.
    :rtype: tuple(numpy.ndarray, float, float)
    """
    # Find the local tangent plane to the model surface.
    n = model.normal(x_m, flag)

    # Signed distance from x to the local tangent plane of
    # the model surface at x_m in the direction of the
    # normal vector n.
    d = (x_m - x).dot(n)

    # Project the covariance matrix M onto the normal
    # to the local tangent plane n.
    var_n = abs_line_project(cov, n)
    x_mle = x + (d / var_n) * cov.dot(n)
    return x_mle, d, var_n


def find_mle(model):
    """
    Find the maximum likelihood point for
    each datum given the current model parameters.

    This is done iteratively using the function mle_estimate.

    :param model: Model object with required attributes.
    :type model: FittableModel

    :returns: MLE data positions, weighted residuals, and weights.
    :rtype: tuple(numpy.ndarray, numpy.ndarray, numpy.ndarray)
    """
    x_mle_arr = np.empty_like(model.data)
    residual_arr = np.empty(len(model.data))
    var_arr = np.empty(len(model.data))

    # For each datum, iterate over the MLE estimates until the
    # global minimum is found.
    for i, (x, cov, flag) in enumerate(
        zip(model.data, model.data_covariances, model.flags)
    ):
        x_mle_arr[i] = model.function(x, flag)
        x_mle_est, residual_arr[i], var_arr[i] = mle_estimate(
            model, x, x_mle_arr[i], cov, flag
        )
        delta_x = x_mle_arr[i] - x
        while np.linalg.norm(delta_x) > model.mle_tolerances[i]:
            x_mle_est, residual_arr[i], var_arr[i] = mle_estimate(
                model, x, x_mle_arr[i], cov, flag
            )
            x_mle_arr[i] = model.function(x_mle_est, flag)
            delta_x = x_mle_arr[i] - x_mle_est
    return x_mle_arr, residual_arr / np.sqrt(var_arr), 1.0 / var_arr


def calculate_jacobian(model):
    """
    Compute the Jacobian matrix of weighted residuals with respect
    to model parameters using central finite differences.

    :param model: Model object.
    :type model: FittableModel

    :modifies: model.jacobian — Populated with partial derivatives.
    """
    n_data = len(model.data)
    n_params = len(model.get_params())
    model.jacobian = np.empty((n_data, n_params))
    diag_delta = np.diag(model.delta_params)
    param_values = model.get_params()

    # Get d (weighted residual) / d parameter by numerical differences
    for i in range(n_params):
        model.set_params(param_values - diag_delta[i])
        _, res_0, _ = find_mle(model)

        model.set_params(param_values + diag_delta[i])
        _, res_1, _ = find_mle(model)

        model.jacobian[:, i] = (res_1 - res_0) / (2.0 * diag_delta[i][i])

    # Reset parameter values
    model.set_params(param_values)


def nonlinear_least_squares_fit(
    model,
    lm_damping=0.0,
    param_tolerance=1.0e-7,
    max_lm_iterations=100,
    param_priors=None,
    param_prior_inv_cov_matrix=None,
    verbose=False,
):
    """
    Function to compute the "best-fit" parameters for a model
    by nonlinear least squares fitting.

    The nonlinear least squares algorithm closely follows the logic in
    Section 23.1 of Bayesian Probability Theory
    (von der Linden et al., 2014; Cambridge University Press).

    :param model: Model with fitting interface.
    :type model: FittableModel

    :param lm_damping: Damping factor for Levenberg-Marquardt updates.
    :type lm_damping: float

    :param param_tolerance: Convergence tolerance based on fractional parameter change.
    :type param_tolerance: float

    :param max_lm_iterations: Maximum number of LM iterations.
    :type max_lm_iterations: int

    :param param_priors: Prior values for the parameters.
    :type param_priors: 1D numpy array

    :param param_prior_inv_cov_matrix: Inverse of the 1 sigma uncertainties for the prior values of the parameters.
    :type param_prior_inv_cov_matrix: 2D numpy array (square)

    :param verbose: If True, print iteration status.
    :type verbose: bool

    :modifies: model — Sets optimized parameters, covariance matrix, weighted residuals, Jacobian, and noise estimates.

    .. note:: The object passed as model must have the following attributes:

        - data [2D numpy.array] - Elements of x[i][j] contain the
            observed position of data point i.
        - data_covariances [3D numpy.array] Elements of cov[i][j][k] contain
            the covariance matrix of data point i.
        - mle_tolerances [numpy.array] - The iterations to find the maximum likelihood
            estimator for each observed data point will stop when mle_tolerances[i] <
            np.linalg.norm(data_mle[i] - model.function(data_mle[i], flag))
        - delta_params [numpy.array] - parameter perturbations used to compute the jacobian

        Must also have the following methods:

        - set_params(self, param_values) -  Function to set parameters.
        - get_params(self) - Function to get current model parameters.
        - function(self, x) - Returns value of model function evaluated at x.
        - normal(self, x) - Returns value of normal to the model function evaluated at x.

        After this function has been performed, the following attributes are added to model:

        - n_dof [int] - Degrees of freedom of the system.
        - data_mle [2D numpy array] - Maximum likelihood estimates of the observed data points
            on the best-fit curve.
        - jacobian [2D numpy array] - d(weighted_residuals)/d(parameter).
        - weighted_residuals [numpy array] - Weighted residuals.
        - weights [numpy array] - 1/(data variances normal to the best fit curve).
        - WSS [float] - Weighted sum of squares residuals.
        - popt [numpy array] - Optimized parameters.
        - pcov [2D numpy array] - Covariance matrix of optimized parameters.
        - noise_variance [float] - Estimate of the variance of the data normal to the curve.

    This function is available as ``burnman.nonlinear_least_squares_fit``.
    """
    model.validate()
    n_data = len(model.data)
    n_params = len(model.get_params())
    model.dof = n_data - n_params

    if not hasattr(model, "flags"):
        model.flags = [None] * n_data

    with_param_priors = False
    if param_priors is not None and param_prior_inv_cov_matrix is not None:
        assert len(param_priors) == n_params
        assert param_prior_inv_cov_matrix.shape == (n_params, n_params)
        with_param_priors = True

    def _update_beta(lmbda):
        # Performs a single Levenberg-Marquardt iteration
        # Step 1: Compute Jacobian matrix of weighted residuals
        # Note that if lmbda = 0, this is a simple Gauss-Newton iteration
        calculate_jacobian(model)

        # Step 2: Compute MLE projections and residuals given the current
        # parameters (does not update parameter values)
        model.data_mle, model.weighted_residuals, model.weights = find_mle(model)

        # Step 2: Build data terms
        current_params = model.get_params()
        J = model.jacobian  # d weighted residuals / d params
        r = model.weighted_residuals
        JTJ = J.T @ J
        JTr = J.T @ r

        # Step 3: Add Gaussian prior if defined
        if with_param_priors:
            prior_residual = current_params - param_priors
            JTJ += param_prior_inv_cov_matrix
            JTr += param_prior_inv_cov_matrix @ prior_residual

        # Step 4: Apply Levenberg-Marquardt update rule
        A = JTJ + lmbda * np.diag(np.diag(JTJ))
        delta_beta = np.linalg.solve(A, JTr)

        # Step 5: Update parameters and compute fractional change
        new_params = current_params - delta_beta
        model.set_params(new_params)

        # set_params may modify the step to satisfy bounds on the problem
        # We therefore need to get the params before
        # calculating the fractional change.
        new_params = model.get_params()

        # In case the new_params object returns a very small value,
        # modify to avoid a pointless comparison:
        mod_params = np.where(
            np.abs(new_params) < param_tolerance, param_tolerance, new_params
        )
        return (current_params - new_params) / mod_params

    for n_it in range(max_lm_iterations):
        try:
            # update the parameters with a LM iteration
            f_delta_beta = _update_beta(lm_damping)
            max_f = np.max(np.abs(f_delta_beta))

            if np.isnan(max_f):
                raise ValueError(
                    "The Levenberg-Marquardt update for "
                    f"Iteration {n_it} was non-numerical."
                )

            if verbose:
                print(f"Iteration {n_it}: max param change = {max_f:.2e}")
            if max_f < param_tolerance:
                break
        except Exception:
            raise Exception(
                f"During non-linear fitting, Iteration {n_it} produced an "
                "exception. This is probably due to numerical failure of "
                "the input model. "
                "Consider imposing bounds on fitting or priors on your "
                "parameter values to prevent this behaviour."
            )

    J = model.jacobian
    r = model.weighted_residuals
    model.WSS = r @ r
    model.popt = model.get_params()
    JTJ = J.T @ J

    if with_param_priors:
        prior_residual = model.popt - param_priors
        model.WSS += prior_residual @ param_prior_inv_cov_matrix @ prior_residual
        JTJ += param_prior_inv_cov_matrix

    model.pcov = np.linalg.inv(JTJ) * model.WSS / model.dof
    model.goodness_of_fit = model.WSS / model.dof
    model.noise_variance = r @ np.diag(1.0 / model.weights) @ r / model.dof

    if verbose:
        print(
            f"Converged in {n_it} iterations"
            if n_it < max_lm_iterations - 1
            else f"Max iterations reached (param tolerance = {param_tolerance:1e})"
        )
        print("\nOptimized parameter values:\n", model.popt)
        print("\nParameter covariance matrix:\n", model.pcov)


def confidence_prediction_bands(model, x_array, confidence_interval, f, flag=None):
    """
    This function calculates the confidence and prediction bands of
    the function f(x) from a best-fit model with uncertainties in its
    parameters as calculated (for example) by
    the function nonlinear_least_squares_fit().

    The values are calculated via the delta method, which estimates
    the variance of f evaluated at x as var(f(x)) = df(x)/dB var(B) df(x)/dB
    where df(x)/dB is the vector of partial derivatives of f(x)
    with respect to B.

    :param model: As modified (for example) by the function
        :func:`burnman.nonlinear_least_squares_fit`.
        Should contain the following functions: get_params, set_params, function, normal
        And attributes: delta_params, pcov, dof, noise_variance
    :type model: object

    :param x_array: Coordinates at which to evaluate the bounds.
    :type x_array: 2D numpy.array

    :param confidence_interval: Probability level of finding the true model
        (confidence bound) or any new data point (probability bound).
        For example, the 95% confidence bounds should be calculated using a
        confidence interval of 0.95.
    :type confidence_interval: float

    :param f: The function defining the variable y=f(x) for which the
        confidence and prediction bounds are desired.
    :type f: function

    :param flag: This (optional) flag is passed to model.function to control how the
        modified position of x is calculated. This value is then used by f(x)
    :type flag: type informed by model object

    :returns: An element of bounds[i][j] gives the lower and upper confidence
        (i=0, i=1) and prediction (i=2, i=3) bounds for the jth data point.
    :rtype: 2D numpy.array
    """

    # Check array dimensions
    n_dimensions = len(model.data[0])
    if len(x_array[0]) != n_dimensions:
        raise Exception(
            "Dimensions of each point must be the same as the "
            "total number of dimensions"
        )

    param_values = model.get_params()
    x_m_0s = np.empty_like(x_array)
    f_m_0s = np.empty_like(x_array[:, 0])
    for i, x in enumerate(x_array):
        x_m_0s[i] = model.function(x, flag)
        f_m_0s[i] = f(x)

    diag_delta = np.diag(model.delta_params)
    dxdbeta = np.empty([len(param_values), len(x_array)])

    for i, value in enumerate(param_values):
        model.set_params(param_values + diag_delta[i])

        for j, x_m_0 in enumerate(x_m_0s):
            x_m_1 = model.function(x_m_0, flag)
            dxdbeta[i][j] = (f(x_m_1) - f_m_0s[j]) / diag_delta[i][i]

    model.set_params(param_values)  # reset params

    variance = np.empty(len(x_array))
    for i, Gprime in enumerate(dxdbeta.T):
        variance[i] = Gprime.T.dot(model.pcov).dot(Gprime)

    critical_value = t.isf(0.5 * (confidence_interval + 1.0), model.dof)

    confidence_half_widths = critical_value * np.sqrt(variance)
    prediction_half_widths = critical_value * np.sqrt(variance + model.noise_variance)

    confidence_bound_0 = f_m_0s - confidence_half_widths
    confidence_bound_1 = f_m_0s + confidence_half_widths
    prediction_bound_0 = f_m_0s - prediction_half_widths
    prediction_bound_1 = f_m_0s + prediction_half_widths

    return np.array(
        [confidence_bound_0, confidence_bound_1, prediction_bound_0, prediction_bound_1]
    )


def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the
    ellipse patch artist.

    :param cov: The 2x2 covariance matrix to base the ellipse on.
    :type cov: numpy.array

    :param pos: The location of the center of the ellipse. Expects a 2-element
        sequence of [x0, y0].
    :type pos: list or numpy.array

    :param nstd: The radius of the ellipse in numbers of standard deviations.
        Defaults to 2 standard deviations.
    :type nstd: float

    :param ax: The axis that the ellipse will be plotted on. Defaults to the
        current axis.
    :type ax: matplotlib.pyplot.axes

    :param kwargs: Additional keyword arguments are passed on to the ellipse patch.

    :returns: The covariance ellipse (already applied to the desired axes object).
    :rtype: matplotlib.patches.Ellipse
    """

    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:, order]

    if ax is None:
        ax = plt.gca()

    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


def corner_plot(popt, pcov, param_names=[], n_std=1.0):
    """
    Creates a corner plot of covariances

    :param popt: Optimized parameters.
    :type popt: numpy.array

    :param pcov: Covariance matrix of the parameters.
    :type pcov: 2D numpy.array

    :param param_names: Parameter names.
    :type param_names: list

    :param n_std: Number of standard deviations for ellipse.
    :type n_std: float

    :returns: ``matplotlib.pyplot.figure`` and list of ``matplotlib.pyplot.Axes``
        objects.
    :rtype: tuple
    """

    if len(pcov[0]) != len(pcov[:, 0]):
        raise Exception("Covariance matrices must be square")

    n_params = len(pcov[0])
    if n_params < 2:
        raise Exception(
            "Covariance matrix must be at least 2x2 for " "a corner plot to be plotted"
        )

    # ellipse plotting is prone to rounding errors, so we scale the plots here
    scaling = 1.0 / np.power(10.0, np.around(np.log10(np.abs(popt)) - 0.5))
    scaling = np.outer(scaling, scaling)

    fig, ax = plt.subplots(n_params - 1, n_params - 1)
    fig.set_size_inches(3.0 * (n_params - 1), 3.0 * (n_params - 1))

    for j in range(n_params - 1):
        for i in range(j):
            fig.delaxes(ax[i][j])

        for i in range(j, n_params - 1):
            ax[i][j].get_xaxis().get_major_formatter().set_useOffset(False)
            ax[i][j].get_yaxis().get_major_formatter().set_useOffset(False)
            ax[i][j].set_box_aspect(1)

            if j > 0:
                ax[i][j].get_yaxis().set_visible(False)
            if i < n_params - 2:
                ax[i][j].get_xaxis().set_visible(False)

            indices = np.array([j, i + 1])
            projected_cov = (pcov * scaling)[indices[:, None], indices]

            scaled_pos = np.array(
                [
                    popt[j] * np.sqrt(scaling[j][j]),
                    popt[i + 1] * np.sqrt(scaling[i + 1][i + 1]),
                ]
            )

            plot_cov_ellipse(
                cov=projected_cov, pos=scaled_pos, nstd=n_std, ax=ax[i][j], color="grey"
            )
            maxx = 1.5 * n_std * np.sqrt(projected_cov[0][0])
            maxy = 1.5 * n_std * np.sqrt(projected_cov[1][1])
            ax[i][j].set_xlim(scaled_pos[0] - maxx, scaled_pos[0] + maxx)
            ax[i][j].set_ylim(scaled_pos[1] - maxy, scaled_pos[1] + maxy)

    if param_names != []:
        for i in range(n_params - 1):
            ax[n_params - 2][i].set_xlabel(
                "{0:s} (x $10^{{{1:d}}}$)".format(
                    param_names[i], -int(np.log10(np.sqrt(scaling[i][i])))
                )
            )

        for j in range(1, n_params):
            ax[j - 1][0].set_ylabel(
                "{0:s} (x $10^{{{1:d}}}$)".format(
                    param_names[j], -int(np.log10(np.sqrt(scaling[j][j])))
                )
            )

    fig.set_layout_engine("tight")

    return fig, ax


def weighted_residual_plot(
    ax,
    model,
    flag=None,
    sd_limit=3,
    cmap=plt.cm.RdYlBu,
    plot_axes=[0, 1],
    scale_axes=[1.0, 1.0],
):
    """
    Creates a plot of the weighted residuals
    The user can choose the projection axes, and scaling to apply to those axes
    The chosen color palette (cmap) is discretised by standard deviation up
    to a cut off value of sd_limit.

    :param ax: Plot.
    :param type: ``matplotlib.pyplot.Axes``

    :param model: A model as used by
        :func:`burnman.nonlinear_least_squares_fit`.
        Must contain the attributes model.data,
        model.weighted_residuals and
        model.flags (if flag is not None).
    :type model: object

    :param flag: String to determine which data to plot.
        Finds matches with model.flags.
    :type flag: str

    :param sd_limit: Data with weighted residuals exceeding this
        limit are plotted in black.
    :type sd_limit: float

    :param cmap: Color palette.
    :type cmap: matplotlib color palette

    :param plot_axes: Data axes to use as plot axes.
    :type plot_axes: list of int

    :param scale_axes: Plot axes are scaled by multiplication
        of the data by these values.
    :type scale_axes: list of float

    :returns: Coloured scatter plot of the weighted residuals in data space.
    :rtype: matplotlib Axes object
    """
    if flag is None:
        mask = range(len(model.data[:, 0]))
    else:
        mask = [i for i, flg in enumerate(model.flags) if flg == flag]

    cmap_cp = copy.copy(cmap)
    cmap_cp.set_under("k")
    cmap_cp.set_over("k")
    bounds = np.linspace(-sd_limit, sd_limit, sd_limit * 2 + 1)
    norm = colors.BoundaryNorm(bounds, cmap_cp.N)

    im = ax.scatter(
        model.data[:, plot_axes[0]][mask] * scale_axes[0],
        model.data[:, plot_axes[1]][mask] * scale_axes[1],
        c=model.weighted_residuals[mask],
        cmap=cmap_cp,
        norm=norm,
        s=50,
    )
    plt.colorbar(im, ax=ax, label="Misfit (standard deviations)")


def extreme_values(weighted_residuals, confidence_interval):
    """
    This function uses extreme value theory to calculate the number of
    standard deviations away from the mean at which we should expect to bracket
    *all* of our n data points at a certain confidence level.

    It then uses that value to identify which (if any) of the data points
    lie outside that region, and calculates the corresponding probabilities
    of finding a data point at least that many standard deviations away.


    :param weighted_residuals: Array of residuals weighted by the square root
        of their variances wr_i = r_i/sqrt(var_i).
    :type weighted_residuals: array of float

    :param confidence_interval: Probability at which all the weighted residuals lie
        within the confidence bounds.
    :type confidence_interval: float

    :returns: Number of standard deviations at which we should expect to encompass
        all data at the user-defined confidence interval, the indices of weighted
        residuals exceeding the confidence_interval defined by the user, and
        the probabilities that the extreme data point of the distribution lies
        further from the mean than the observed position wr_i for each i in
        the "indices" output array.
    :rtype: tuple of (float, numpy.array, numpy.array)
    """

    n = len(weighted_residuals)
    mean = norm.isf(1.0 / n)
    # good approximation for > 10 data points
    scale = 0.8 / np.power(np.log(n), 1.0 / 2.0)
    # good approximation for > 10 data points
    c = 0.33 / np.power(np.log(n), 3.0 / 4.0)

    # We now need a 1-tailed probability from the given confidence_interval
    # p_total = 1. - confidence_interval = p_upper + p_lower - p_upper*p_lower
    # p_total = 1. - confidence_interval = 2p - p^2, therefore:
    p = 1.0 - np.sqrt(confidence_interval)
    confidence_bound = genextreme.isf(p, c, loc=mean, scale=scale)

    indices = [
        i for i, r in enumerate(weighted_residuals) if np.abs(r) > confidence_bound
    ]
    # Convert back to 2-tailed probabilities
    probabilities = 1.0 - np.power(
        genextreme.sf(np.abs(weighted_residuals[indices]), c, loc=mean, scale=scale)
        - 1.0,
        2.0,
    )

    return confidence_bound, indices, probabilities


def plot_residuals(ax, weighted_residuals, n_bins=None, flags=[]):
    if flags is []:
        flags = [""] * len(weighted_residuals)
        list_flags = [""]
    else:
        list_flags = list(set(flags))

    if n_bins is None:
        try:  # Only works for recent versions of numpy
            bin_heights, bin_bounds = np.histogram(
                weighted_residuals, bins="auto", density=True
            )
            n_bins = len(bin_heights)
        except:
            n_bins = 11.0

    mask = [i for i, f in enumerate(flags)]
    for flag in list_flags:
        binwidth = np.ptp(weighted_residuals) / n_bins
        dmin = min(weighted_residuals) - binwidth
        dmax = max(weighted_residuals) + binwidth
        bins = np.linspace(dmin, dmax, n_bins)
        bin_heights, bin_bounds = np.histogram(
            weighted_residuals[mask], bins=bins, density=True
        )

        normalisation = float(len(weighted_residuals[mask])) / float(
            len(weighted_residuals)
        )
        bin_centers = (bin_bounds[:-1] + bin_bounds[1:]) / 2.0
        bin_heights = bin_heights * normalisation
        bin_widths = bin_bounds[1] - bin_bounds[0]
        plt.bar(bin_centers, bin_heights, width=bin_widths, label=flag, alpha=0.2)
        mask = [i for i, f in enumerate(flags) if f != flag and i in mask]

        x = np.linspace(bin_bounds[0], bin_bounds[-1], 1001)
        ax.plot(x, norm.pdf(x) * normalisation)

    ax.set_title("Residual plot versus expected normal distribution")
    ax.set_xlabel("Number of standard deviations from the mean")
    ax.set_ylabel("Probability")
    ax.legend(loc="upper right")
