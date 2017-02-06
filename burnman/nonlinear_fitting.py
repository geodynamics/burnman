import numpy as np
from scipy.stats import t
import itertools
import copy

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def nonlinear_least_squares_fit(model,
                                lm_damping = 0.,
                                param_tolerance = 1.e-7,
                                max_lm_iterations = 100,
                                verbose = False):

    """
    Function to compute the "best-fit" parameters for a model
    by nonlinear least squares fitting.

    The nonlinear least squares algorithm closely follows the logic in 
    Section 23.1 of Bayesian Probability Theory
    (von der Linden et al., 2014; Cambridge University Press).

    Parameters
    ----------
    model : class instance
        Must have the following attributes:
            data : 2D numpy array. 
                Elements of x[i][j] contain the observed position of 
                data point i

            data_covariance : 3D numpy array
                Elements of cov[i][j][k] contain the covariance matrix
                of data point i

            mle_tolerances : numpy array
                The iterations to find the maximum likelihood estimator 
                for each observed data point will stop when mle_tolerances[i] <  
                np.linalg.norm(data_mle[i] - model.function(data_mle[i], flag))

            delta_params : numpy array
                parameter perturbations used to compute the jacobian

        Must also contain the following functions:
            set_params(self, param_values):
                Function to set parameters

            get_params(self):
                Function to get current model parameters

            function(self, x):
                Returns value of model function evaluated at x

            normal(self, x):
                Returns value of normal to the model function 
                evaluated at x

    lm_damping : float (optional, default: 0)
        Levenberg-Marquardt parameter for least squares minimization

    param_tolerance : float (optional, default: 1.e-5)
        Levenberg-Marquardt iterations are terminated when 
        the maximum fractional change in any of the parameters 
        during an iteration drops below this value

    max_lm_iterations : integer (optional, default: 100)
        Maximum number of Levenberg-Marquardt iterations

    verbose : bool
        Print some information to standard output
        

    Attributes added to model
    ----------
    n_dof : integer
        Degrees of freedom of the system
    data_mle : 2D numpy array
        Maximum likelihood estimates of the observed data points 
        on the best-fit curve
    jacobian : 2D numpy array
        d(weighted_residuals)/d(parameter)
    weighted_residuals : numpy array
        Weighted residuals
    weights : numpy array
        1/(data variances normal to the best fit curve)
    WSS : float
        Weighted sum of squares residuals
    popt : numpy array
        Optimized parameters
    pcov : 2D numpy array
        Covariance matrix of optimized parameters
    noise_variance : float
        Estimate of the variance of the data normal to the curve

    This function is available as ``burnman.nonlinear_least_squares_fit``.
    """
    
    def _mle_estimate(x, x_m, cov, flag):
        n = model.normal(x_m, flag)
        var_n = abs_line_project(cov, n)
        d = (x_m - x).dot(n)
        x_mle = x + d*((n.dot(cov)).T)/var_n
        return x_mle, d, var_n
    
    def _find_mle():
        x_mle_arr = np.empty_like(model.data)
        residual_arr = np.empty(n_data)
        var_arr = np.empty(n_data)
        for i, (x, cov, flag) in enumerate(zip(*[model.data, model.data_covariance, model.flags])):
            x_mle_arr[i] = model.function(x, flag)
            x_mle_est, residual_arr[i], var_arr[i] = _mle_estimate(x, x_mle_arr[i], cov, flag)
            delta_x = x_mle_arr[i] - x
        
            while np.linalg.norm(delta_x) > model.mle_tolerances[i]:
                x_mle_est, residual_arr[i], var_arr[i] = _mle_estimate(x, x_mle_arr[i], cov, flag)
                x_mle_arr[i] = model.function(x_mle_est, flag)
                delta_x = x_mle_arr[i] - x_mle_est

        return x_mle_arr, residual_arr/np.sqrt(var_arr), 1./var_arr

    def calculate_jacobian():
        model.jacobian = np.empty((n_data, n_params))
        diag_delta = np.diag(model.delta_params)
        param_values = model.get_params()
        for prm_i, value in enumerate(param_values):

            model.set_params(param_values - diag_delta[prm_i])
            x_mle_arr, residual_arr_0, weights_0 = _find_mle()

            model.set_params(param_values + diag_delta[prm_i])
            x_mle_arr, residual_arr_1, weights_1 = _find_mle()
            
            model.jacobian[:,prm_i] = (residual_arr_1 - residual_arr_0)/(2.*diag_delta[prm_i][prm_i])
        model.set_params(param_values) # reset params

    def _update_beta(lmbda):
        # Performs a Levenberg-Marquardt iteration
        # Note that if lambda = 0, this is a simple Gauss-Newton iteration
        calculate_jacobian()
        model.data_mle, model.weighted_residuals, model.weights = _find_mle()

        J = model.jacobian # this the weighted Jacobian
        JTJ = J.T.dot(J)
        delta_beta = np.linalg.inv(JTJ + lmbda*np.diag(JTJ)).dot(J.T).dot(model.weighted_residuals)
        new_params = model.get_params() - delta_beta
        f_delta_beta = delta_beta/new_params

        model.set_params(new_params)
        
        return f_delta_beta

    n_data = len(model.data)
    n_params = len(model.get_params())
    n_dimensions = len(model.data[:,0])
    model.dof = n_data - n_params

    
    if not hasattr(model, 'flags'):
        model.flags = [None] * n_data

    for n_it in range(max_lm_iterations):
        f_delta_beta = _update_beta(lm_damping)
        max_f = np.max(np.abs(f_delta_beta))
        if verbose == True:
            print('Iteration {0:d}: {1}. Max change in param: {2}'.format(n_it, model.get_params(), max_f))
        if max_f < param_tolerance:
            break
        
    J = model.jacobian
    r = model.weighted_residuals
    model.WSS = r.dot(r.T)
        
    model.popt = model.get_params()
    model.pcov = np.linalg.inv(J.T.dot(J))*r.dot(r.T)/model.dof
    
    # Estimate the noise variance normal to the curve
    model.noise_variance = r.dot(np.diag(1./model.weights)).dot(r.T)/model.dof
    
    if verbose == True:
        if n_it == max_lm_iterations - 1:
            print('Max iterations ({0:d}) reached (param tolerance = {1:1e})'.format(max_lm_iterations, param_tolerance))
        else:
            print('Converged in {0:d} iterations'.format(n_it))
        print('\nOptimised parameter values:')
        print(model.popt)
        print('\nParameter covariance matrix:')
        print(model.pcov)
        print('')


    
def orthogonal_distance_confidence_prediction_bands(model, x_array, confidence_interval, projection_axis=[], flag=None):
    """
    This function calculates the confidence and prediction bands of the orthogonal distance 
    from a best-fit model with uncertainties in its parameters as calculated (for example) 
    by the function nonlinear_least_squares_fit().

    The values are calculated via the delta method, which estimates the variance of a function f 
    evaluated at x as var(f,x) = df(x)/dB var(B) df(x)/dB
    where df(x)/dB is the vector of partial derivatives of f(x) with respect to B 


    Parameters
    ----------
    model : class instance
        As modified (for example) by the function nonlinear_least_squares_fit().
        Should contain the following attributes:
            function, normal, delta_params, pcov, dof, noise_variance

    x_array : 2D numpy array
        coordinates at which to evaluate the bounds

    confidence_interval : float
        Probability level of finding the true model (confidence bound) or any new 
        data point (probability bound). For example, the 95% confidence bounds 
        should be calculated using a confidence interval of 0.95.

    projection_axis : list or numpy array
        Axis along which the confidence bounds should be projected
        For example, we are often interested in plotting the confidence and prediction
        bounds as a function of a single coordinate direction.

    Output
    ------
    bounds : 3D numpy array
        An element of bounds[i][j][k] gives the lower and upper confidence (i=0, i=1) and
        prediction (i=2, i=3) bounds for the jth data point on the kth coordinate axis.
    """
    
    # Check array dimensions
    n_dimensions = len(model.data[0])
    if len(x_array[0]) != n_dimensions:
        raise Exception('Dimensions of each point must be the same as the total number of dimensions')


    param_values = model.get_params()
    normals = np.empty_like(x_array)
    x_m_0 = np.empty_like(x_array)
    for i, x in enumerate(x_array):
        normals[i] = model.normal(x, flag)
        x_m_0[i] = model.function(x, flag)
            
    diag_delta = np.diag(model.delta_params)
    dxdbeta = np.empty([len(param_values), len(x_array)])      
    for i, value in enumerate(param_values):
        model.set_params(param_values + diag_delta[i])

        for j, x in enumerate(x_m_0):
            dxdbeta[i][j] = (model.function(x, flag) - x).dot(normals[j])/diag_delta[i][i]

    model.set_params(param_values) # reset params
        

    variance = np.empty(len(x_array))
    for i, Gprime in enumerate(dxdbeta.T):
        variance[i] = Gprime.T.dot(model.pcov).dot(Gprime)

    critical_value = t.isf(0.5*(confidence_interval + 1.), model.dof)
        
    confidence_half_widths = critical_value*np.sqrt(variance)
    prediction_half_widths = critical_value*np.sqrt(variance + model.noise_variance)

    if projection_axis == []:
        projection_axis = normals[0]

    norms = np.empty(len(x_array))
    for i, normal in enumerate(normals):
        norms[i] = projection_axis.dot(normal)

    axes = np.array([projection_axis,]*len(x_array))
        
    projected_confidence_half_widths = np.array([confidence_half_widths/norms,]*n_dimensions).T * axes
    projected_prediction_half_widths = np.array([prediction_half_widths/norms,]*n_dimensions).T * axes

        
    confidence_bound_0 = x_m_0 - projected_confidence_half_widths
    confidence_bound_1 = x_m_0 + projected_confidence_half_widths
    prediction_bound_0 = x_m_0 - projected_prediction_half_widths
    prediction_bound_1 = x_m_0 + projected_prediction_half_widths
    
    return np.array([confidence_bound_0, confidence_bound_1,
                     prediction_bound_0, prediction_bound_1])
    
def confidence_prediction_bands(model, f, x_array, confidence_interval, flag=None):
    """
    This function calculates the confidence and prediction bands of the function f
    from a best-fit model with uncertainties in its parameters as calculated (for example) 
    by the function nonlinear_least_squares_fit().

    The values are calculated via the delta method, which estimates the variance of f 
    evaluated at x as var(f(x)) = df(x)/dB var(B) df(x)/dB
    where df(x)/dB is the vector of partial derivatives of f(x) with respect to B 


    Parameters
    ----------
    model : class instance
        As modified (for example) by the function nonlinear_least_squares_fit().
        Should contain the following attributes:
            function, normal, delta_params, pcov, dof, noise_variance

    x_array : 2D numpy array
        coordinates at which to evaluate the bounds

    confidence_interval : float
        Probability level of finding the true model (confidence bound) or any new 
        data point (probability bound). For example, the 95% confidence bounds 
        should be calculated using a confidence interval of 0.95.


    Output
    ------
    bounds : 2D numpy array
        An element of bounds[i][j] gives the lower and upper confidence (i=0, i=1) and
        prediction (i=2, i=3) bounds for the jth data point.
    """
    
    # Check array dimensions
    n_dimensions = len(model.data[0])
    if len(x_array[0]) != n_dimensions:
        raise Exception('Dimensions of each point must be the same as the total number of dimensions')

        
    param_values = model.get_params()
    normals = np.empty_like(x_array)
    x_m_0s = np.empty_like(x_array)
    f_m_0s = np.empty_like(x_array[:,0])
    for i, x in enumerate(x_array):
        normals[i] = model.normal(x, flag)
        x_m_0s[i] = model.function(x, flag)
        f_m_0s[i] = f(x)
            
    diag_delta = np.diag(model.delta_params)
    dxdbeta = np.empty([len(param_values), len(x_array)])

    for i, value in enumerate(param_values):
        model.set_params(param_values + diag_delta[i])

        for j, x_m_0 in enumerate(x_m_0s):
            x_m_1 = x_m_0 + ((model.function(x_m_0, flag) - x_m_0).dot(normals[j]))*normals[j]
            dxdbeta[i][j] = (f(x_m_1) - f_m_0s[j])/diag_delta[i][i]

    model.set_params(param_values) # reset params
    
    variance = np.empty(len(x_array))
    for i, Gprime in enumerate(dxdbeta.T):
        variance[i] = Gprime.T.dot(model.pcov).dot(Gprime)

    critical_value = t.isf(0.5*(confidence_interval + 1.), model.dof)
        
    confidence_half_widths = critical_value*np.sqrt(variance)
    prediction_half_widths = critical_value*np.sqrt(variance + model.noise_variance)
        
    confidence_bound_0 = f_m_0s - confidence_half_widths
    confidence_bound_1 = f_m_0s + confidence_half_widths
    prediction_bound_0 = f_m_0s - prediction_half_widths
    prediction_bound_1 = f_m_0s + prediction_half_widths
    

    return np.array([confidence_bound_0, confidence_bound_1,
                     prediction_bound_0, prediction_bound_1])

def normalised(a, order=2, axis=-1):
    l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)[0][0]

def abs_line_project(M, n):
    n = normalised(n)
    return n.dot(M).dot(n.T)

def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    """
    Plots an `nstd` sigma error ellipse based on the specified covariance
    matrix (`cov`). Additional keyword arguments are passed on to the 
    ellipse patch artist.

    Parameters
    ----------
        cov : The 2x2 covariance matrix to base the ellipse on
        pos : The location of the center of the ellipse. Expects a 2-element
            sequence of [x0, y0].
        nstd : The radius of the ellipse in numbers of standard deviations.
            Defaults to 2 standard deviations.
        ax : The axis that the ellipse will be plotted on. Defaults to the 
            current axis.
        Additional keyword arguments are pass on to the ellipse patch.

    Returns
    -------
        A matplotlib ellipse artist
    """
    def eigsorted(cov):
        vals, vecs = np.linalg.eigh(cov)
        order = vals.argsort()[::-1]
        return vals[order], vecs[:,order]

    if ax is None:
        ax = plt.gca()
        
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))

    # Width and height are "full" widths, not radius
    width, height = 2 * nstd * np.sqrt(vals)
    ellip = Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)

    ax.add_artist(ellip)
    return ellip


def corner_plot(popt, pcov, param_names=[]):
    """
    Creates a corner plot of covariances

    Parameters
    ----------


    Returns
    -------
    fig : Instance of matplotlib.pyplot.figure

    """
    
    if len(pcov[0]) != len(pcov[:,0]):
        raise Exception('Covariance matrices must be square')

    n_params = len(pcov[0])
    if n_params < 2:
        raise Exception('Covariance matrix must be at least 2x2 for a corner plot to be plotted')

    
    # ellipse plotting is prone to rounding errors, so we scale the plots here
    scaling = 1./np.power(10., np.around(np.log10(np.abs(popt)) - 0.5))
    scaling = np.outer(scaling, scaling)

    fig = plt.figure()
    
    for i in range(n_params):
        for j in range(i+1, n_params):
            indices = np.array([i, j])
            projected_cov = (pcov*scaling)[indices[:, None], indices]

            scaled_pos = np.array([popt[i]*np.sqrt(scaling[i][i]),
                                   popt[j]*np.sqrt(scaling[j][j])])

            nstd = 1.
            ax = fig.add_subplot(n_params-1, n_params-1, (j-1)*(n_params-1)+(i+1))
            plot_cov_ellipse(cov=projected_cov, pos=scaled_pos,
                             nstd=nstd, ax=ax, color='grey')
            maxx = 1.5*nstd*np.sqrt(projected_cov[0][0])
            maxy = 1.5*nstd*np.sqrt(projected_cov[1][1])
            ax.set_xlim(scaled_pos[0]-maxx, scaled_pos[0]+maxx)
            ax.set_ylim(scaled_pos[1]-maxy, scaled_pos[1]+maxy)

    if param_names != []:
        for i in range(n_params-1):
            ax = fig.add_subplot(n_params-1, n_params-1, (n_params-2)*(n_params-1)+(i+1))
            ax.set_xlabel('{0:s} (x 10^{1:d})'.format(param_names[i], -int(np.log10(np.sqrt(scaling[i][i])))))
        for j in range(n_params-1):
            ax = fig.add_subplot(n_params-1, n_params-1, j*(n_params-1)+1)
            ax.set_ylabel('{0:s} (x 10^{1:d})'.format(param_names[j+1], -int(np.log10(np.sqrt(scaling[j+1][j+1])))))
            
    return fig
