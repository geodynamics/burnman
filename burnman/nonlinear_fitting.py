import numpy as np
from scipy.stats import t
import itertools
    
def nonlinear_least_squares_fit(model,
                                mle_tolerance,
                                lm_damping = 0.,
                                param_tolerance = 1.e-7,
                                max_lm_iterations = 100,
                                verbose = False):

    """
    Base class for optimal nonlinear least squares fitting.

    The nonlinear least squares algorithm closely follows the logic in 
    Section 23.1 of Bayesian Probability Theory
    (von der Linden et al., 2014; Cambridge University Press).

    Initialisation inputs
    ---------------------
    x : 2D numpy array. 
        Elements of x[i][j] contain the observed position of 
        data point i

    cov : 3D numpy array
        Elements of cov[i][j][k] contain the covariance matrix
        of data point i

    mle_tolerance : float
        

    model : class instance
        Must contain the following functions:
            set_params(self, param_values):
                Function to set parameters

            get_params(self):
                Function to get current model parameters

            function(self, x):
                Returns value of model function evaluated at x

            (self, x):
                Returns value of model function evaluated at x

            normal(self, x):
                Returns value of normal to the model function 
                evaluated at x

    lm_damping : float (optional, default: 0)

    param_tolerance : float (optional, default: 1.e-7)

    max_lm_iterations : integer (optional, default: 100)

    verbose : bool
        

    Attributes
    ----------
    As above, plus 
    n_dimensions : integer
        Number of dimensions
    n_data : integer
        Number of data points
    n_params : integer
        Number of fitting params
    n_dof : integer
        Degrees of freedom of the system
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

    This class is available as ``burnman.NonlinearLeastSquaresFit``.
    """
    def normalised(a, order=2, axis=-1):
        l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
        l2[l2==0] = 1
        return a / np.expand_dims(l2, axis)[0][0]

    def abs_line_project(M, n):
        n = normalised(n)
        return n.dot(M).dot(n.T)

    def _mle_estimate(x, x_m, cov):
        n = model.normal(x_m)
        var_n = abs_line_project(cov, n)
        d = (x_m - x).dot(n)
        x_mle = x + d*((n.dot(cov)).T)/var_n
        return x_mle, d, var_n
    
    def _find_mle():
        x_mle_arr = np.empty_like(model.data)
        residual_arr = np.empty(n_data)
        var_arr = np.empty(n_data)
        for i, (x, cov) in enumerate(zip(*[model.data, model.data_covariance])):
            x_mle_arr[i] = model.function(x)
            x_mle_est, residual_arr[i], var_arr[i] = _mle_estimate(x, x_mle_arr[i], cov)
            delta_x = x_mle_arr[i] - x
        
            while np.linalg.norm(delta_x) > mle_tolerance:
                x_mle_est, residual_arr[i], var_arr[i] = _mle_estimate(x, x_mle_arr[i], cov)
                x_mle_arr[i] = model.function(x_mle_est)
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
        
        return None

    def _update_beta(lmbda):
        # Performs a Levenberg-Marquardt iteration
        # Note that if lambda = 0, this is a simple Gauss-Newton iteration
        calculate_jacobian()
        model.x_mle, model.weighted_residuals, model.weights = _find_mle()

        J = model.jacobian # this the weighted Jacobian
        JTJ = J.T.dot(J)
        delta_beta = np.linalg.inv(JTJ + lmbda*np.diag(JTJ)).dot(J.T).dot(model.weighted_residuals)
        f_delta_beta = delta_beta/model.delta_params

        model.set_params(model.get_params() - delta_beta)
        
        return f_delta_beta

    n_data = len(model.data)
    n_params = len(model.get_params())
    n_dimensions = len(model.data[:,0])
    model.dof = n_data - n_params
    
    n_it = -1
    tol_achieved = False
    while tol_achieved == False and n_it < max_lm_iterations:
        n_it += 1
        delta_beta = _update_beta(lm_damping)
        tol_achieved = np.min(np.abs(delta_beta)) < param_tolerance
        
    if verbose == True:
        print 'Converged in {0:d} iterations'.format(n_it)

    J = model.jacobian
    r = model.weighted_residuals
    model.WSS = r.dot(r.T)
        
    model.popt = model.get_params()
    model.pcov = np.linalg.inv(J.T.dot(J))*r.dot(r.T)/model.dof
    
    # Estimate the noise variance normal to the curve
    model.noise_variance = r.dot(np.diag(1./model.weights)).dot(r.T)/model.dof

    
def orthogonal_distance_confidence_prediction_bands(model, x_array, confidence_interval, projection_axis=[]):
    '''
    The delta method states that the variance of a function f with parameters B will be
    f'(B, x) Var(B) f'(B, x)
    where f' is the vector of partial derivatives of the function with respect to B 
    '''
    
    # Check array dimensions
    n_dimensions = len(model.data[0])
    if len(x_array[0]) != n_dimensions:
        raise Exception('Dimensions of each point must be the same as the total number of dimensions')


    param_values = model.get_params()
    normals = np.empty_like(x_array)
    x_m_0 = np.empty_like(x_array)
    for i, x in enumerate(x_array):
        normals[i] = model.normal(x)
        x_m_0[i] = model.function(x)
            
    diag_delta = np.diag(model.delta_params)
    dxdbeta = np.empty([len(param_values), len(x_array)])      
    for i, value in enumerate(param_values):
        model.set_params(param_values + diag_delta[i])

        for j, x in enumerate(x_m_0):
            dxdbeta[i][j] = (model.function(x) - x).dot(normals[j])/diag_delta[i][i]

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
    
def confidence_prediction_bands(model, function, x_array, confidence_interval):
    '''
    The delta method states that the variance of a function f with parameters B will be
    f'(B, x) Var(B) f'(B, x)
    where f' is the vector of partial derivatives of the function with respect to B 
    '''
    
    # Check array dimensions
    n_dimensions = len(model.data[0])
    if len(x_array[0]) != n_dimensions:
        raise Exception('Dimensions of each point must be the same as the total number of dimensions')

        
    param_values = model.get_params()
    normals = np.empty_like(x_array)
    x_m_0s = np.empty_like(x_array)
    f_m_0s = np.empty_like(x_array[:,0])
    for i, x in enumerate(x_array):
        normals[i] = model.normal(x)
        x_m_0s[i] = model.function(x)
        f_m_0s[i] = function(x)
            
    diag_delta = np.diag(model.delta_params)
    dxdbeta = np.empty([len(param_values), len(x_array)])

    for i, value in enumerate(param_values):
        model.set_params(param_values + diag_delta[i])

        for j, x_m_0 in enumerate(x_m_0s):
            x_m_1 = x_m_0 + ((model.function(x_m_0) - x_m_0).dot(normals[j]))*normals[j]
            dxdbeta[i][j] = (function(x_m_1) - f_m_0s[j])/diag_delta[i][i]

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



def attribute_function(mineral, attributes, powers=[]):
    if type(attributes) is str:
        attributes = [attributes]
    if powers == []:
        powers = [1. for a in attributes]
    def f(x):
        P, T, V = x
        mineral.set_state(P, T)
        value = 1.
        for a, p in zip(*[attributes, powers]):
            value *= np.power(getattr(mineral, a), p)
        return value
    return f
