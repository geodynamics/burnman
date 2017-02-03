import numpy as np
from scipy.stats import t

class NonlinearLeastSquaresFit():
    '''
    This is the base class for optimal 
    nonlinear least squares fits
    The algorithm closely follows the logic in 
    Section 23.1 of Bayesian Probability Theory
    (von der Linden et al., 2014; Cambridge University Press)
    '''

    def __init__(self, x, cov, 
                 mle_tolerance,
                 model,
                 lm_damping = 0.,
                 param_tolerance = 1.e-7,
                 max_lm_iterations = 100.,
                 verbose = False):
        
        self.x_arr = x
        self.cov_arr = cov
        self.mle_tolerance = mle_tolerance
        self.model = model
        self.lm_damping = lm_damping
        self.max_lm_iterations = max_lm_iterations
        self.param_tolerance = param_tolerance
        self.verbose = verbose

        self.n_dimensions = len(x[0])
        self.n_data = len(x)
        self.n_params = len(model.get_params())
        self.dof = self.n_data - self.n_params
        
        self.calculate_best_fit()

    def normalised(self, a, order=2, axis=-1):
        l2 = np.atleast_1d(np.linalg.norm(a, order, axis))
        l2[l2==0] = 1
        return a / np.expand_dims(l2, axis)[0][0]
    
    def abs_line_project(self, M, n):
        n = self.normalised(n)
        return n.dot(M).dot(n.T)

    def _mle_estimate(self, x, x_m, cov):
        n = self.model.normal(x_m)
        var_n = self.abs_line_project(cov, n)
        d = (x_m - x).dot(n)
        x_mle = x + d*((n.dot(cov)).T)/var_n
        return x_mle, d, var_n

    def _find_mle(self):
        x_mle_arr = np.empty_like(self.x_arr)
        residual_arr = np.empty(self.n_data)
        var_arr = np.empty(self.n_data)

        for i, (x, cov) in enumerate(zip(*[self.x_arr, self.cov_arr])):
            x_mle_arr[i] = self.model.function(x)
            x_mle_est, residual_arr[i], var_arr[i] = self._mle_estimate(x, x_mle_arr[i], cov)
            delta_x = x_mle_arr[i] - x
        
            while np.linalg.norm(delta_x) > self.mle_tolerance:
                x_mle_est, residual_arr[i], var_arr[i] = self._mle_estimate(x, x_mle_arr[i], cov)
                x_mle_arr[i] = self.model.function(x_mle_est)
                delta_x = x_mle_arr[i] - x_mle_est

        return x_mle_arr, residual_arr/np.sqrt(var_arr), 1./var_arr

    def calculate_jacobian(self):
        self.jacobian = np.empty((self.n_data, self.n_params))

        diag_delta = np.diag(self.model.delta_params)
        param_values = self.model.get_params()
        for i, value in enumerate(param_values):

            self.model.set_params(param_values - diag_delta[i])
            x_mle_arr, residual_arr_0, weights_0 = self._find_mle()

            self.model.set_params(param_values + diag_delta[i])
            x_mle_arr, residual_arr_1, weights_1 = self._find_mle()
            
            self.jacobian[:,i] = (residual_arr_1 - residual_arr_0)/(2.*diag_delta[i][i])
        self.model.set_params(param_values) # reset params
        return None

    def _update_beta(self, lmbda):
        # Performs a Levenberg-Marquardt iteration
        # Note that if lambda = 0, this is a simple Gauss-Newton iteration
        self.calculate_jacobian()
        self.x_mle, self.weighted_residuals, self.weights = self._find_mle()

        J = self.jacobian # this the weighted Jacobian
        JTJ = J.T.dot(J)
        delta_beta = np.linalg.inv(JTJ + lmbda*np.diag(JTJ)).dot(J.T).dot(self.weighted_residuals)
        f_delta_beta = delta_beta/self.model.delta_params

        self.model.set_params(self.model.get_params() - delta_beta)
        
        return f_delta_beta
    
    def calculate_best_fit(self):
        n_it = -1
        tol_achieved = False
        while tol_achieved == False and n_it < self.max_lm_iterations:
            n_it += 1
            delta_beta = self._update_beta(self.lm_damping)
            tol_achieved = np.min(np.abs(delta_beta)) < self.param_tolerance
        if self.verbose == True:
            print 'Converged in {0:d} iterations'.format(n_it)

        J = self.jacobian
        r = self.weighted_residuals
        self.WSS = r.dot(r.T)
        
        self.popt=self.model.get_params()
        self.pcov = np.linalg.inv(J.T.dot(J))*r.dot(r.T)/self.dof

        # Estimate the noise variance normal to the curve
        self.noise_variance = r.dot(np.diag(1./self.weights)).dot(r.T)/self.dof

    def confidence_prediction_bands(self, x_array, confidence_interval, projection_axis=[]):
        '''
        The delta method states that the variance of a function f with parameters B will be
        f'(B, x) Var(B) f'(B, x)
        where f' is the vector of partial derivatives of the function with respect to B 
        '''
        
        # Check array dimensions
        if len(x_array[0]) != self.n_dimensions:
            raise Exception('Dimensions of each point must be the same as the total number of dimensions')

        
        param_values = self.model.get_params()
        normals = np.empty_like(x_array)
        x_m_0 = np.empty_like(x_array)
        for i, x in enumerate(x_array):
            normals[i] = self.model.normal(x)
            x_m_0[i] = self.model.function(x)
            
        diag_delta = np.diag(self.model.delta_params)
        dxdbeta = np.empty([self.n_params, len(x_array)])      
        for i, value in enumerate(param_values):
            self.model.set_params(param_values + diag_delta[i])

            for j, x in enumerate(x_m_0):
                dxdbeta[i][j] = (self.model.function(x) - x).dot(normals[j])/diag_delta[i][i]

        self.model.set_params(param_values) # reset params
        

        variance = np.empty(len(x_array))
        for i, Gprime in enumerate(dxdbeta.T):
            variance[i] = Gprime.T.dot(self.pcov).dot(Gprime)

        critical_value = t.isf(confidence_interval, self.dof)
        
        confidence_half_widths = critical_value*np.sqrt(variance)
        prediction_half_widths = critical_value*np.sqrt(variance + self.noise_variance)

        if projection_axis == []:
            projection_axis = normals[0]

        norms = np.empty(len(x_array))
        for i, normal in enumerate(normals):
            norms[i] = projection_axis.dot(normal)

        axes = np.array([projection_axis,]*len(x_array))
        
        projected_confidence_half_widths = np.array([confidence_half_widths/norms,]*self.n_dimensions).T * axes
        projected_prediction_half_widths = np.array([prediction_half_widths/norms,]*self.n_dimensions).T * axes

        
        confidence_bound_0 = x_m_0 - projected_confidence_half_widths
        confidence_bound_1 = x_m_0 + projected_confidence_half_widths
        prediction_bound_0 = x_m_0 - projected_prediction_half_widths
        prediction_bound_1 = x_m_0 + projected_prediction_half_widths


        return np.array([confidence_bound_0, confidence_bound_1,
                         prediction_bound_0, prediction_bound_1])
