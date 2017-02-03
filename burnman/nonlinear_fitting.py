import numpy as np

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
        return x_mle, d/np.sqrt(var_n)

    def _find_mle(self):
        x_mle_arr = np.empty_like(self.x_arr)
        residual_arr = np.empty(self.n_data)

        for i, (x, cov) in enumerate(zip(*[self.x_arr, self.cov_arr])):
            x_mle_arr[i] = self.model.function(x)
            x_mle_est, residual_arr[i] = self._mle_estimate(x, x_mle_arr[i], cov)
            delta_x = x_mle_arr[i] - x
        
            while np.linalg.norm(delta_x) > self.mle_tolerance:
                x_mle_est, residual_arr[i] = self._mle_estimate(x, x_mle_arr[i], cov)
                x_mle_arr[i] = self.model.function(x_mle_est)
                delta_x = x_mle_arr[i] - x_mle_est
            
        return x_mle_arr, residual_arr

    def calculate_jacobian(self):
        self.jacobian = np.empty((self.n_data, self.n_params))

        diag_delta = np.diag(self.model.delta_params)
        param_values = self.model.get_params()
        for i, value in enumerate(param_values):

            self.model.set_params(param_values - diag_delta[i])
            x_mle_arr, residual_arr_0 = self._find_mle()

            self.model.set_params(param_values + diag_delta[i])
            x_mle_arr, residual_arr_1 = self._find_mle()
            
            self.jacobian[:,i] = (residual_arr_1 - residual_arr_0)/diag_delta[i][i]
            
        return None

    def _update_beta(self, lmbda):
        # Performs a Levenberg-Marquardt iteration
        # Note that if lambda = 0, this is a simple Gauss-Newton iteration
        self.calculate_jacobian()
        self.x_mle, self.residuals = self._find_mle()
        
        J = self.jacobian
        delta_beta = np.linalg.inv(J.T.dot(J) + lmbda*np.diag(J.T.dot(J))).dot(J.T).dot(self.residuals)
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
        r = self.residuals
        self.popt=self.model.get_params()
        self.noise_variance = r.dot(r.T)/self.dof
        self.pcov = np.linalg.inv(J.T.dot(J))*self.noise_variance
        self.WSS = r.dot(r.T)

