# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np

from . import nonlinear_fitting
from .tools import flatten, unit_normalize

def fit_PTp_data(mineral, fit_params, flags, data, data_covariances=[], mle_tolerances=[], param_tolerance=1.e-5, max_lm_iterations=50, verbose=True):
    """
    Given a mineral of any type, a list of fit parameters
    and a set of P-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    Parameters
    ----------
    mineral : mineral
        Mineral for which the parameters should be optimized.

    fit_params : list of strings
        List of dictionary keys contained in mineral.params
        corresponding to the variables to be optimized
        during fitting. Initial guesses are taken from the existing
        values for the parameters

    flags : string or list of strings
        Attribute names for the property to be fit for the whole 
        dataset or each datum individually (e.g. 'V')

    data : numpy array of observed P-T-property values

    data_covariances : numpy array of P-T-property covariances (optional)
        If not given, all covariance matrices are chosen 
        such that C00 = 1, otherwise Cij = 0
        In other words, all data points have equal weight, 
        with all error in the pressure

    Returns
    -------
    model : instance of fitted model
        Fitting-related attributes are as follows:
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
    """

    class Model(object):
        def __init__(self, mineral, data, data_covariances, flags, fit_params, guessed_params, delta_params, mle_tolerances):
            self.m = mineral
            self.data = data
            self.data_covariances = data_covariances
            self.flags = flags
            self.fit_params = fit_params
            self.set_params(guessed_params)
            self.delta_params = delta_params
            self.mle_tolerances = mle_tolerances
            
        def set_params(self, param_values):
            i=0
            for param in self.fit_params:
                if isinstance(self.m.params[param], float):
                    self.m.params[param] = param_values[i]
                    i += 1
                else:
                    for j in range(len(self.m.params[param])):
                        self.m.params[param][j] = param_values[i]
                        i += 1

        def get_params(self):
            params = []
            for i, param in enumerate(self.fit_params):
                params.append(self.m.params[param])
            return np.array(flatten([mineral.params[prm] for prm in fit_params]))

        def function(self, x, flag):
            P, T, p = x
            self.m.set_state(P, T)
            return np.array([P, T, getattr(self.m, flag)])

        def normal(self, x, flag):
            P, T, p = x
            
            if flag == 'V':
                self.m.set_state(P, T)
                dPdp = -self.m.K_T/self.m.V
                dpdT = self.m.alpha*self.m.V
            elif flag == 'H':
                self.m.set_state(P, T)
                dPdp = 1./((1.-T*self.m.alpha)*self.m.V)
                dpdT = self.m.molar_heat_capacity_p
            elif flag == 'S':
                self.m.set_state(P, T)
                dPdp = -1./(self.m.alpha*self.m.V)
                dpdT = self.m.molar_heat_capacity_p/T
            elif flag == 'gibbs':
                self.m.set_state(P, T)
                dPdp = 1./self.m.V
                dpdT = -self.S
            else:
                dP = 1.e5
                dT = 1.
                dPdp = (2.*dP)/(self.function([P+dP, T, 0.], flag)[2] - self.function([P-dP, T, 0.], flag)[2])
                dpdT = (self.function([P, T+dT, 0.], flag)[2] - self.function([P, T-dT, 0.], flag)[2])/(2.*dT)
            dPdT = -dPdp*dpdT    
            n = np.array([-1., dPdT, dPdp])
            return unit_normalize(n)

        
    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:,0]))

    # Apply mle tolerances if they dont exist
    if mle_tolerances == []:
        mineral.set_state(1.e5, 300.)
        mle_tolerance_factor = 1.e-5
        mle_tolerances = np.empty(len(flags))
        for i, flag in enumerate(flags):
            if flag in ['gibbs', 'enthalpy', 'H', 'helmholtz']:
                mle_tolerances[i] = 1. # 1 J
            else:
                mle_tolerances[i] = mle_tolerance_factor*getattr(mineral, flag)
        
    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and p)
    covariances_defined = True
    if data_covariances == []:
        covariances_defined = False
        data_covariances = np.zeros((len(data[:,0]), len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][0][0] = 1.
    
    guessed_params = np.array(flatten([mineral.params[prm] for prm in fit_params]))
    model = Model(mineral = mineral,
                  data = data,
                  data_covariances = data_covariances,
                  flags = flags,
                  fit_params = fit_params,
                  guessed_params = guessed_params,
                  delta_params = guessed_params*1.e-5,
                  mle_tolerances = mle_tolerances)

    nonlinear_fitting.nonlinear_least_squares_fit(model, max_lm_iterations = max_lm_iterations, param_tolerance=param_tolerance, verbose=verbose)

    if verbose == True and covariances_defined == True:
        confidence_interval = 0.9
        confidence_bound, indices, probabilities = nonlinear_fitting.extreme_values(model.weighted_residuals, confidence_interval)
        if indices != []:
            print('The function nonlinear_fitting.extreme_values(model.weighted_residuals, confidence_interval) '
                  'has determined that there are {0:d} data points which have residuals which are not '
                  'expected at the {1:.1f}% confidence level (> {2:.1f} s.d. away from the model fit).\n'
                  'Their indices and the probabilities of finding such extreme values are:'.format(len(indices), confidence_interval*100., confidence_bound))
            for i, idx in enumerate(indices):
                print('[{0:d}]: {1:.4f} ({2:.1f} s.d. from the model)'.format(idx, probabilities[i], np.abs(model.weighted_residuals[idx])))
            print('You might consider removing them from your fit, '
                  'or increasing the uncertainties in their measured values.\n')
        
    return model


def fit_PTV_data(mineral, fit_params, data, data_covariances=[], param_tolerance=1.e-5, max_lm_iterations=50, verbose=True):
    """
    A simple alias for the fit_PTp_data for when all the data is volume data
    """
        
    return fit_PTp_data(mineral=mineral, flags='V',
                        data=data, data_covariances=data_covariances, 
                        fit_params=fit_params, param_tolerance=param_tolerance, max_lm_iterations=max_lm_iterations, verbose=verbose)
