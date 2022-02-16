# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import
from __future__ import print_function

import numpy as np

from . import nonlinear_fitting
from ..utils.misc import flatten
from ..utils.math import unit_normalize
from .nonlinear_fitting import nonlinear_least_squares_fit


def fit_PTp_data(mineral, fit_params, flags, data, data_covariances=[],
                 mle_tolerances=[], param_tolerance=1.e-5,
                 delta_params=None, bounds=None,
                 max_lm_iterations=50, verbose=True):
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

    data : 2D numpy array of observed X-P-T-property values

    data_covariances : 3D numpy array of X-P-T-property covariances (optional)
        If not given, all covariance matrices are chosen
        such that all data points have equal weight,
        with all error in the pressure.

    mle_tolerances : numpy array (optional)
        Tolerances for termination of the maximum likelihood iterations.

    param_tolerance : float (optional)
        Fractional tolerance for termination of the nonlinear optimization.

    delta_params : numpy array (optional)
        Initial values for the change in parameters.

    bounds : 2D numpy array (optional)
        Minimum and maximum bounds for the parameters. The shape must be
        (n_parameters, 2).

    max_lm_iterations : integer (default : 50)
        Maximum number of Levenberg-Marquardt iterations.

    verbose : boolean (default : True)
        Whether to print detailed information about the optimization to screen.

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
        def __init__(self, mineral, data, data_covariances, flags, fit_params,
                     mle_tolerances, delta_params=None, bounds=None):
            self.m = mineral
            self.data = data
            self.data_covariances = data_covariances
            self.flags = flags
            self.fit_params = fit_params
            self.mle_tolerances = mle_tolerances
            if delta_params is None:
                self.delta_params = self.get_params()*1.e-5 + 1.e-10
            else:
                self.delta_params = delta_params
            self.bounds = bounds

        def set_params(self, param_values):
            i = 0

            if self.bounds is not None:
                param_values = np.clip(param_values,
                                       self.bounds[:, 0],
                                       self.bounds[:, 1])

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
            return np.array(flatten([mineral.params[prm]
                                     for prm in fit_params]))

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
                dpdT = -self.m.S
            else:
                dP = 1.e5
                dT = 1.
                dPdp = (2.*dP)/(self.function([P+dP, T, 0.], flag)[2]
                                - self.function([P-dP, T, 0.], flag)[2])
                dpdT = (self.function([P, T+dT, 0.], flag)[2]
                        - self.function([P, T-dT, 0.], flag)[2])/(2.*dT)
            dPdT = -dPdp*dpdT
            n = np.array([-1., dPdT, dPdp])
            return unit_normalize(n)

    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:, 0]))

    if len(flags) != len(data):
        raise Exception(f'The number of flags (n = {len(flags)}) must be equal '
                        f'to the number of data (n = {len(data)}).')

    # Apply mle tolerances if they dont exist
    if mle_tolerances == []:
        mineral.set_state(1.e5, 300.)
        mle_tolerance_factor = 1.e-5
        mle_tolerances = np.empty(len(flags))
        for i, flag in enumerate(flags):
            if flag in ['gibbs', 'enthalpy', 'H', 'helmholtz']:
                mle_tolerances[i] = 1.  # 1 J
            else:
                mle_tolerances[i] = mle_tolerance_factor * getattr(mineral,
                                                                   flag)

    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and p)
    covariances_defined = True
    if data_covariances == []:
        covariances_defined = False
        data_covariances = np.zeros((len(data[:, 0]),
                                     len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][0][0] = 1.

    model = Model(mineral=mineral,
                  data=data,
                  data_covariances=data_covariances,
                  flags=flags,
                  fit_params=fit_params,
                  delta_params=delta_params,
                  mle_tolerances=mle_tolerances,
                  bounds=bounds)

    nonlinear_least_squares_fit(model,
                                max_lm_iterations=max_lm_iterations,
                                param_tolerance=param_tolerance,
                                verbose=verbose)

    if verbose is True and covariances_defined is True:
        confidence_interval = 0.9
        d = nonlinear_fitting.extreme_values(model.weighted_residuals,
                                             confidence_interval)
        confidence_bound, indices, probabilities = d
        if indices != []:
            print('The function nonlinear_fitting.extreme_values'
                  '(model.weighted_residuals, confidence_interval) '
                  f'has determined that there are {len(indices):d} data points'
                  ' which have residuals which are not expected at the '
                  f'{confidence_interval*100.:.1f}% confidence level '
                  f'(> {confidence_bound:.1f} s.d. away from the model fit).\n'
                  'Their indices and the probabilities of finding '
                  'such extreme values are:')
            for i, idx in enumerate(indices):
                print(f'[{idx:d}]: {probabilities[i]:.4f} '
                      f'({np.abs(model.weighted_residuals[idx]):.1f} s.d. '
                      'from the model)')
            print('You might consider removing them from your fit, '
                  'or increasing the uncertainties in their '
                  'measured values.\n')

    return model


def fit_PTV_data(mineral, fit_params,
                 data, data_covariances=[],
                 delta_params=None, bounds=None,
                 param_tolerance=1.e-5, max_lm_iterations=50,
                 verbose=True):
    """
    A simple alias for the fit_PTp_data for when all the data is volume data
    """

    return fit_PTp_data(mineral=mineral, flags='V',
                        data=data, data_covariances=data_covariances,
                        fit_params=fit_params, param_tolerance=param_tolerance,
                        delta_params=delta_params, bounds=bounds,
                        max_lm_iterations=max_lm_iterations, verbose=verbose)


def fit_XPTp_data(solution, fit_params, flags, data, data_covariances=[],
                  mle_tolerances=[], param_tolerance=1.e-5,
                  delta_params=None, bounds=None,
                  max_lm_iterations=50, verbose=True):
    """
    Given a symmetric solution, a list of fit parameters
    and a set of P-T-property points and (optional) uncertainties,
    this function returns a list of optimized parameters
    and their associated covariances, fitted using the
    scipy.optimize.curve_fit routine.

    Parameters
    ----------
    solution : Solution
        Solution for which the parameters should be optimized.

    fit_params : list of lists
        List of lists corresponding to the variables to be optimized
        during fitting. Each list is either of length two or three.
        The first item of length-2 lists should be a
        dictionary key contained in one of the endmember
        mineral.params, and the second item should be the index of
        the endmember in the solution (indexing starts from 0).
        The first item of length-3 lists should be one of 'E', 'S' or
        'V' (the excess energies, entropies or volumes in each binary).
        The second two items should be the indices of the pair of
        endmembers bounding the binary, in ascending order
        (indexing starts from 0). Initial guesses are taken from the existing
        values for the parameters.

    flags : string or list of strings
        Attribute names for the property to be fit for the whole
        dataset or each datum individually (e.g. 'V')

    data : 2D numpy array of observed X-P-T-property values

    data_covariances : 3D numpy array of X-P-T-property covariances (optional)
        If not given, all covariance matrices are chosen
        such that all data points have equal weight,
        with all error in the pressure.

    mle_tolerances : numpy array (optional)
        Tolerances for termination of the maximum likelihood iterations.

    param_tolerance : float (optional)
        Fractional tolerance for termination of the nonlinear optimization.

    delta_params : numpy array (optional)
        Initial values for the change in parameters.

    bounds : 2D numpy array (optional)
        Minimum and maximum bounds for the parameters. The shape must be
        (n_parameters, 2).

    max_lm_iterations : integer (default : 50)
        Maximum number of Levenberg-Marquardt iterations.

    verbose : boolean (default : True)
        Whether to print detailed information about the optimization to screen.

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
            fit_params : list of lists
                The parameter lists in their original form
            fit_params_strings : list of strings
                The parameters in user-readable string form
            popt : numpy array
                Optimized parameters
            pcov : 2D numpy array
                Covariance matrix of optimized parameters
            noise_variance : float
                Estimate of the variance of the data normal to the curve
    """

    class Model(object):
        def __init__(self, solution, data, data_covariances, flags, fit_params,
                     mle_tolerances, delta_params=None, bounds=None):
            self.m = solution
            self.data = data
            self.data_covariances = data_covariances
            self.flags = flags
            self.fit_params = fit_params
            self.fit_params_strings = []
            for p in fit_params:
                if isinstance(p, list):
                    csv_list_mbrs = ",".join([str(i) for i in p[1:]])
                    self.fit_params_strings.append(f'{p[0]} ({csv_list_mbrs})')
                else:
                    self.fit_params_strings.append(p)

            self.mle_tolerances = mle_tolerances
            if delta_params is None:
                self.delta_params = self.get_params()*1.e-5 + 1.e-10
            else:
                self.delta_params = delta_params
            self.bounds = bounds

        def set_params(self, param_values):
            # fit_params is a list of lists
            # if the list has length 2, the first item should be an integer
            # indicating the endmember number in the solution
            # if the list has length 3, the first two items should be endmember
            # numbers, and the third should be the interaction parameter type
            # (E, S or V).
            i = 0

            if self.bounds is not None:
                param_values = np.clip(param_values,
                                       self.bounds[:, 0],
                                       self.bounds[:, 1])

            for param in self.fit_params:
                value = param_values[i]
                if len(param) == 2:
                    key, imbr = param
                    if isinstance(self.m.endmembers[imbr][0].params[key],
                                  float):
                        self.m.endmembers[imbr][0].params[key] = value
                        i += 1
                    else:
                        n_values = len(self.m.endmembers[imbr][0].params[key])
                        for j in range(n_values):
                            self.m.endmembers[imbr][0].params[key][j] = value
                            i += 1
                elif len(param) == 3:
                    key, imbr, jmbr = param
                    ai = self.m.solution_model.alphas[imbr]
                    aj = self.m.solution_model.alphas[jmbr]
                    if key == 'E':
                        self.m.solution_model.We[imbr, jmbr] = 2.*value/(ai*aj)
                    if key == 'S':
                        self.m.solution_model.Ws[imbr, jmbr] = 2.*value/(ai*aj)
                    if key == 'V':
                        self.m.solution_model.Wv[imbr, jmbr] = 2.*value/(ai*aj)

                    i += 1
                else:
                    raise Exception('param length must be two or three')

        def get_params(self):
            params = []
            for param in self.fit_params:
                if len(param) == 2:
                    key, imbr = param
                    value = self.m.endmembers[imbr][0].params[key]
                    if isinstance(value, float):
                        params.append(value)
                    else:
                        params.extend(list(value))

                elif len(param) == 3:
                    key, imbr, jmbr = param
                    ai = self.m.solution_model.alphas[imbr]
                    aj = self.m.solution_model.alphas[jmbr]
                    if key == 'E':
                        params.append(self.m.solution_model.We[imbr, jmbr]
                                      * (ai * aj)/2.)
                    if key == 'S':
                        params.append(self.m.solution_model.Ws[imbr, jmbr]
                                      * (ai * aj)/2.)
                    if key == 'V':
                        params.append(self.m.solution_model.Wv[imbr, jmbr]
                                      * (ai * aj)/2.)
                else:
                    raise Exception('param length must be two or three')
            return np.array(params)

        def function(self, x, flag):
            self.m.set_composition(x[:self.m.n_endmembers])
            P, T, p = x[self.m.n_endmembers:]
            self.m.set_state(P, T)

            f = np.copy(x)
            f[-1] = getattr(self.m, flag)
            return f

        def normal(self, x, flag):
            self.m.set_composition(x[:self.m.n_endmembers])
            P, T, p = x[self.m.n_endmembers:]

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
                dpdT = -self.m.S
            else:
                dP = 1.e5
                dT = 1.
                xP0 = np.copy(x)
                xP1 = np.copy(x)
                xT0 = np.copy(x)
                xT1 = np.copy(x)
                xP0[-3] = xP1[-3] - dP
                xP1[-3] = xP1[-3] + dP
                xT0[-2] = xP1[-2] - dT
                xT1[-2] = xP1[-2] + dT

                dPdp = (2.*dP)/(self.function(xP1, flag)[2]
                                - self.function(xP0, flag)[2])
                dpdT = (self.function(xT1, flag)[2]
                        - self.function(xT0, flag)[2])/(2.*dT)
            dPdT = -dPdp*dpdT
            n = np.zeros(len(x))
            n[-3:] = np.array([-1., dPdT, dPdp])
            return unit_normalize(n)

    # If only one property flag is given, assume it applies to all data
    if type(flags) is str:
        flags = np.array([flags] * len(data[:, 0]))

    if len(flags) != len(data):
        raise Exception(f'The number of flags (n = {len(flags)}) must be equal '
                        f'to the number of data (n = {len(data)}).')

    # Apply mle tolerances if they dont exist
    if mle_tolerances == []:
        solution.set_state(1.e5, 300.)
        mle_tolerance_factor = 1.e-5
        mle_tolerances = np.empty(len(flags))
        for i, flag in enumerate(flags):
            if flag in ['gibbs', 'enthalpy', 'H', 'helmholtz']:
                mle_tolerances[i] = 1.  # 1 J
            else:
                mle_tolerances[i] = mle_tolerance_factor*getattr(solution,
                                                                 flag)

    # If covariance matrix is not given, apply unit weighting to all pressures
    # (with zero errors on T and property)
    covariances_defined = True
    if data_covariances == []:
        covariances_defined = False
        nX = solution.n_endmembers
        data_covariances = np.zeros((len(data[:, 0]),
                                     len(data[0]), len(data[0])))
        for i in range(len(data_covariances)):
            data_covariances[i][nX][nX] = 1.

    model = Model(solution=solution,
                  data=data,
                  data_covariances=data_covariances,
                  flags=flags,
                  fit_params=fit_params,
                  mle_tolerances=mle_tolerances,
                  delta_params=delta_params,
                  bounds=bounds)

    nonlinear_least_squares_fit(model,
                                max_lm_iterations=max_lm_iterations,
                                param_tolerance=param_tolerance,
                                verbose=verbose)

    if verbose is True and covariances_defined is True:
        confidence_interval = 0.9
        v_extreme = nonlinear_fitting.extreme_values(model.weighted_residuals,
                                                     confidence_interval)
        confidence_bound, indices, probabilities = v_extreme
        if indices != []:
            print('The function nonlinear_fitting.extreme_values'
                  '(model.weighted_residuals, confidence_interval) '
                  f'has determined that there are {len(indices):d} '
                  'data points which have residuals which are not '
                  f'expected at the {confidence_interval*100.:.1f}% '
                  'confidence level '
                  f'(> {confidence_bound:.1f} s.d. away from the model fit).\n'
                  'Their indices and the probabilities of '
                  'finding such extreme values are:')
            for i, idx in enumerate(indices):
                print(f'[{idx:d}]: {probabilities[i]:.4f} '
                      f'({np.abs(model.weighted_residuals[idx]):.1f} s.d. '
                      'from the model)')
            print('You might consider removing them from your fit, '
                  'or increasing the uncertainties in '
                  'their measured values.\n')

    return model
