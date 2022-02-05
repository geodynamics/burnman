# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
eos_fitting_solution
--------------------

This script fits parameters from the the Stixrude and Lithgow-Bertelloni (2011)
equation of state for solid solutions to mocked-up olivine volume data.
"""

import numpy as np
import matplotlib.pyplot as plt
from numpy import random

import burnman
from burnman.tools.misc import pretty_print_values
from burnman.optimize.eos_fitting import fit_XPTp_data
from burnman.optimize.nonlinear_fitting import plot_residuals, extreme_values
from burnman.optimize.nonlinear_fitting import corner_plot, weighted_residual_plot

if __name__ == "__main__":

    # Solution to optimise
    # (along with tweaks to initial properties if necessary)
    solution = burnman.minerals.SLB_2011.mg_fe_olivine()
    solution.set_state(1.e5, 300.)

    print(solution.endmember_names)

    # Fit parameters
    fit_params = [['V_0', 0],
                  ['V_0', 1],
                  ['K_0', 0],
                  ['K_0', 1],
                  ['Kprime_0', 0],
                  ['Kprime_0', 1],
                  ['V', 0, 1]]

    delta_params = np.array([1.e-8, 1.e-8, 1.e7, 1.e7, 1.e-1, 1.e-1, 1.e-8])
    bounds = np.array([[0, np.inf],
                       [0, np.inf],
                       [0, np.inf],
                       [0, np.inf],
                       [3.5, 6.],
                       [3.5, 6.],
                       [-np.inf, np.inf]])

    # make up some data
    n_data = 100
    data = []
    data_covariances = []

    f_Verror = 1.e-3

    random.seed(10)
    for i in range(n_data):
        x_fa = random.random()
        P = random.random() * 1.e10
        T = random.random() * 1000. + 300.
        X = [1.-x_fa, x_fa]
        solution.set_composition(X)
        solution.set_state(P, T)
        f = (1. + (random.normal() - 0.5)*f_Verror)
        V = solution.V * f

        data.append([1.-x_fa, x_fa, P, T, V])
        data_covariances.append(np.zeros((5, 5)))
        data_covariances[-1][4, 4] = np.power(solution.V*f_Verror, 2.)

    # add one awful data point based on the last one
    data.append([1.-x_fa, x_fa, P, T, V + 2.e-7])
    data_covariances.append(np.zeros((5, 5)))
    data_covariances[-1][4, 4] = np.power(solution.V*f_Verror, 2.)

    data = np.array(data)
    data_covariances = np.array(data_covariances)
    param_tolerance = 1.e-8
    flags = 'V'

    confidence_interval = 0.95
    remove_outliers = True
    good_data_confidence_interval = 0.99
    param_tolerance = 1.e-5

    properties_for_data_comparison_plots = [('V', 1.e6, 'Volume (cm^3/mol)')]

    print('Starting to fit user-defined data. Please be patient.')
    fitted_eos = fit_XPTp_data(solution=solution,
                               flags=flags,
                               fit_params=fit_params,
                               data=data,
                               data_covariances=data_covariances,
                               delta_params=delta_params,
                               bounds=bounds,
                               param_tolerance=param_tolerance,
                               verbose=False)

    # Print the optimized parameters
    print('Optimized equation of state:')
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov,
                        fitted_eos.fit_params_strings)
    print('\nParameters:')
    print(fitted_eos.popt)
    print('\nFull covariance matrix:')
    print(fitted_eos.pcov)
    print('\nGoodness of fit:')
    print(fitted_eos.goodness_of_fit)
    print('\n')

    # Create a plot of the residuals
    fig, ax = plt.subplots()
    plot_residuals(ax=ax,
                   weighted_residuals=fitted_eos.weighted_residuals,
                   flags=fitted_eos.flags)
    plt.show()

    val_extreme = extreme_values(fitted_eos.weighted_residuals,
                                 good_data_confidence_interval)
    confidence_bound, indices, probabilities = val_extreme

    if indices != [] and remove_outliers is True:
        print(f'Removing {len(indices):d} outliers'
              f' (at the {good_data_confidence_interval*100.:.1f}% '
              'confidence interval) and refitting. '
              'Please wait just a little longer.')

        mask = [i for i in range(len(fitted_eos.weighted_residuals))
                if i not in indices]
        data = data[mask]
        data_covariances = data_covariances[mask]
        fitted_eos = fit_XPTp_data(solution=solution,
                                   flags=flags,
                                   fit_params=fit_params,
                                   data=data,
                                   data_covariances=data_covariances,
                                   param_tolerance=param_tolerance,
                                   verbose=False)

    # Print the optimized parameters
    print('Optimized equation of state:')
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov,
                        fitted_eos.fit_params_strings)
    print('\nParameters:')
    print(fitted_eos.popt)
    print('\nFull covariance matrix:')
    print(fitted_eos.pcov)
    print('\nGoodness of fit:')
    print(fitted_eos.goodness_of_fit)
    print('\n')

    # Create a plot of the residuals
    fig, ax = plt.subplots()
    plot_residuals(ax=ax,
                   weighted_residuals=fitted_eos.weighted_residuals,
                   flags=fitted_eos.flags)
    plt.show()

    # Create a corner plot of the covariances
    fig, ax_array = corner_plot(popt=fitted_eos.popt,
                                pcov=fitted_eos.pcov,
                                param_names=fitted_eos.fit_params_strings)
    plt.show()

    # Create plots for the weighted residuals of each type of measurement
    enum_prps = enumerate(properties_for_data_comparison_plots)
    for i, (material_property, scaling, name) in enum_prps:
        fig, ax = plt.subplots()
        weighted_residual_plot(ax=ax,
                               model=fitted_eos,
                               flag=material_property,
                               sd_limit=3,
                               cmap=plt.cm.RdYlBu,
                               plot_axes=[2, 3],
                               scale_axes=[1.e-9, 1.])
        ax.set_title(f'Weighted residual plot for {name:s}')
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Temperature (K)')
        plt.show()
