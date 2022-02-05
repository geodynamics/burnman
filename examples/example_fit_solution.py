# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_fit_solution
--------------------

This example demonstrates BurnMan's functionality to fit data
on solutions to a solution model and equations of state
of the user's choice.

The example in this file deals with simple PVT fitting,
but could easily be extended to fitting of other data.

teaches:
- least squares fitting for solution data

"""
import numpy as np
import matplotlib.pyplot as plt
from numpy import random

import burnman
from burnman.tools.misc import pretty_print_values
from burnman.optimize.eos_fitting import fit_XPTp_data
from burnman.optimize.nonlinear_fitting import plot_residuals, extreme_values
from burnman.optimize.nonlinear_fitting import corner_plot
from burnman.optimize.nonlinear_fitting import weighted_residual_plot


if __name__ == "__main__":

    # Set np.array printing precision to be low
    # (for more readable covariance matrices)
    np.set_printoptions(precision=1)

    # First, let's create a solution to optimise.
    # In this case, we choose a solution model that exists in
    # the BurnMan repository; the Mg-Fe olivine from
    # the Stixrude and Lithgow-Bertelloni dataset
    solution = burnman.minerals.SLB_2011.mg_fe_olivine()
    solution.set_state(1.e5, 300.)

    print('Names of endmembers in the olivine solution:')
    print(solution.endmember_names)
    print('')

    # Fit parameters are provided via a list of lists.
    # The first element of each list is a string that corresponds
    # either to one of the keys in an endmember parameter dictionary,
    # or to an excess property for a binary join in the solution.
    # The next parameters correspond to the indices of the endmembers
    # to which the parameter corresponds.

    # Here, we choose to fit he standard state volume, isothermal
    # bulk modulus and its first derivative for both endmembers.
    # Endmember 0 is forsterite, and Endmember 1 is fayalite.
    # We also choose to fit the excess volume on the binary join.
    fit_params = [['V_0', 0],
                  ['V_0', 1],
                  ['K_0', 0],
                  ['K_0', 1],
                  ['Kprime_0', 0],
                  ['Kprime_0', 1],
                  ['V', 0, 1]]

    # Next, we make some synthetic data
    n_data = 100
    data = []
    data_covariances = []

    # All of the data is volume data, so we declare a single
    # flag. For data of different types, you can also pass a
    # list of flags to the solver.
    flags = 'V'

    # For this example, we add some Gaussian noise
    # to the volumes of olivines on the binary between
    # 0-10 GPa and 300-1300 K
    f_Verror = 1.e-3

    # Choose a specific seed for the random number generator
    # so that this example is reproducible.
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

    # Here, we add one awful data point in the middle of the domain
    # We do this to demonstrate the semi-automatic removal of bad data
    # using extreme value theory.
    solution.set_composition([0.5, 0.5])
    solution.set_state(5.e9, 800.)
    data.append([0.5, 0.5, 5.e9, 800., solution.V + 3.e-7])
    data_covariances.append(np.zeros((5, 5)))
    data_covariances[-1][4, 4] = np.power(solution.V*f_Verror, 2.)

    data = np.array(data)
    data_covariances = np.array(data_covariances)

    # Finally, we choose some initial step sizes for the optimizer.
    delta_params = np.array([1.e-8, 1.e-8, 1.e7, 1.e7, 1.e-1, 1.e-1, 1.e-8])

    # And some bounds. For this example the bounds are not necessary,
    # but if the data is somewhat shaky it can be useful to provide some
    # guidance for the least squares minimizer.
    bounds = np.array([[0, np.inf],
                       [0, np.inf],
                       [0, np.inf],
                       [0, np.inf],
                       [3.5, 6.],
                       [3.5, 6.],
                       [-np.inf, np.inf]])

    param_tolerance = 1.e-5

    # Finally, some post-processing options
    confidence_interval = 0.95
    remove_outliers = True
    good_data_confidence_interval = 0.99
    properties_for_data_comparison_plots = [('V', 1.e6, 'Volume (cm^3/mol)')]

    # The following line fits the parameters to the data we defined above.
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
