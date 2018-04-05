# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['axes.edgecolor'] = 'black'


# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import burnman
from fitting_functions import read_fitting_file

if __name__ == "__main__":

    """
    First, please create a data file. This file should be in one of the two following formats:
    Type, P, T, property
    Type, P, T, property, Perr, Terr, property_err
    Type, P, T, property, cov_PP, cov_TT, cov_pp, cov_PT, cov_Pp, cov_Tp
    where
    err means standard error, and 
    cov is the covariance matrix of the data observations
    PLEASE REMEMBER THAT cov_PP = Perr^2.

    Type is a string describing the property value *as used in burnman*.
    The property strings could be any of the following:

    helmholtz
    gibbs
    H
    S
    V
    molar_heat_capacity_p
    molar_heat_capacity_v
    p_wave_velocity
    s_wave_velocity
    K_T
    K_S

    Make sure that *all* parameters are in SI units.
    """
    
    # Input file
    filename = 'test.dat'

    # Mineral to optimise (along with tweaks to initial properties if necessary)
    mineral = burnman.minerals.SLB_2011.periclase()
    mineral.set_state(1.e5, 300.)
    mineral.params['F_0'] = mineral.params['F_0'] - mineral.H

    # Fit parameters
    fit_params =  ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0', 'Debye_0', 'F_0']

    # Pressure and temperature sections to plot through the models
    pressures = np.linspace(1.e5, 100.e9, 101)
    temperature_sections = [300., 2000.]
    
    temperatures = np.linspace(300., 2000., 101)
    pressure_sections = [1.e5]

    # Properties to plot which have data in the input file
    properties_for_data_comparison_plots = [('V', 1.e6, 'Volume (cm^3/mol)'),
                                            ('H', 1.e-3, 'Enthalpy (kJ/mol)')]
    
    # Properties to plot along with confidence interval
    properties_for_confidence_plots = [('p_wave_velocity', 1.e-3, 'P wave velocity (km/s)'),
                                       ('K_T', 1.e-9, 'Bulk modulus (GPa)'),
                                       ('alpha', 1., 'Thermal expansion (/K)'),
                                       (['alpha', 'K_T'], 1.e-6, 'Thermal pressure (MPa/K)')]
    confidence_interval = 0.95
    remove_outliers = True
    good_data_confidence_interval = 0.9
    param_tolerance = 1.e-5

    # That's it for user inputs. Now just sit back and watch the plots appear...
    flags, data, data_covariances = read_fitting_file(filename)
    
    list_flags = list(set(flags))
    
    print('Starting to fit user-defined data. Please be patient.')    
    fitted_eos = burnman.eos_fitting.fit_PTp_data(mineral = mineral,
                                                  flags = flags,
                                                  fit_params = fit_params,
                                                  data = data,
                                                  data_covariances = data_covariances,
                                                  param_tolerance = param_tolerance,
                                                  verbose = False)

    # Print the optimized parameters
    print('Optimized equation of state:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('\nParameters:')
    print(fitted_eos.popt)
    print('\nFull covariance matrix:')
    print(fitted_eos.pcov)
    print('\nGoodness of fit:')
    print(fitted_eos.goodness_of_fit)
    print('\n')

    # Create a plot of the residuals
    fig, ax = plt.subplots()
    burnman.nonlinear_fitting.plot_residuals(ax=ax,
                                             weighted_residuals=fitted_eos.weighted_residuals,
                                             flags=fitted_eos.flags)
    plt.show()

    confidence_bound, indices, probabilities = burnman.nonlinear_fitting.extreme_values(fitted_eos.weighted_residuals, good_data_confidence_interval)
    if indices != [] and remove_outliers == True:
        print('Removing {0:d} outliers (at the {1:.1f}% confidence interval) and refitting. Please wait just a little longer.'.format(len(indices), good_data_confidence_interval*100.))
        
        mask = [i for i in range(len(fitted_eos.weighted_residuals)) if i not in indices]
        flags = [flag for i, flag in enumerate(flags) if i not in indices]
        data = data[mask]
        data_covariances = data_covariances[mask]  
        fitted_eos = burnman.eos_fitting.fit_PTp_data(mineral = mineral,
                                                      flags = flags,
                                                      fit_params = fit_params,
                                                      data = data,
                                                      data_covariances = data_covariances,
                                                      param_tolerance = param_tolerance,
                                                      verbose = False)

    # Print the optimized parameters
    print('Optimized equation of state:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('\nParameters:')
    print(fitted_eos.popt)
    print('\nFull covariance matrix:')
    print(fitted_eos.pcov)
    print('\nGoodness of fit:')
    print(fitted_eos.goodness_of_fit)
    print('\n')
        

    # Create a plot of the residuals
    fig, ax = plt.subplots()
    burnman.nonlinear_fitting.plot_residuals(ax=ax,
                                             weighted_residuals=fitted_eos.weighted_residuals,
                                             flags=fitted_eos.flags)
    plt.show()
    
    # Create a corner plot of the covariances
    fig, ax_array = burnman.nonlinear_fitting.corner_plot(popt=fitted_eos.popt,
                                                          pcov=fitted_eos.pcov,
                                                          param_names=fitted_eos.fit_params)
    plt.show()

    # Create plots for the weighted residuals of each type of measurement
    for i, (material_property, scaling, name) in enumerate(properties_for_data_comparison_plots):
        fig, ax = plt.subplots()
        burnman.nonlinear_fitting.weighted_residual_plot(ax=ax,
                                                         model=fitted_eos,
                                                         flag=material_property,
                                                         sd_limit=3,
                                                         cmap=plt.cm.RdYlBu,
                                                         plot_axes=[0, 1],
                                                         scale_axes=[1.e-9, 1.])
        ax.set_title('Weighted residual plot for {0:s}'.format(name))
        ax.set_xlabel('Pressure (GPa)')
        ax.set_ylabel('Temperature (K)')
        plt.show()


        flag_mask = [i for i, flag in enumerate(flags) if flag==material_property]
    
        if temperature_sections != []:
            for T in temperature_sections:
                PTVs = np.array([pressures, [T]*len(pressures), mineral.evaluate(['V'], pressures, [T]*len(pressures))[0]]).T

                # Plot confidence bands on the volumes 
                cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                                 x_array=PTVs,
                                                                                 confidence_interval=confidence_interval,
                                                                                 f=burnman.tools.attribute_function(mineral, material_property),
                                                                                 flag='V')
            
                plt.plot(PTVs[:,0] / 1.e9, (cp_bands[0] + cp_bands[1])/2.*scaling, label='Optimised fit at {0:.0f} K'.format(T))
                plt.plot(PTVs[:,0] / 1.e9, (cp_bands[0])*scaling, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
                plt.plot(PTVs[:,0] / 1.e9, (cp_bands[1])*scaling, linestyle='--', color='r')
        
            plt.errorbar(fitted_eos.data[:,0][flag_mask] / 1.e9, fitted_eos.data[:,2][flag_mask]*scaling,
                         xerr=np.sqrt(fitted_eos.data_covariances.T[0][0][flag_mask]) / 1.e9,
                         yerr=np.sqrt(fitted_eos.data_covariances.T[2][2][flag_mask])*scaling,
                         linestyle='None', marker='o', label='Data')
            
            plt.plot(fitted_eos.data_mle[:,0][flag_mask] / 1.e9, fitted_eos.data_mle[:,2][flag_mask]*scaling, marker='o', markersize=2, color='k', linestyle='None', label='Maximum likelihood estimates')
            plt.ylabel('{0:s}'.format(name))
            plt.xlabel('Pressure (GPa)')
            plt.legend(loc='upper right')
            plt.title('Data comparison for fitted equation of state as a function of pressure')
            plt.show()


        if pressure_sections != []:
            for P in pressure_sections:
                PTVs = np.array([[P]*len(temperatures), temperatures, mineral.evaluate(['V'], [P]*len(temperatures), temperatures)[0]]).T

                # Plot confidence bands on the volumes 
                cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                                 x_array=PTVs,
                                                                                 confidence_interval=confidence_interval,
                                                                                 f=burnman.tools.attribute_function(mineral, material_property),
                                                                                 flag='V')
                
                plt.plot(PTVs[:,1], (cp_bands[0] + cp_bands[1])/2.*scaling, label='Optimised fit at {0:.0f} GPa'.format(P/1.e9))
                plt.plot(PTVs[:,1], (cp_bands[0])*scaling, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
                plt.plot(PTVs[:,1], (cp_bands[1])*scaling, linestyle='--', color='r')
        
            plt.errorbar(fitted_eos.data[:,1][flag_mask], fitted_eos.data[:,2][flag_mask]*scaling,
                         xerr=np.sqrt(fitted_eos.data_covariances.T[1][1][flag_mask]),
                         yerr=np.sqrt(fitted_eos.data_covariances.T[2][2][flag_mask])*scaling,
                         linestyle='None', marker='o', label='Data')
            
            plt.plot(fitted_eos.data_mle[:,1][flag_mask], fitted_eos.data_mle[:,2][flag_mask]*scaling, marker='o', markersize=2, color='k', linestyle='None', label='Maximum likelihood estimates')
            plt.ylabel('{0:s}'.format(name))
            plt.xlabel('Temperature (K)')
            plt.legend(loc='upper right')
            plt.title('Data comparison for fitted equation of state as a function of temperature')
            plt.show()

        


    # We can also look at the uncertainty in other properties
    # For example, let's look at the uncertainty in P wave velocities, bulk modulus, thermal expansion and thermal pressure
    def closest_factors(n):
        d = np.int(np.floor(np.sqrt(n)))
        for i in reversed(range(1, d+1)):
            if (n % i) == 0:
                return i, int(n/i)

    nj, ni = closest_factors(len(properties_for_confidence_plots))

    if temperature_sections != []:
        fig = plt.figure()
        for T in temperature_sections:
            PTVs = np.array([pressures, [T]*len(pressures), mineral.evaluate(['V'], pressures, [T]*len(pressures))[0]]).T

            for i, (material_property, scaling, name) in enumerate(properties_for_confidence_plots):
                ax = fig.add_subplot(ni, nj, i+1)
            
                # Plot the confidence bands for the various material properties
                cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                                 x_array=PTVs,
                                                                                 confidence_interval=confidence_interval,
                                                                                 f=burnman.tools.attribute_function(mineral, material_property),
                                                                                 flag='V')
                ax.plot(PTVs[:,0]/1.e9, (cp_bands[0] + cp_bands[1])/2*scaling, label='Best fit at {0:.0f} K'.format(T))
                ax.plot(PTVs[:,0]/1.e9, (cp_bands[0])*scaling, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
                ax.plot(PTVs[:,0]/1.e9, (cp_bands[1])*scaling, linestyle='--', color='r')
                
                plt.ylabel(name)
                plt.xlabel('Pressure (GPa)')
        
        plt.legend(loc='upper right')
        plt.show()



    if pressure_sections != []:
        fig = plt.figure()
        for P in pressure_sections:
            PTVs = np.array([[P]*len(temperatures), temperatures, mineral.evaluate(['V'], [P]*len(temperatures), temperatures)[0]]).T

            for i, (material_property, scaling, name) in enumerate(properties_for_confidence_plots):
                ax = fig.add_subplot(ni,nj, i+1)
            
                # Plot the confidence bands for the various material properties
                cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                                 x_array=PTVs,
                                                                                 confidence_interval=confidence_interval,
                                                                                 f=burnman.tools.attribute_function(mineral, material_property),
                                                                                 flag='V')
                ax.plot(PTVs[:,1], (cp_bands[0] + cp_bands[1])/2*scaling, label='Best fit at {0:.0f} GPa'.format(P/1.e9))
                ax.plot(PTVs[:,1], (cp_bands[0])*scaling, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
                ax.plot(PTVs[:,1], (cp_bands[1])*scaling, linestyle='--', color='r')
                plt.ylabel(name)
                plt.xlabel('Temperature (K)')

        plt.legend(loc='upper right')
        plt.show()

