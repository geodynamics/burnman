# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_fit_eos
----------------

This example demonstrates BurnMan's functionality to fit data to
an EoS of the user's choice. 

The first example deals with simple PVT fitting. 
The second example illustrates how powerful it can be to 
provide non-PVT constraints to the same fitting problem.

teaches:
- least squares fitting

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman

if __name__ == "__main__":
    
    print('Least squares equation of state fitting\n')

    print('1) Fitting to room temperature PV data\n')
    # Now let's fit some EoS data
    # Let's just take a bit of data from Andrault et al. (2003) on stishovite
    PV = np.array([[0.0001, 46.5126, 0.0061],
                   [1.168, 46.3429, 0.0053],
                   [2.299, 46.1756, 0.0043],
                   [3.137, 46.0550, 0.0051],
                   [4.252, 45.8969, 0.0045],
                   [5.037, 45.7902, 0.0053],
                   [5.851, 45.6721, 0.0038],
                   [6.613, 45.5715, 0.0050],
                   [7.504, 45.4536, 0.0041],
                   [8.264, 45.3609, 0.0056],
                   [9.635, 45.1885, 0.0042],
                   [11.69, 44.947, 0.002],
                   [17.67, 44.264, 0.002],
                   [22.38, 43.776, 0.003],
                   [29.38, 43.073, 0.009],
                   [37.71, 42.278, 0.008],
                   [46.03, 41.544, 0.017],
                   [52.73, 40.999, 0.009],
                   [26.32, 43.164, 0.006],
                   [30.98, 42.772, 0.005],
                   [34.21, 42.407, 0.003],
                   [38.45, 42.093, 0.004],
                   [43.37, 41.610, 0.004],
                   [47.49, 41.280, 0.007]])
    
    PTV = np.array([PV.T[0]*1.e9,
                    PV.T[0]*0. + 298.15,
                    burnman.tools.molar_volume_from_unit_cell_volume(PV.T[1], 2.)]).T

    nul = 0.*PTV.T[0]
    PTV_covariances = np.array([[0.03*PTV.T[0], nul, nul], [nul, nul, nul], [nul, nul, burnman.tools.molar_volume_from_unit_cell_volume(PV.T[2], 2.)]]).T

    # Here's where we fit the data
    # The mineral parameters are automatically updated during fitting
    stv = burnman.minerals.HP_2011_ds62.stv()
    params = ['V_0', 'K_0', 'Kprime_0']
    fitted_eos = burnman.eos_fitting.fit_PTV_data(stv, params, PTV, PTV_covariances, verbose=False)

    # Print the optimized parameters
    print('Optimized equation of state for stishovite:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('')
    
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, params)
    plt.show()
        
    # Finally, let's plot our equation of state
    T = 298.15
    pressures = np.linspace(1.e5, 60.e9, 101)
    volumes = np.empty_like(pressures)

    PTVs = np.empty((len(pressures), 3))
    for i, P in enumerate(pressures):
        stv.set_state(P, T)
        PTVs[i] = [P, T, stv.V]

    # Plot the 95% confidence and prediction bands 
    cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model = fitted_eos,
                                                                     x_array = PTVs,
                                                                     confidence_interval = 0.95,
                                                                     f=burnman.tools.attribute_function(stv, 'V'),
                                                                     flag='V')
    
    plt.plot(PTVs[:,0]/1.e9, cp_bands[0] * 1.e6, linestyle='--', color='r', label='95% confidence bands')
    plt.plot(PTVs[:,0]/1.e9, cp_bands[1] * 1.e6, linestyle='--', color='r')
    plt.plot(PTVs[:,0]/1.e9, cp_bands[2] * 1.e6, linestyle='--', color='b', label='95% prediction bands')
    plt.plot(PTVs[:,0]/1.e9, cp_bands[3] * 1.e6, linestyle='--', color='b')
        
    plt.plot(PTVs[:,0] / 1.e9, PTVs[:,2] * 1.e6,
             label='Optimized fit for stishovite')
    plt.errorbar(PTV[:,0] / 1.e9, PTV[:,2] * 1.e6,
                 xerr=PTV_covariances.T[0][0] / 1.e9,
                 yerr=PTV_covariances.T[2][2] * 1.e6,
                 linestyle='None', marker='o', label='Andrault et al. (2003)')

    plt.ylabel("Volume (cm^3/mol)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Stishovite EoS (room temperature)")
    plt.show()

    
    # Plot the 95% confidence and prediction bands for the bulk modulus
    cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model = fitted_eos,
                                                                     x_array = PTVs,
                                                                     confidence_interval = 0.95,
                                                                     f=burnman.tools.attribute_function(stv, 'K_T'),
                                                                     flag='V')
    plt.plot(PTVs[:,0]/1.e9, (cp_bands[0] + cp_bands[1])/2.e9, color='b', label='Best fit')
    plt.plot(PTVs[:,0]/1.e9, (cp_bands[0])/1.e9, linestyle='--', color='r', label='95% confidence band')
    plt.plot(PTVs[:,0]/1.e9, (cp_bands[1])/1.e9, linestyle='--', color='r')
    plt.ylabel("Bulk modulus (GPa)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Stishovite EoS; uncertainty in bulk modulus (room temperature)")
    plt.show()
    
    print('2) Fitting to high pressure PTV data only\n')
    
    # First, let's read in the PVT equation of state data for MgO from Dewaele et al., (2000).
    D2000_T, D2000_Terr, D2000_Pta, D2000_P, D2000_Perr, D2000_V, D2000_Verr = np.loadtxt('../burnman/data/input_fitting/PVT_MgO_Dewaele_et_al_2000.dat', unpack=True)
    D2000_P = D2000_P * 1.e9
    D2000_Perr = D2000_Perr * 1.e9
    D2000_V = burnman.tools.molar_volume_from_unit_cell_volume(D2000_V, 4.)
    D2000_Verr = burnman.tools.molar_volume_from_unit_cell_volume(D2000_Verr, 4.)


    PTV_data = np.array([D2000_P, D2000_T, D2000_V]).T
    nul = 0.*PTV_data[:,0]
    
    Pcov = np.power(D2000_Perr, 2.)
    Tcov = np.power(D2000_Terr, 2.)
    Vcov = np.power(D2000_Verr, 2.)
    
    PTV_covariances = np.array([[Pcov, nul, nul], [nul, Tcov, nul], [nul, nul, Vcov]]).T

    
    # Here's where we fit the data
    # We choose to fit to the SLB_2011 equation of state with four parameters (V_0, K_0, K'_0, grueneisen_0)
    per_opt = burnman.minerals.SLB_2011.periclase()
    params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0']
    param_tolerance = 1.e-3 # fairly low resolution fitting (corresponding to 0.1% variation)
    fitted_eos = burnman.eos_fitting.fit_PTV_data(mineral = per_opt,
                                                  fit_params = params,
                                                  data = PTV_data,
                                                  data_covariances = PTV_covariances, 
                                                  param_tolerance = param_tolerance,
                                                  verbose = False)

    
    # We're done! That wasn't too painful, was it?!

    
    # Print the optimized parameters
    print('Optimized equation of state:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('')
    
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    plt.show()

    # Now let's plot the residuals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    burnman.nonlinear_fitting.weighted_residual_plot(ax=ax, model=fitted_eos, flag='V', sd_limit=3,
                                                     cmap=plt.cm.RdYlBu, plot_axes=[0, 1],
                                                     scale_axes=[1.e-9, 1.])
    ax.set_title('Weighted residual plot for volumes')
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Temperature (K)')
    plt.show()

    # Let's also plot a comparison of heat capacities
    per_SLB = burnman.minerals.SLB_2011.periclase()
    per_HP = burnman.minerals.HP_2011_ds62.per()

    temperatures = np.linspace(200., 2000., 101)
    pressures = np.array([298.15] * len(temperatures))
    plt.plot(temperatures, per_HP.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='HP')
    plt.plot(temperatures, per_SLB.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='SLB')

    temperatures = np.linspace(200., 1200., 101)
    plt.plot(temperatures, per_opt.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], label='Optimised fit')
    
    plt.legend(loc='lower right')
    plt.xlim(0., temperatures[-1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat capacity (J/K/mol)')
    plt.show()

    print('Fitting an equation of state using only high pressure P-T-V data can produce some surprising '
          '(incorrect) thermodynamic predictions. We must use additional data to provide strong '
          'constraints on all parameters. In this case, we would like to improve our constraints on '
          'grueneisen_0 and q_0.\n\n'
          'How should we improve on this "optimal" PVT fit?\n')
    

    # 2) Now let's be a bit more exotic and fit enthalpy *and* PVT data together

    print('2) Fitting to PTV *and* PT-enthalpy data\n')
    
    # Create the mineral instance that we'll optimise
    per_opt = burnman.minerals.SLB_2011.periclase() 

    # Let's use some low pressure volume data to supplement the high pressure data:
    H1976_data = np.loadtxt('../burnman/data/input_fitting/Hazen_1976_TV_periclase.dat')
    H1976_T = H1976_data[:,0] + 273.15
    H1976_Terr = np.abs(H1976_data[:,0])*0.001 # assume 0.1% temperature error
    H1976_P = np.array([1.e5] * len(H1976_data[:,0]))
    H1976_Perr = np.array([0.] * len(H1976_data[:,0]))
    H1976_V = burnman.tools.molar_volume_from_unit_cell_volume(np.power(H1976_data[:,1], 3.), 4.)
    H1976_Verr = burnman.tools.molar_volume_from_unit_cell_volume(3.*np.power(H1976_data[:,1], 2.)*H1976_data[:,2], 4.)


    PTV_data = np.array([np.concatenate([D2000_P, H1976_P]),
                         np.concatenate([D2000_T, H1976_T]),
                         np.concatenate([D2000_V, H1976_V])]).T
    nul = 0.*PTV_data[:,0]
    
    Pcov = np.power(np.concatenate([D2000_Perr, H1976_Perr]), 2.)
    Tcov = np.power(np.concatenate([D2000_Terr, H1976_Terr]), 2.)
    Vcov = np.power(np.concatenate([D2000_Verr, H1976_Verr]), 2.)
    
    PTV_covariances = np.array([[Pcov, nul, nul], [nul, Tcov, nul], [nul, nul, Vcov]]).T

    
    # Let's load some enthalpy data as extra constraints.
    # We'll use the data that we used in example_fit_data.py:
    TH_data = np.loadtxt('../burnman/data/input_fitting/Victor_Douglas_1963_deltaH_MgO.dat')
    per_opt.set_state(1.e5, 300.)
    per_opt.params['F_0'] = per_opt.params['F_0'] - per_opt.H # necessary to fit the enthalpy relative to 298.15
    
    PTH_data = np.array([TH_data[:,0]*0. + 1.e5, TH_data[:,0], TH_data[:,2]*4.184]).T
    nul = TH_data[:,0]*0.
    PTH_covariances = np.array([[nul, nul, nul], [nul, TH_data[:,1], nul], [nul, nul, np.power(TH_data[:,2]*4.184*0.0004, 2.)]]).T

    # Ok, so we've now loaded the new enthalpy data. Let's combine it with the volume data we loaded earlier and create some flags:
    flags = ['H' for P in PTH_data[:,0]]
    flags.extend(['V' for P in PTV_data[:,0]])

    PTp_data = np.concatenate([PTH_data, PTV_data])
    PTp_covariances = np.concatenate([PTH_covariances, PTV_covariances])

    V_mask = [i for i, flag in enumerate(flags) if flag == 'V']
    
    # And now we can fit our data!
    fit_params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0', 'Debye_0', 'F_0']
    fitted_eos = burnman.eos_fitting.fit_PTp_data(mineral = per_opt,
                                                  flags = flags,
                                                  fit_params = fit_params,
                                                  data = PTp_data,
                                                  data_covariances = PTp_covariances,
                                                  param_tolerance = param_tolerance,
                                                  verbose = False)


    # Print the optimized parameters
    print('Optimized equation of state:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('')


    
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    plt.show()

    # And a plot of the revised heat capacities
    per_SLB = burnman.minerals.SLB_2011.periclase()
    per_HP = burnman.minerals.HP_2011_ds62.per()
    temperatures = np.linspace(200., 2000., 101)
    pressures = np.array([298.15] * len(temperatures))
    plt.plot(temperatures, per_HP.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='HP')
    plt.plot(temperatures, per_SLB.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='SLB')
    plt.plot(temperatures, per_opt.evaluate(['molar_heat_capacity_p'], pressures, temperatures)[0], label='Optimised fit')

    plt.legend(loc='lower right')
    plt.xlim(0., temperatures[-1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat capacity (J/K/mol)')
    plt.show()

    print('Hmmm, looks promising. The heat capacities don\'t blow up any more. '
          'Let\'s check the new weighted residuals...\n')
    
    # Plot the new residuals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    burnman.nonlinear_fitting.weighted_residual_plot(ax=ax, model=fitted_eos, flag='V', sd_limit=3,
                                                     cmap=plt.cm.RdYlBu, plot_axes=[0, 1],
                                                     scale_axes=[1.e-9, 1.])
    ax.set_title('Weighted residual plot for volumes')
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Temperature (K)')
    plt.show()

    
    # Create a plot of the residuals
    fig, ax = plt.subplots()
    burnman.nonlinear_fitting.plot_residuals(ax=ax,
                                             weighted_residuals=fitted_eos.weighted_residuals,
                                             flags=fitted_eos.flags)
    plt.show()

    print('Eugh. There are two volume points with very large weighted residuals. '
          'You might be interested in asking what the chances are that those extreme '
          'residuals belong in our data set.\n'
          'Here we use extreme value theory to calculate '
          'how many of our data have residuals that are greater than we would expect for a '
          'data set of this size at a given confidence interval. We assume that the errors '
          'provided by the user are approximately correct '
          '(i.e. the variance of weighted residuals should be similar to 1).\n')

    # Let's pick a confidence interval of 90%. Essentially we're now calculating the
    # number of standard deviations from the mean that would bracket an *entire* data set of
    # this size 90% of the time.
    good_data_confidence_interval = 0.9
    confidence_bound, indices, probabilities = burnman.nonlinear_fitting.extreme_values(fitted_eos.weighted_residuals, good_data_confidence_interval)

    print('There are {0:d} outliers (at the {1:.1f}% confidence interval). Their indices and probabilities are:'.format(len(indices), good_data_confidence_interval*100.))
    for i, idx in enumerate(indices):
        print('[{0:d}]: {1:.2f}% ({2:.1f} s.d. from the model)'.format(idx, probabilities[i]*100., np.abs(fitted_eos.weighted_residuals[idx])))
    print('')
    print('As we expected, those two data points are way outside what we should expect. '
          'Least squares fitting is very sensitive to data with large residuals, '
          'so we should always consider adjusting uncertainties or removing suspicious data.\n'
          'Here, the offending data points are so far outside reasonable values that we remove them altogether.')

    mask = [i for i in range(len(fitted_eos.weighted_residuals)) if i not in indices]
    flags = [flag for i, flag in enumerate(flags) if i not in indices]
    PTp_data = PTp_data[mask]
    PTp_covariances = PTp_covariances[mask]  
    fitted_eos = burnman.eos_fitting.fit_PTp_data(mineral = per_opt,
                                                  flags = flags,
                                                  fit_params = fit_params,
                                                  data = PTp_data,
                                                  data_covariances = PTp_covariances,
                                                  param_tolerance = 1.e-5, # higher resolution fitting now that we removed those outliers
                                                  verbose = False)

    # Print the optimized parameters
    print('Optimized equation of state (outliers removed):')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('\nGoodness of fit:')
    print(fitted_eos.goodness_of_fit)
    print('\n')

    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    plt.show()

    # Now let's plot the new residuals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    burnman.nonlinear_fitting.weighted_residual_plot(ax=ax, model=fitted_eos, flag='V', sd_limit=3,
                                                     cmap=plt.cm.RdYlBu, plot_axes=[0, 1],
                                                     scale_axes=[1.e-9, 1.])
    ax.set_title('Weighted residual plot for volumes')
    ax.set_xlabel('Pressure (GPa)')
    ax.set_ylabel('Temperature (K)')
    plt.show()

    # Create a plot of the residuals
    fig, ax = plt.subplots()
    burnman.nonlinear_fitting.plot_residuals(ax=ax,
                                             weighted_residuals=fitted_eos.weighted_residuals,
                                             flags=fitted_eos.flags)
    plt.show()
    
    print('Hurrah! That looks much better! The patches of positive and negative weighted residuals '
          'in P-T space have completely disappeared, and the weighted residuals are now distributed '
          'more evenly about zero. The uncertainties on all the parameters have got much smaller, '
          'and several have moved outside the previous 1-s.d. uncertainty bounds. '
          'Cool, huh?! Let\'s finish this example by looking at some pretty plots '
          'characterising our optimised equation of state.')
    
    # A plot of the volumes at 1 bar
    temperatures = np.linspace(1., 3000., 101)
    pressures = np.array([1.e5] * 101)
    plt.plot(temperatures, per_HP.evaluate(['V'], pressures, temperatures)[0]*1.e6, linestyle='--', label='HP')
    plt.plot(temperatures, per_SLB.evaluate(['V'], pressures, temperatures)[0]*1.e6, linestyle='--', label='SLB')
    plt.plot(temperatures, per_opt.evaluate(['V'], pressures, temperatures)[0]*1.e6, label='Optimised fit')
    
    DS1997_data = np.loadtxt('../burnman/data/input_fitting/Dubrovinsky_Saxena_1997_TV_periclase.dat')
    plt.errorbar(DS1997_data[:,0], DS1997_data[:,1], yerr=DS1997_data[:,2], linestyle='None', label='Dubrovinsky and Saxena (1997)')
    plt.errorbar(H1976_T, H1976_V*1.e6, yerr=H1976_Verr*1.e6,
                 marker='o', markersize=5, linestyle='None', label='Hazen (1976)')

    plt.legend(loc='upper left')
    plt.xlim(0., temperatures[-1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Volume (cm^3/mol)')
    plt.title('Periclase volumes at 1 bar')
    plt.show()
    
    # Here we plot our equation of state, along with the 95% confidence intervals for the volume
    temperature_sections = [298.15, 2000.]
    confidence_interval = 0.95
    flag_mask = [i for i, flag in enumerate(fitted_eos.flags) if flag=='V']
    
    pressures = np.linspace(1.e5, 100.e9, 101)
    for T in temperature_sections:
        PTVs = np.array([pressures, [T]*len(pressures), per_opt.evaluate(['V'], pressures, [T]*len(pressures))[0]]).T
        
        # Plot confidence bands on the volumes 
        cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                         x_array=PTVs,
                                                                         confidence_interval=confidence_interval,
                                                                         f=burnman.tools.attribute_function(per_opt, 'V'),
                                                                         flag='V')
        
        plt.plot(PTVs[:,0] / 1.e9, (cp_bands[0] + cp_bands[1])/2.*1.e6, label='Optimised fit at {0:.0f} K'.format(T))
        plt.plot(PTVs[:,0] / 1.e9, (cp_bands[0])*1.e6, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
        plt.plot(PTVs[:,0] / 1.e9, (cp_bands[1])*1.e6, linestyle='--', color='r')
        
        plt.errorbar(fitted_eos.data[:,0][flag_mask] / 1.e9, fitted_eos.data[:,2][flag_mask]*1.e6,
                     xerr=np.sqrt(fitted_eos.data_covariances.T[0][0][flag_mask]) / 1.e9,
                     yerr=np.sqrt(fitted_eos.data_covariances.T[2][2][flag_mask])*1.e6,
                     linestyle='None', marker='o', label='Data')
        
    plt.plot(fitted_eos.data_mle[:,0][flag_mask] / 1.e9, fitted_eos.data_mle[:,2][flag_mask]*1.e6, marker='o', markersize=2, color='k', linestyle='None', label='Maximum likelihood estimates')
    plt.ylabel('Volume (cm^3/mol)')
    plt.xlabel('Pressure (GPa)')
    plt.legend(loc='upper right')
    plt.title('Data comparison for fitted equation of state as a function of pressure')
    plt.show()

    # We can also look at the uncertainty in other properties
    # For example, let's look at the uncertainty in P wave velocities, bulk modulus, thermal expansion and thermal pressure
    properties_for_confidence_plots = [('p_wave_velocity', 1.e-3, 'P wave velocity (km/s)'),
                                       ('K_T', 1.e-9, 'Bulk modulus (GPa)'),
                                       ('alpha', 1., 'Thermal expansion (/K)'),
                                       (['alpha', 'K_T'], 1.e-6, 'Thermal pressure (MPa/K)')]
    
    def closest_factors(n):
        d = np.int(np.floor(np.sqrt(n)))
        for i in reversed(range(1, d+1)):
            if (n % i) == 0:
                return i, int(n/i)
            
    fig = plt.figure()
    nj, ni = closest_factors(len(properties_for_confidence_plots))
    ax = [fig.add_subplot(ni, nj, i+1) for i in range(len(properties_for_confidence_plots))]
    
    for T in temperature_sections:
        PTVs = np.array([pressures, [T]*len(pressures), per_opt.evaluate(['V'], pressures, [T]*len(pressures))[0]]).T
        
        for i, (material_property, scaling, name) in enumerate(properties_for_confidence_plots):
            
            # Plot the confidence bands for the various material properties
            cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                             x_array=PTVs,
                                                                             confidence_interval=confidence_interval,
                                                                             f=burnman.tools.attribute_function(per_opt, material_property),
                                                                             flag='V')
            ax[i].plot(PTVs[:,0]/1.e9, (cp_bands[0] + cp_bands[1])/2*scaling, label='Best fit at {0:.0f} K'.format(T))
            ax[i].plot(PTVs[:,0]/1.e9, (cp_bands[0])*scaling, linestyle='--', color='r', label='{0:.1f}% confidence bands'.format(confidence_interval*100))
            ax[i].plot(PTVs[:,0]/1.e9, (cp_bands[1])*scaling, linestyle='--', color='r')
            
            ax[i].set_xlabel('Pressure (GPa)')
            ax[i].set_ylabel(name)
            
    ax[i].legend(loc='upper right')
    plt.show()

    
    # Finally, let's look at the effect of uncertainties in the equation of state on gibbs free energy at high pressure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    for T in [298.15, 2000.]:
        temperatures = pressures*0. + T
        PTVs = np.array([pressures, temperatures, per_opt.evaluate(['V'], pressures, temperatures)[0]]).T
        scaling = 1.e3
        
        # Plot the confidence bands for the gibbs free energy
        cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                         x_array=PTVs,
                                                                         confidence_interval=confidence_interval,
                                                                         f=burnman.tools.attribute_function(per_opt, 'gibbs'),
                                                                         flag='V')
        ax.plot(PTVs[:,0]/1.e9, (cp_bands[0] - cp_bands[1])/2/scaling, label='95% confidence half width at {0:.0f} K'.format(T))
          
    ax.set_xlabel('Pressure (GPa)')  
    ax.set_ylabel('Gibbs free energy uncertainty (kJ/mol)')
    ax.legend(loc='lower right')
    plt.show()

    
