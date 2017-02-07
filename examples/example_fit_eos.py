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
    print('1) Fitting to high pressure PTV data only\n')
    
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
    fitted_eos = burnman.tools.fit_PTV_data(per_opt, params, PTV_data, PTV_covariances, verbose=False)

    
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


    sd_limit = 3 # plot discrete colors up to this integer limit 
    cmap = plt.cm.RdYlBu
    cmap.set_under('k')
    cmap.set_over('k')
    bounds = np.linspace(-sd_limit, sd_limit, sd_limit*2+1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    im = ax.scatter(PTV_data[:,0]/1.e9, PTV_data[:,1], c=fitted_eos.weighted_residuals, cmap=cmap, norm=norm, s=50)
    fig.colorbar(im, ax=ax)
    plt.show()
    
    per_SLB = burnman.minerals.SLB_2011.periclase()
    per_HP = burnman.minerals.HP_2011_ds62.per()

    temperatures = np.linspace(200., 2000., 101)
    pressures = np.array([298.15] * len(temperatures))
    plt.plot(temperatures, per_HP.evaluate(['heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='HP')
    plt.plot(temperatures, per_SLB.evaluate(['heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='SLB')

    temperatures = np.linspace(200., 1200., 101)
    plt.plot(temperatures, per_opt.evaluate(['heat_capacity_p'], pressures, temperatures)[0], label='Optimised fit')
    

    plt.legend(loc='lower right')
    plt.xlim(0., temperatures[-1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat capacity (J/K/mol)')
    plt.show()

    print('Fitting an equation of state using only high pressure P-T-V data can produce some surprising '
          '(read incorrect) thermodynamic predictions. We must use additional data to provide strong '
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
    fitted_eos = burnman.tools.fit_PTp_data(mineral = per_opt,
                                            p_flags = flags,
                                            fit_params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0', 'F_0'],
                                            data = PTp_data,
                                            data_covariances = PTp_covariances,
                                            verbose = False)


    
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
    plt.show()
    
    # Print the optimized parameters
    print('Optimized equation of state:')
    burnman.tools.pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print('')
    
    # Create a corner plot of the covariances
    fig=burnman.nonlinear_fitting.corner_plot(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    plt.show()

    # Now let's plot the volume residuals
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    sd_limit = 3 # plot discrete colors up to this integer limit 
    cmap = plt.cm.RdYlBu
    cmap.set_under('k')
    cmap.set_over('k')
    bounds = np.linspace(-sd_limit, sd_limit, sd_limit*2+1)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    
    im = ax.scatter(PTV_data[:,0]/1.e9, PTV_data[:,1], c=fitted_eos.weighted_residuals[V_mask], cmap=cmap, norm=norm, s=50)
    fig.colorbar(im, ax=ax)
    plt.show()
    
    # Plot models
    
    per_SLB = burnman.minerals.SLB_2011.periclase()
    per_HP = burnman.minerals.HP_2011_ds62.per()
    temperatures = np.linspace(200., 2000., 101)
    pressures = np.array([298.15] * len(temperatures))
    plt.plot(temperatures, per_HP.evaluate(['heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='HP')
    plt.plot(temperatures, per_SLB.evaluate(['heat_capacity_p'], pressures, temperatures)[0], linestyle='--', label='SLB')
    plt.plot(temperatures, per_opt.evaluate(['heat_capacity_p'], pressures, temperatures)[0], label='Optimised fit')

    plt.legend(loc='lower right')
    plt.xlim(0., temperatures[-1])
    plt.xlabel('Temperature (K)')
    plt.ylabel('Heat capacity (J/K/mol)')
    plt.show()



    print('Hurrah! That looks much better! Note that not only have the errors on grueneisen_0 and q_0 got smaller, '
          'but they have also moved well-outside the previous 1-s.d. uncertainty bounds. '
          'Plus, now the heat capacities don\'t blow up at high temperature. Cool, huh?!')



    
    
    # Here we plot our equation of state, along with the 95% confidence intervals for the volume
    
    pressures = np.linspace(1.e5, 100.e9, 101)
    for T in [298.15, 2000.]:
        temperatures = pressures*0. + T
        PTVs = np.array([pressures, temperatures, per_opt.evaluate(['V'], pressures, temperatures)[0]]).T

        # Plot the 95% confidence bands on the volumes 
        cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                         x_array=PTVs,
                                                                         confidence_interval=0.95,
                                                                         f=burnman.tools.attribute_function(per_opt, 'V'),
                                                                         flag='V')
        plt.plot(PTVs[:,0]/1.e9, (cp_bands[0])*1.e6, linestyle='--', color='r', label='95% confidence bands')
        plt.plot(PTVs[:,0]/1.e9, (cp_bands[1])*1.e6, linestyle='--', color='r')


        plt.plot(pressures/1.e9, per_opt.evaluate(['V'], pressures, temperatures)[0]* 1.e6, label='Optimised fit at {0:.0f} K'.format(T))
        
    plt.errorbar(PTp_data[:,0][V_mask] / 1.e9, PTp_data[:,2][V_mask] * 1.e6,
                 xerr=np.sqrt(PTp_covariances.T[0][0][V_mask]) / 1.e9,
                 yerr=np.sqrt(PTp_covariances.T[2][2][V_mask]) * 1.e6,
                 linestyle='None', marker='o', label='Dewaele et al. (2000)')

    plt.scatter(fitted_eos.data_mle[:,0][V_mask] / 1.e9, fitted_eos.data_mle[:,2][V_mask] * 1.e6, s=20, label='Maximum likelihood estimates')
    plt.ylabel("Volume (cm^3/mol)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Periclase EoS")
    plt.show()


    # We can also look at the uncertainty in other properties
    # For example, let's look at the uncertainty in P wave velocities, bulk modulus, thermal expansion and thermal pressure
    fig = plt.figure()
    for T in [298.15, 2000.]:
        temperatures = pressures*0. + T
        PTVs = np.array([pressures, temperatures, per_opt.evaluate(['V'], pressures, temperatures)[0]]).T

        for i, (material_property, scaling, name) in enumerate([('p_wave_velocity', 1.e3, 'P wave velocity (km/s)'),
                                                                ('K_T', 1.e9, 'Bulk modulus (GPa)'),
                                                                ('alpha', 1., 'Thermal expansion (/K)'),
                                                                (['alpha', 'K_T'], 1.e6, 'Thermal pressure (MPa/K)')]):
            ax = fig.add_subplot(2, 2, i+1)
            
            # Plot the 95% confidence bands for the various material properties
            cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                             x_array=PTVs,
                                                                             confidence_interval=0.95,
                                                                             f=burnman.tools.attribute_function(per_opt, material_property),
                                                                             flag='V')
            ax.plot(PTVs[:,0]/1.e9, (cp_bands[0] + cp_bands[1])/2/scaling, label='Best fit at {0:.0f} K'.format(T))
            ax.plot(PTVs[:,0]/1.e9, (cp_bands[0])/scaling, linestyle='--', color='r', label='95% confidence band')
            ax.plot(PTVs[:,0]/1.e9, (cp_bands[1])/scaling, linestyle='--', color='r')
            plt.ylabel(name)
            plt.xlabel('Pressure (GPa)')
            
    plt.legend(loc='upper right')
    plt.show()

    
    # Finally, let's look at the effect of uncertainties in the equation of state on gibbs free energy at high pressure
    fig = plt.figure()
    for T in [298.15, 2000.]:
        temperatures = pressures*0. + T
        PTVs = np.array([pressures, temperatures, per_opt.evaluate(['V'], pressures, temperatures)[0]]).T

        ax = fig.add_subplot(1, 1, 1)
        scaling = 1.e3
        
        # Plot the 95% confidence bands for the gibbs free energy
        cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(model=fitted_eos,
                                                                         x_array=PTVs,
                                                                         confidence_interval=0.95,
                                                                         f=burnman.tools.attribute_function(per_opt, 'gibbs'),
                                                                         flag='V')
        ax.plot(PTVs[:,0]/1.e9, (cp_bands[0] - cp_bands[1])/2/scaling, label='95% confidence half width at {0:.0f} K'.format(T))
        plt.ylabel('Gibbs free energy uncertainty (kJ/mol)')
        plt.xlabel('Pressure (GPa)')
            
    plt.legend(loc='lower right')
    plt.show()

    
