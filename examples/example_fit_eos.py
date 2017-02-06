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

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman

if __name__ == "__main__":
    
    print('Least squares equation of state fitting\n')
    print('1) Fitting to PTV data only\n')
    
    # First, let's read in the PVT equation of state data for MgO from Dewaele et al., (2000).
    T, Terr, Pta, P, Perr, V, Verr = np.loadtxt('../burnman/data/input_fitting/PVT_MgO_Dewaele_et_al_2000.dat', unpack=True)
    PTV = np.array([P*1.e9, T, burnman.tools.molar_volume_from_unit_cell_volume(V, 4.)]).T
    nul = 0.*PTV.T[0]

    Pcov = np.power(Perr*1.e9, 2.)
    Tcov = np.power(Terr, 2.)
    Vcov = np.power(burnman.tools.molar_volume_from_unit_cell_volume(Verr, 4.), 2.)
    
    PTV_covariance = np.array([[Pcov, nul, nul], [nul, Tcov, nul], [nul, nul, Vcov]]).T

    
    # Here's where we fit the data
    # We choose to fit to the SLB_2011 equation of state with four parameters (V_0, K_0, K'_0, grueneisen_0)
    per_opt = burnman.minerals.SLB_2011.periclase()
    params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0']
    fitted_eos = burnman.tools.fit_PTV_data(per_opt, params, PTV, PTV_covariance, verbose=False)

    
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
    im = ax.scatter(PTV[:,0]/1.e9, PTV[:,1], c=fitted_eos.weighted_residuals, cmap=plt.cm.RdYlBu, s=50)
    fig.colorbar(im, ax=ax)
    max_abs_r = np.max(np.abs(fitted_eos.weighted_residuals))
    im.set_clim(-max_abs_r, max_abs_r)
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

    print('Clearly a PVT fit of a PVT equation of state can produce some surprising thermodynamic predictions. '
          'This is because PVT alone may not provide strong constraints on all parameters. '
          'In this case, grueneisen_0 and q_0 have very large covariances, which is undesirable.\n\n'
          'Can we improve on this "optimal" PVT fit?\n')
    

    # 2) Now let's be a bit more exotic and fit enthalpy *and* PVT data together

    print('2) Fitting to PTV *and* PT-enthalpy data\n')
    
    # Create the mineral instance that we'll optimise
    per_opt = burnman.minerals.SLB_2011.periclase() 

    
    # Let's load some enthalpy data as extra constraints.
    # We'll use the data that we used in example_fit_data.py:
    TH_data = np.loadtxt('../burnman/data/input_fitting/Victor_Douglas_1963_deltaH_MgO.dat')
    per_opt.set_state(1.e5, 298.15)
    PTH_data = np.array([TH_data[:,0]*0. + 1.e5, TH_data[:,0], TH_data[:,2]*4.184 + per_opt.H]).T
    nul = TH_data[:,0]*0.
    PTH_covariances = np.array([[nul, nul, nul], [nul, TH_data[:,1], nul], [nul, nul, np.power(TH_data[:,2]*4.184*0.0004, 2.)]]).T


    # Just for completeness, let's reload our PVT data:
    T, Terr, Pta, P, Perr, V, Verr = np.loadtxt('../burnman/data/input_fitting/PVT_MgO_Dewaele_et_al_2000.dat', unpack=True)
    PTV_data = np.array([P*1.e9, T, burnman.tools.molar_volume_from_unit_cell_volume(V, 4.)]).T
    nul = 0.*T

    Pcov = np.power(Perr*1.e9, 2.)
    Tcov = np.power(Terr, 2.)
    Vcov = np.power(burnman.tools.molar_volume_from_unit_cell_volume(Verr, 4.), 2.)
    
    PTV_covariances = np.array([[Pcov, nul, nul], [nul, Tcov, nul], [nul, nul, Vcov]]).T

    # Ok, so we've now got two sets of data, PTH and PTV. Let's combine them and create some flags:
    flags = ['H' for P in PTH_data[:,0]]
    flags.extend(['V' for P in PTV_data[:,0]])

    PTp_data = np.concatenate([PTH_data, PTV_data])
    PTp_covariances = np.concatenate([PTH_covariances, PTV_covariances])

    # And now we can fit our data!
    fitted_eos = burnman.tools.fit_PTp_data(mineral = per_opt,
                                            p_flags = flags,
                                            fit_params = ['V_0', 'K_0', 'Kprime_0', 'grueneisen_0', 'q_0', 'F_0'],
                                            PTp = PTp_data,
                                            PTp_covariances = PTp_covariances,
                                            verbose = False)

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
    mask = [i for i, flag in enumerate(flags) if flag=='V']
    im = ax.scatter(PTp_data[:,0][mask]/1.e9, PTp_data[mask][:,1], c=fitted_eos.weighted_residuals[mask], cmap=plt.cm.RdYlBu, s=50)
    fig.colorbar(im, ax=ax)
    max_abs_r = np.max(np.abs(fitted_eos.weighted_residuals[mask]))
    im.set_clim(-max_abs_r, max_abs_r)
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

        # Plot the 95% confidence bands 
        cp_bands = burnman.nonlinear_fitting.orthogonal_distance_confidence_prediction_bands(fitted_eos, PTVs, 0.95, [], 'V')
        plt.plot(cp_bands[0][:,0]/1.e9, cp_bands[0][:,2] * 1.e6, linestyle='--', color='r', label='95% confidence bands')
        plt.plot(cp_bands[1][:,0]/1.e9, cp_bands[1][:,2] * 1.e6, linestyle='--', color='r')
        
        plt.plot(PTVs[:,0] / 1.e9, PTVs[:,2] * 1.e6,
                 label='Optimized fit for periclase at {0:.0f} K'.format(T))
        
    plt.errorbar(PTp_data[:,0][mask] / 1.e9, PTp_data[:,2][mask] * 1.e6,
                 xerr=np.sqrt(PTp_covariances.T[0][0][mask]) / 1.e9,
                 yerr=np.sqrt(PTp_covariances.T[2][2][mask]) * 1.e6,
                 linestyle='None', marker='o', label='Dewaele et al. (2000)')

    plt.scatter(fitted_eos.data_mle[:,0][mask] / 1.e9, fitted_eos.data_mle[:,2][mask] * 1.e6, s=20, label='Maximum likelihood estimates')
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
            cp2_bands = burnman.nonlinear_fitting.confidence_prediction_bands(fitted_eos,
                                                                              burnman.tools.attribute_function(per_opt, material_property),
                                                                              PTVs, 0.95, 'V')
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0] + cp2_bands[1])/2/scaling, label='Best fit at {0:.0f} K'.format(T))
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0])/scaling, linestyle='--', color='r', label='95% confidence band')
            ax.plot(PTVs[:,0]/1.e9, (cp2_bands[1])/scaling, linestyle='--', color='r')
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
        cp2_bands = burnman.nonlinear_fitting.confidence_prediction_bands(fitted_eos,
                                                                          burnman.tools.attribute_function(per_opt, 'gibbs'),
                                                                          PTVs, 0.95, 'V')
        ax.plot(PTVs[:,0]/1.e9, (cp2_bands[0] - cp2_bands[1])/2/scaling, label='95% confidence half width at {0:.0f} K'.format(T))
        plt.ylabel('Gibbs free energy uncertainty (kJ/mol)')
        plt.xlabel('Pressure (GPa)')
            
    plt.legend(loc='lower right')
    plt.show()

    
