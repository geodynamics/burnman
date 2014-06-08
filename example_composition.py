# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows how to create different minerals, how to compute seismic
velocities, and how to compare them to a seismic reference model.

Currently BurnMan can input mineral compositions in 4 different ways: 

1. Two minerals mixed in simple mole fractions. Can be chosen from the BurnMan
   libraries or from user defined minerals (see example_user_input_material)
2. Two minerals mixed in simple mole fractions with user-defined Fe
   partitioning
3. The user can input wt% of each cation (Mg, Fe and Si) and BurnMan will
   calculate Fe partitioning along a P, T profile (see
   example_partition_coef.py)
4. A mixture of three or more  minerals. 

In compositions 2,3 and 4 of the above inputs, BurnMan will mix the mineral
physical paremeters of end member minerals (pure Mg and Fe) of the user's
choice using either volumetric (moduli) or molar averaging (all others) at
room pressure and temperature (see example_user_input_material.py for
information on these parameters) .

To turn a method of mineral creation "on" the first if statement above the
method must be set to True, with all others set to False.

Note: These minerals can include a spin transition in (Mg,Fe)O, see
example_spintransition.py for explanation of how to implement this

requires:
- geotherms
- seismic models
- compute seismic velocities
- composite mineral helpers

teaches:
- creating minerals
- seismic comparison

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
    #INPUT for method
    """ choose 'slb2' (finite-strain 2nd order sheer modulus, stixrude and 
    	lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus, matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shearl modulus, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature
    (your choice in geotherm will not matter in this case))"""
    
    method = 'slb3' 
        
    
    # To compute seismic velocities and other properties, we need to supply
    # burnman with a list of minerals (phases) and their molar abundances. Minerals
    # are classes found in burnman.minerals and are derived from
    # burnman.minerals.material.
    # Here are a few ways to define phases and molar_abundances:
    
    #Example 1: two simple fixed minerals
    if True:
        amount_perovskite = 0.95
        rock = burnman.Composite ( [amount_perovskite, 1.0-amount_perovskite], \
                                       [minerals.SLB_2011.mg_perovskite(), \
                                            minerals.SLB_2011.periclase()] )
    
    #Example 2: specify fixed iron content
    if False:
        amount_perovskite = 0.95
        rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite], \
                                      [minerals.SLB_2011.mg_fe_perovskite(0.2), \
                                           minerals.SLB_2011.ferropericlase(0.2)] )
    
    #Example 3: input weight percentages
    #See comments in example_partition_coef.py for references to 
    #partition coefficent calculation

    if False:
        weight_percents = {'Mg':0.213, 'Fe': 0.08, 'Si':0.27, 'Ca':0., 'Al':0.}
        Kd_0 = .59 #Fig 5 Nakajima et al 2012

        phase_fractions,relative_molar_percent = burnman.\
        							calculate_phase_percents(weight_percents)
        iron_content = lambda p,t: burnman.calculate_partition_coefficient\
        							(p,t,relative_molar_percent,Kd_0)

        rock = burnman.Composite ( [ (minerals.mg_fe_perovskite_pt_dependent\
        							(iron_content,0), phase_fractions['pv'] ),
                                    (minerals.ferropericlase_pt_dependent(iron_content,1),\
                                    phase_fractions['fp'] ) ] )
        
    #Example 4: three materials
    if False:
        rock = burnman.Composite ( [ (minerals.SLB_2011.fe_perovskite(), 0.7),
                                     (minerals.SLB_2011.ferropericlase(0.5), 0.2) ,
                                     (minerals.SLB_2011.stishovite(), 0.1 ) ] )
    
    
    #seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.prem()
    number_of_points = 20 #set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    temperature = burnman.geotherm.brown_shankland(seis_p)
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    rock.debug_print()
    
    moduli_list = burnman.calculate_moduli(rock, seis_p, temperature)
    moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.voigt_reuss_hill())
    mat_vp, mat_vs, mat_vphi = burnman.compute_velocities(moduli)
    mat_K = np.array([m.K for m in moduli])
    mat_G = np.array([m.G for m in moduli])
    mat_rho = np.array([m.rho for m in moduli])
    [rho_err,vphi_err,vs_err]=burnman.compare_chifactor\
    	([mat_vs,mat_vphi,mat_rho],[seis_vs,seis_vphi,seis_rho])
        
    
    # PLOTTING
    
    # plot vs
    plt.subplot(2,2,1)
    plt.plot(seis_p/1.e9,mat_vs/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label='computation')
    plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4,label='reference')
    plt.title("Vs (km/s)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.ylim(5.1,7.6)
    plt.legend(loc='lower right')
    plt.text(40,7.3,"misfit= %3.3f" % vs_err)
    
    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p/1.e9,mat_vphi/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,seis_vphi/1.e3,color='k',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4)
    plt.title("Vphi (km/s)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.ylim(7,12)
    plt.text(40,11.5,"misfit= %3.3f" % vphi_err)
    
    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p/1.e9,mat_rho/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,seis_rho/1.e3,color='k',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4)
    plt.title("density (kg/m^3)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.text(40,4.3,"misfit= %3.3f" % rho_err)
    plt.xlabel("Pressure (GPa)")
    
    
    # plot geotherm
    plt.subplot(2,2,4)
    plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4)
    plt.title("Geotherm (K)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.xlabel("Pressure (GPa)")
    
    plt.savefig("output_figures/example_composition.png")
    plt.show()
