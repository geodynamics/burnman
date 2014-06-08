# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows you how to create two materials from wt% determines the
optimum mixing between the two to match the seismic model of your choice.
Currently it compares two end member meteorite groups among the chondrites:
carbonaceous and enstatite. Velocities are calculated for each set of minerals
and plotted for comparison.

requires:
- geotherms
- seismic models
- compute seismic velocities
- creating minerals

teaches:
- weight percent materials

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals


if __name__ == "__main__":    

    ###Input Model 1
    #INPUT for method
    """ choose 'slb2' (finite-strain 2nd order shear modulus,       
        stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus,
        matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shear modulus,
        matas et al. 2007)
    or 'bm2' (birch-murnaghan 2nd order, if you choose to ignore temperature
       (your choice in geotherm will not matter in this case))
    or 'bm3' (birch-murnaghan 3rd order, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case))"""    
    method = 'slb3' 
    
    
    #Carbonaceous chondrites From Mcdonough 2003    
    weight_percents_pyro = {'Mg':0.228, 'Fe': 0.0626, 'Si':0.21, 'Ca':0., 'Al':0.} 
    phase_fractions_pyro,relative_molar_percent_pyro = \
    burnman.calculate_phase_percents(weight_percents_pyro)
    #input initial distribution coefficent. See example_partition_coefficient.py 
    #for explanation
    Kd_0 = .5
    iron_content = lambda p,t: burnman.calculate_partition_coefficient\
    (p,t,relative_molar_percent_pyro, Kd_0)
    pyrolite = burnman.Composite( [phase_fractions_pyro['pv'], phase_fractions_pyro['fp']], \
                                      [minerals.SLB_2005.mg_fe_perovskite_pt_dependent(iron_content,0),\
                                           minerals.SLB_2005.ferropericlase_pt_dependent(iron_content,1)] )
    
    #input pressure range for first model. 
    seis_p_1 = np.arange(25e9, 125e9, 5e9)
    
    #input your geotherm. 
    
    temperature_1 = burnman.geotherm.brown_shankland(seis_p_1)
    
    ##Now onto the second model parameters
    
    ##input second method
    
    #Input composition of model enstatite chondrites. See example_composition for ways of inputting composition.
        
    #Enstatite chondrites: Javoy 2009 Table 6 PLoM
    weight_percents_enst = {'Mg':0.213, 'Fe': 0.0721, 'Si':0.242, 'Ca':0., 'Al':0.} 
    phase_fractions_enst,relative_molar_percent_enst = burnman.\
    calculate_phase_percents(weight_percents_enst)
    iron_content = lambda p,t: burnman.calculate_partition_coefficient\
    (p,t,relative_molar_percent_enst, Kd_0)
    enstatite = burnman.Composite ([phase_fractions_enst['pv'], phase_fractions_enst['fp']], \
                                       [minerals.SLB_2005.mg_fe_perovskite_pt_dependent(iron_content,0), \
                                            minerals.SLB_2005.ferropericlase_pt_dependent(iron_content,1)] )
    
    
    #input second pressure range. Same as the first for comparison
    seis_p_2 = seis_p_1
    
    #input your geotherm. 
    
    temperature_2 = burnman.geotherm.brown_shankland(seis_p_2)
    
    #Now we'll calculate the models. 
    
    pyrolite.set_method(method)
    
    
    
    mat_rho_pyro, mat_vp_pyro, mat_vs_pyro, mat_vphi_pyro, mat_K_pyro, mat_G_pyro = \
        burnman.velocities_from_rock(pyrolite, seis_p_1, temperature_1, \
        burnman.averaging_schemes.VoigtReussHill())
    
    print "Calculations are done for:"
    pyrolite.debug_print()
    
    enstatite.set_method(method)
    
    
    mat_rho_enst, mat_vp_enst, mat_vs_enst, mat_vphi_enst, mat_K_enst, mat_G_enst = \
        burnman.velocities_from_rock(enstatite, seis_p_2, temperature_2, \
        burnman.averaging_schemes.VoigtReussHill())

    print "Calculations are done for:"
    enstatite.debug_print()

    
    ##let's create PREM for reference
    s=burnman.seismic.PREM()
    depths = map(s.depth, seis_p_1) 
    pressures, rho_prem, vp_prem, vs_prem, v_phi_prem = s.evaluate_all_at(depths)
    
    
    ##Now let's plot the comparison.
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,mat_vs_pyro/1.e3,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4)
    plt.plot(seis_p_2/1.e9,mat_vs_enst/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p_1/1.e9,vs_prem/1.e3,color='k',linestyle='-',marker='x',\
    markerfacecolor='k',markersize=4)
    plt.title("Vs (km/s)")
    
    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,mat_vphi_pyro/1.e3,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4)
    plt.plot(seis_p_2/1.e9,mat_vphi_enst/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p_1/1.e9,v_phi_prem/1.e3,color='k',linestyle='-',marker='x',\
    markerfacecolor='k',markersize=4)
    plt.title("Vphi (km/s)")
    
    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p_1/1.e9,mat_rho_pyro/1.e3,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4,label="C-chondrite")
    plt.plot(seis_p_2/1.e9,mat_rho_enst/1.e3,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label="enstatite")
    plt.plot(seis_p_1/1.e9,rho_prem/1.e3,color='k',linestyle='-',marker='x',\
    markerfacecolor='k',markersize=4,label="PREM")
    plt.title("density (kg/m^3)")
    plt.legend(loc='lower right')
    plt.xlabel("Pressure (GPa)")
    
    
    #plot percent differences
    mat_vs_enst_interp = np.interp(seis_p_1, seis_p_2,mat_vs_enst)
    mat_vphi_enst_interp = np.interp(seis_p_1, seis_p_2,mat_vphi_enst)
    mat_rho_enst_interp = np.interp(seis_p_1, seis_p_2,mat_rho_enst)
    
    per_diff_vs = 100*(mat_vs_pyro - mat_vs_enst_interp)/mat_vs_pyro
    per_diff_vphi = 100*(mat_vphi_pyro - mat_vphi_enst_interp)/mat_rho_pyro
    per_diff_rho = 100*(mat_rho_pyro - mat_rho_enst_interp)/mat_rho_pyro
    
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,per_diff_vs,color='g',linestyle='-',marker='o',\
    markerfacecolor='k',markersize=4,label="vs")
    plt.plot(seis_p_1/1.e9,per_diff_vphi,color='c',linestyle='-',marker='+',\
    markerfacecolor='k',markersize=4,label="vphi")
    plt.plot(seis_p_1/1.e9,per_diff_rho,color='m',linestyle='-',marker='x',\
    markerfacecolor='k',markersize=4,label="density")
    
    plt.title("percent difference")
    plt.legend(loc='center right')
    plt.xlabel("Pressure (GPa)")
    plt.savefig("output_figures/example_compare_enstpyro.png") 
    plt.show()
