# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Calculates and plots two models for different minerals or methods and plots
the results. Calculates basic percent difference statistics as well.

requires:
- geotherms
- creating minerals

teaches:
- compute seismic velocities and compare
"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
    ###Input Model 1
    
    #INPUT for method_1
    """ choose 'slb2' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus, matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shearl modulus, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))"""
    
    method_1 = 'mgd3' 
    
    
    #Input composition of model 1. See example_composition for potential choices. We'll just choose something simple here
        
    amount_perovskite_1 = 0.95
    rock_1 = burnman.composite ( ((minerals.Murakami_fe_perovskite(), amount_perovskite_1),
                                  (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite_1)) ) 
    rock_1.set_method(method_1)
    
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    
    seis_p_1 = np.arange(25e9, 125e9, 5e9)
    
    #input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.
    
    geotherm = burnman.geotherm.brown_shankland
    temperature_1 = [geotherm(p) for p in seis_p_1]
    
    ##Now onto the second model parameters
    
    ##input second method
    method_2 = 'slb3' 
    
    
    #Input composition of model 2. See example_composition for potential choices. We'll just choose something simple here
        
    amount_perovskite_2 = 0.95
    rock_2 = burnman.composite ( ((minerals.Murakami_fe_perovskite(), amount_perovskite_2),
                                  (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite_2)) ) 
    rock_2.set_method(method_2)
    
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    
    seis_p_2 = np.arange(25e9, 125e9, 5e9)
    
    #input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.
    
    geotherm = burnman.geotherm.brown_shankland
    temperature_2 = [geotherm(p) for p in seis_p_2]
    
    
    #Now we'll calculate the models. 
    
    
    print "Calculations are done for:"
    for ph in rock_1.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.calculate_velocities(seis_p_1, temperature_1, rock_1)    
    
    print "Calculations are done for:"
    for ph in rock_2.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_mu_2 = burnman.calculate_velocities(seis_p_2, temperature_2, rock_2)    
    
    
    ##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='mgd')
    plt.plot(seis_p_2/1.e9,mat_vs_2/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='slb')
    plt.title("Vs (km/s)")
    
    
    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,mat_vphi_1/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='mgd')
    plt.plot(seis_p_2/1.e9,mat_vphi_2/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='slb')
    plt.title("Vphi (km/s)")
    
    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p_1/1.e9,mat_rho_1/1.e3,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='mgd')
    plt.plot(seis_p_2/1.e9,mat_rho_2/1.e3,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='slb')
    plt.title("density (kg/m^3)")
    plt.legend(loc='upper left')
    
    
    #plot percent differences
    mat_vs_2_interp = np.interp(seis_p_1, seis_p_2,mat_vs_2)
    mat_vphi_2_interp = np.interp(seis_p_1, seis_p_2,mat_vphi_2)
    mat_rho_2_interp = np.interp(seis_p_1, seis_p_2,mat_rho_2)
    
    per_diff_vs = 100*(mat_vs_1 - mat_vs_2_interp)/mat_vs_1
    per_diff_vphi = 100*(mat_vphi_1 - mat_vphi_2_interp)/mat_rho_1
    per_diff_rho = 100*(mat_rho_1 - mat_rho_2_interp)/mat_rho_1
    
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,per_diff_vs,color='g',linestyle='-',marker='o',markerfacecolor='k',markersize=4,label='V_s')
    plt.plot(seis_p_1/1.e9,per_diff_vphi,color='c',linestyle='-',marker='+',markerfacecolor='k',markersize=4,label='V_phi')
    plt.plot(seis_p_1/1.e9,per_diff_rho,color='m',linestyle='-',marker='x',markerfacecolor='k',markersize=4,label='density')
    
    plt.title("percent difference")
    plt.legend(loc='center left')
    
    plt.show()
