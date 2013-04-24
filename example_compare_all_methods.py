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
    #Input composition. 

    amount_perovskite = 0.95
    rock = burnman.composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite),
                            (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite) ) )

    seis_p = np.arange(25e9, 125e9, 5e9)

    T0 = 1500.0
    
    #Now we'll calculate the models. 
    
    rock.set_method('mgd3')
    temperature = burnman.geotherm.self_consistent(seis_p, T0, rock)    
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())

    rock.set_method('slb2')
    temperature = burnman.geotherm.self_consistent(seis_p, T0, rock)    
    
    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_mu_2 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())

    rock.set_method('slb3')
    temperature = burnman.geotherm.self_consistent(seis_p, T0, rock)    

    mat_rho_3, mat_vp_3, mat_vs_3, mat_vphi_3, mat_K_3, mat_mu_3 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())

    rock.set_method('bm2')
    temperature = burnman.geotherm.self_consistent(seis_p, T0, rock)    

    mat_rho_4, mat_vp_4, mat_vs_4, mat_vphi_4, mat_K_4, mat_mu_4 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())

    rock.set_method('bm3')
    temperature = burnman.geotherm.self_consistent(seis_p, T0, rock)   
    mat_rho_5, mat_vp_5, mat_vs_5, mat_vphi_5, mat_K_5, mat_mu_5 = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())
        
    ##Now let's plot the comparison. You can conversely just output to a data file 
    #(see example_woutput.py)
    
    plt.subplot(2,2,1)
    plt.plot(seis_p/1.e9,mat_vs_1/1.e3,color='r',linestyle='-',marker='^', \
    markerfacecolor='r',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_2/1.e3,color='k',linestyle='-',marker='v',\
    markerfacecolor='k',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_3/1.e3,color='b',linestyle='-',marker='x',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_4/1.e3,color='g',linestyle='-',marker='o',\
    markerfacecolor='g',markersize=4)
    plt.plot(seis_p/1.e9,mat_vs_5/1.e3,color='y',linestyle='-',marker='*',\
    markerfacecolor='y',markersize=4)
    plt.title("Vs (km/s)")
    
    
    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p/1.e9,mat_vphi_1/1.e3,color='r',linestyle='-',marker='^',\
    markerfacecolor='r',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_2/1.e3,color='k',linestyle='-',marker='v',\
    markerfacecolor='k',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_3/1.e3,color='b',linestyle='-',marker='x',\
    markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_4/1.e3,color='g',linestyle='-',marker='o',\
    markerfacecolor='g',markersize=4)
    plt.plot(seis_p/1.e9,mat_vphi_5/1.e3,color='y',linestyle='-',marker='*',\
    markerfacecolor='y',markersize=4)

    plt.title("Vphi (km/s)")
    
    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p/1.e9,mat_rho_1/1.e3,color='r',linestyle='-',marker='^',\
    markerfacecolor='r',markersize=4,label='mgd')
    plt.plot(seis_p/1.e9,mat_rho_2/1.e3,color='k',linestyle='-',marker='v',\
    markerfacecolor='k',markersize=4,label='slb')
    plt.plot(seis_p/1.e9,mat_rho_3/1.e3,color='b',linestyle='-',marker='x',\
    markerfacecolor='b',markersize=4,label='slb3')
    plt.plot(seis_p/1.e9,mat_rho_4/1.e3,color='g',linestyle='-',marker='o',\
    markerfacecolor='g',markersize=4,label='bm2')
    plt.plot(seis_p/1.e9,mat_rho_5/1.e3,color='y',linestyle='-',marker='*',\
    markerfacecolor='y',markersize=4,label='bm3')
    plt.title("density (kg/m^3)")
    plt.legend(loc='upper left')
    
    
    plt.savefig("output_figures/example_compare_all_methods.png")    
    plt.show()
