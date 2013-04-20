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
from burnman.materials import composite
import matplotlib.image as mpimg

if __name__ == "__main__":    
    
    ###Input Model 1
    
    #INPUT for method_1
    """ choose 'slb' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)"""
    
    method = 'slb'
    
    
    #Input composition of model 1. See example_composition for potential choices. We'll just choose something simple here
        
    amount_perovskite_1 = 1.0
    rock_1 = composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite_1),
                            (minerals.Murakami_fe_periclase(), 1.0-amount_perovskite_1) ) )
    rock_1.set_method(method)
 
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    
    seis_p_1 = np.arange(28.5e9, 133e9, 4.75e9)
    #seismic model for comparison:
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    depths = map(seismic_model.depth, seis_p_1)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)    
    geotherm = burnman.geotherm.brown_shankland
    temperature_bs = [geotherm(p) for p in seis_p_1]
   ##Now onto the second model parameters
    
    
    
    #Input composition of model 2. See example_composition for potential choices. We'll just choose something simple here
        
    amount_perovskite_2 = 0.95
    rock_2 = composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite_2),
                            (minerals.Murakami_fe_periclase(), 1.0-amount_perovskite_2) ) )
    rock_2.set_method(method) 
	#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    
    seis_p_2 = seis_p_1
    
    #input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.
    
    
    ##Now onto the third model parameters



    #Input composition of model 3. See example_composition for potential choices. We'll just choose something simple here

    amount_perovskite_3 = 0.0
    rock_3 = composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite_3),
                            (minerals.Murakami_fe_periclase(), 1.0-amount_perovskite_3) ) )
    rock_3.set_method(method) 
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

    seis_p_3 = seis_p_1 

    #input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

   


    #Input composition of model 4. See example_composition for potential choices. We'll just choose something simple here

    amount_perovskite_4 = 0.8
    rock_4 = composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite_4),
                            (minerals.Murakami_fe_periclase(), 1.0-amount_perovskite_4) ) )
    rock_4.set_method(method) 

    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

    seis_p_4 = seis_p_1 

    #input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

 
    #Now we'll calculate the models. 
    

    T0 = 1600.
    temperature_1 = burnman.geotherm.self_consistent(seis_p_1, T0, rock_1)
    
    
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.calculate_velocities(seis_p_1, temperature_bs, rock_1)    
    

    temperature_2 = burnman.geotherm.self_consistent(seis_p_1, T0, rock_2)
    
    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_mu_2 = burnman.calculate_velocities(seis_p_2, temperature_bs, rock_2)    
    #mat_rho_t, mat_vp_t, mat_vs_t, mat_vphi_t, mat_K_t, mat_mu_t = burnman.calculate_velocities(seis_p_2, temperature_bs, rock_2.phases[0]) 

    temperature_3 = burnman.geotherm.self_consistent(seis_p_1, T0, rock_3)
   
    mat_rho_3, mat_vp_3, mat_vs_3, mat_vphi_3, mat_K_3, mat_mu_3 = burnman.calculate_velocities(seis_p_3, temperature_bs, rock_3)


    temperature_4 = burnman.geotherm.self_consistent(seis_p_1, T0, rock_4)
  
    mat_rho_4, mat_vp_4, mat_vs_4, mat_vphi_4, mat_K_4, mat_mu_4 = burnman.calculate_velocities(seis_p_4, temperature_bs, rock_4)
    
    chi_den_murk = 0.
    for i in range(len(mat_rho_2)):
	chi_den_murk = chi_den_murk + pow(mat_rho_2[i] - seis_rho[i],2.)/seis_rho[i]
    	
    chi_vs_murk = 0.
    for i in range(len(mat_rho_2)):
    	chi_vs_murk = chi_vs_murk + pow(mat_vs_2[i] - seis_vs[i],2.)/seis_vs[i]
    
    chi_vp_murk = 0.
    for i in range(len(mat_rho_2)):
    	chi_vp_murk = chi_vp_murk + pow(mat_vp_2[i] - seis_vp[i],2.)/seis_vp[i]
    
    chi_den_py = 0.
    for i in range(len(mat_rho_4)):
    	chi_den_py = chi_den_py + pow(mat_rho_4[i] - seis_rho[i],2.)/seis_rho[i]
    	
    chi_vs_py = 0.
    for i in range(len(mat_rho_4)):
    	chi_vs_py = chi_vs_py + pow(mat_vs_4[i] - seis_vs[i],2.)/seis_vs[i]
    
    chi_vp_py = 0.
    for i in range(len(mat_rho_4)):
    	chi_vp_py = chi_vp_py + pow(mat_vp_4[i] - seis_vp[i],2.)/seis_vp[i]
    	
    
 	
    ##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)
    print "Chi squared Vs: ", chi_vs_murk
    print "Chi squared Vs pyrolite: ", chi_vs_py
    print "	"
    print "Chi squared Vp: ", chi_vp_murk
    print "Chi squared Vp pyrolite: ", chi_vp_py  
    print "	"  
    print "Chi squared Rho: ", chi_den_murk 
    print "Chi squared Rho pyrolite: ", chi_den_py
    den = chi_den_murk/chi_den_py
    vs = chi_vs_murk/chi_vs_py
    vp = chi_vp_murk/chi_vp_py
    #plt.subplot(2,2,2)
    #plt.plot(seis_p_1/1.e9,mat_vs_1,color='b',linestyle='-',label='pv + 4% Al')
    #plt.plot(seis_p_2/1.e9,mat_vs_2,color='r',linestyle='-',label='Xpv = 0.95')
    #plt.plot(seis_p_3/1.e9,mat_vs_3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4,label='ferropericlase Xpv = 0.8')
    #plt.plot(seis_p_4/1.e9,mat_vs_4,color='g',linestyle='-',label='pyrolite')
    #plt.plot(seis_p_1/1.e9,seis_vs,color='k',linestyle='',marker='o',markerfacecolor='w',markersize=4,label='PREM')
    #plt.title("Vs (km/s)")
    #plt.legend(loc='upper left')
    #plt.ylim(5,7.6) 
    #plt.xlim(29,131)
 
    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_2/1.e9,mat_vs_2/1.e3,color='r',linestyle='-')
    plt.plot(seis_p_3/1.e9,mat_vs_3/1.e3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
    plt.plot(seis_p_4/1.e9,mat_vs_4/1.e3,color='g',linestyle='-')
    plt.plot(seis_p_1/1.e9,seis_vs,color='k',linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.title("Vs (km/s)")
    plt.text(40,7.3,"chi^2 murk/py= %3.3f" % vs)
    plt.ylim(5,7.6) 
    plt.xlim(29,131)
    
    # plot Vp
    plt.subplot(2,2,3)
    plt.title("Vp (km/s)")
    plt.plot(seis_p_1/1.e9,mat_vp_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_2/1.e9,mat_vp_2/1.e3,color='r',linestyle='-')
    plt.plot(seis_p_3/1.e9,mat_vp_3/1.e3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
    plt.plot(seis_p_4/1.e9,mat_vp_4/1.e3,color='g',linestyle='-') 
    plt.plot(seis_p_1/1.e9,seis_vp,color='k', linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.text(50,9.7,"chi^2 murk/py= %3.3f" % vp)
    plt.ylim(9.25,14.0)    
    plt.xlim(29,131)

    # plot density
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,mat_rho_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_2/1.e9,mat_rho_2/1.e3,color='r',linestyle='-')
    plt.plot(seis_p_3/1.e9,mat_rho_3/1.e3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
    plt.plot(seis_p_4/1.e9,mat_rho_4/1.e3,color='g',linestyle='-')
    plt.plot(seis_p_1/1.e9,seis_rho/1.e3,color='k', linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.text(40,4.3,"chi^2 murk/py= %3.3f" % den)
    plt.title("Density (kg/m^3)")
    plt.xlim(29,131)    
    
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,temperature_bs,color='k',linestyle='-',label='brown_shankland')
    plt.plot(seis_p_1/1.e9,temperature_4,color='k',linestyle='-',label='self-consistent pyrolite')
    plt.ylim(1600,3100)
    plt.xlim(29,131)
    plt.title("temperature")
    plt.legend(loc='upper left')
    
    #plt.savefig("murakami_fig4.png") 
    plt.show()


    plt.subplot(1,1,1)
    plt.ylim(5,7.6) 
    plt.xlim(25,135)
    fig1 = mpimg.imread('input_figures/murakami_vs_cmp.png')
    plt.imshow(fig1, extent=[25,135,5.0,7.6], aspect='auto')
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='b',linestyle='--',marker='o',markerfacecolor='b',markersize=3)
    plt.plot(seis_p_2/1.e9,mat_vs_2/1.e3,color='r',linestyle='--',marker='o',markerfacecolor='r',markersize=3)
    #plt.plot(seis_p_1/1.e9,(0.8*mat_vs_1+0.2*mat_vs_3)/1.e3,color='m')
    plt.plot(seis_p_3/1.e9,mat_vs_3/1.e3,color='k',linestyle='--',marker='o',markerfacecolor='k',markersize=3)
    plt.plot(seis_p_4/1.e9,mat_vs_4/1.e3,color='g',linestyle='--',marker='o',markerfacecolor='g',markersize=3)
    plt.plot(seis_p_1/1.e9,seis_vs/1.e3,color='k',linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.title("Vs (km/s)")
    plt.savefig("murakami_fig4.png")
    plt.show()
