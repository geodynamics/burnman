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
import numpy

if __name__ == "__main__":    
    
    ###Input Model 1
    
    #INPUT for method_1
    """ choose 'slb' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)"""
    
    method = 'bm3' 
    
    
    # Input for model M_pyrolite as defined in Matas et al. 2007 table 3. Molar proportions are converted to atomic fractions    

    #weight_percents = {'Mg':0.297882, 'Fe': 0.0489699, 'Si':0.1819, 'Ca':0.0228576, 'Al':0.0116446}
    #phase_fractions,relative_molar_percent = burnman.calculate_phase_percents(weight_percents)
    rock = burnman.composite( ( (minerals.SLB2011_mg_perovskite(),0.5*1.00),
                            (minerals.SLB2011_fe_perovskite(), .9*0.0 ),
                            (minerals.SLB2011_periclase(), 0.5*1.00 ),
                            (minerals.SLB2011_wuestite(), 0.1*0.0 )))
    rock2 = burnman.composite( ( (minerals.mg_perovskite(),.62 ),
                            (minerals.fe_perovskite(), .078/1.5 ),
                            (minerals.periclase(), .28 ),
                            (minerals.wuestite(), .034/1.5 ))) 
    #KD is 2... doesn't match either
    #rock2 = burnman.composite( ( (minerals.Matas_mg_perovskite(),.3574 ),
    #                        (minerals.Matas_fe_perovskite(), .0536 ),
    #                        (minerals.Matas_periclase(), .1656 ),
    #                        (minerals.Matas_wuestite(), .0124 )))
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    rock.set_method(method) 
    rock2.set_method(method)
    seis_p_1 = np.arange(28e9, 128e9, 4.8e9)
    temperature_bs = burnman.geotherm.brown_shankland(seis_p_1)
 
 
    #Now we'll calculate the models. 
    T0 = 0.


    temperature_1 = burnman.geotherm.adiabatic(seis_p_1, T0, rock)
 
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.velocities_from_rock(rock,seis_p_1, temperature_1)  
    temperature_2 = burnman.geotherm.adiabatic(seis_p_1, T0, rock2)  
    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_mu_2 = burnman.velocities_from_rock(rock2,seis_p_1, temperature_2) 
    


    prem=burnman.seismic.prem()
    depths=map(prem.depth,seis_p_1)
    prem_p, prem_rho, prem_vp, prem_vs, prem_vphi = prem.evaluate_all_at(depths) 

    ##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)

    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,prem_vs/1.e3,color='g',linestyle='-',label='ak135')
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='b',linestyle='-',label='rock')
   # plt.plot(seis_p_1/1.e9,mat_vs_2/1.e3,color='r',linestyle='-',label='rock2')
    #plt.plot(depths,prem_vs/1.e3,color='r')
    plt.title("Vs (km/s)")
    plt.ylim(5,7.6) 
#    plt.xlim(29,131)
    plt.legend(loc='upper left')
 
    # plot Vp
    plt.subplot(2,2,3)
    plt.title("Vp (km/s)")
    plt.plot(seis_p_1/1.e9,prem_vp/1.e3,color='g',linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_vp_1/1.e3,color='b', linestyle='-')
   # plt.plot(seis_p_1/1.e9,mat_vp_2/1.e3,color='r', linestyle='-')
    plt.ylim(6.25,14.0)    
    plt.xlim(29,131)

    # plot density
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,prem_rho/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_rho_1/1.e3,color='g', linestyle='-')
    #plt.plot(seis_p_1/1.e9,mat_rho_2/1.e3,color='r', linestyle='-')
    plt.title("Density (kg/m^3)")
    plt.xlim(29,131)    
    
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,temperature_1,color='k',linestyle='-',label='adiabat')
    plt.ylim(1600,3100)
    plt.xlim(29,131)
    plt.title("temperature")
    plt.legend(loc='upper left')
    
    plt.savefig("output_figures/reproduce_matas.png") 
    plt.show()

