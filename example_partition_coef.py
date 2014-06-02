# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows how to vary the distribution coefficient of the
perovskite/ferropericlase system. The user sets Kd_0 and BurnMan scales Kd as
a function of P and T adopting the formalism of Nakajima et al.,
2012. Specifically we adopt equation 5 of Nakajima et al., 2012 with DeltaV_0
= 0.2 cc/mol, and calculating the partition coefficient of Fe in each phase
from stoichiometry.

This example will calculate mineral input parameters from Mg and Fe endmembers
from Stixrude and Lithgow-bertelloni, 2005 with weighting determined by the
calculated partition coefficients. Finally, output plots of X_Fe pv and X_Fe
fp our output as well as the user's choice of geotherm

requires:
- geotherms
-input distribution coefficient Kd_0

teaches:
- creating varying proportions of Fe and its effect on seismic velocities


"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
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
        

    
    #input weight percentages and distribution coefficient 
    #See comments in burnman/composition.py for references to 
    #partition coefficent calculation

    weight_percents = {'Mg':0.213, 'Fe': 0.08, 'Si':0.27, 'Ca':0., 'Al':0.}
    phase_fractions,relative_molar_percent = \
    burnman.calculate_phase_percents(weight_percents)
    
    Kd_0 = .5 #Fig 5 Nakajima et al 2012
    
    iron_content = lambda p,t: \
    burnman.calculate_partition_coefficient(p,t,relative_molar_percent,Kd_0)

    rock = burnman.composite ( [ (minerals.SLB_2005.mg_fe_perovskite_pt_dependent(iron_content,1),\
    								phase_fractions['pv'] ), 
    								(minerals.SLB_2005.ferropericlase_pt_dependent(iron_content,0),\
    								phase_fractions['fp'] ) ] )
            
    #seismic model for comparison:
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() 
    #(see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    temperature = burnman.geotherm.brown_shankland(seis_p)

    part_coef_pv=[0 for x in seis_p]
    part_coef_fp=[0 for x in seis_p]
    
    for t in range(0,number_of_points):
        part_coef_fp[t],part_coef_pv[t] = burnman.calculate_partition_coefficient \
        (seis_p[t],temperature[t],relative_molar_percent,Kd_0)
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    rock.debug_print()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())
    
    [rho_err,vphi_err,vs_err]=burnman.compare_chifactor\
    ([mat_vs,mat_vphi,mat_rho],[seis_vs,seis_vphi,seis_rho])
     

    # PLOTTING
    
    #plot Fe_content_pv vs P
    plt.subplot(2,2,1)
    plt.plot(seis_p/1e9,part_coef_pv,color='b',linestyle='-',marker='o',\
    markerfacecolor='b',markersize=4,label='Perovskite')
    plt.title("Fe content pv")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.legend(loc='lower right')
    
    
    
    # plot Fe_content_fp
    plt.subplot(2,2,2)
    plt.plot(seis_p/1e9,part_coef_fp,color='b',linestyle='-',marker='o',\
    markerfacecolor='g',markersize=4,label='Ferropericlase')
    plt.title("Fe content fp")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.legend(loc='lower right')

    
    
    # plot Geotherm
    plt.subplot(2,2,3)
    plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',\
    markerfacecolor='r',markersize=4)
    plt.title("Geotherm (K)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.xlabel("Pressure (GPa)")
    plt.text(40,2400,"Kd_0 = %3.3f" % Kd_0)

    plt.show()
    
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
    
    plt.show()
