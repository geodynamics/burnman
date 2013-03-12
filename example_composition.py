# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows how to create different minerals, how to compute seismic
velocities, and how to compare them to a seismic reference model.

requires:
- geotherms
- seismic models
- compute seismic velocities

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
    """ choose 'slb' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005) 
    or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)"""
    
    method = 'slb' 
        
    
    # To compute seismic velocities and other properties, we need to supply
    # burnman with a list of minerals (phases) and their molar abundances. Minerals
    # are classes found in burnman.minerals and are derived from
    # burnman.minerals.material.
    # Here are a few ways to define phases and molar_abundances:
    
    #Example 1: two simple fixed minerals
    if True:
        amount_perovskite = 0.95
        rock = burnman.composite ( ( (minerals.Murakami_fe_perovskite(), amount_perovskite),
                                     (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite)) )
    
    #Example 2: specify fixed iron content
    if False:
        amount_perovskite = 0.95
        rock = burnman.composite( ( (minerals.mg_fe_perovskite(0.8), amount_perovskite), 
                                    (minerals.ferropericlase(0.8), 1.0-amount_perovskite) ) )
    
    #Example 3: input weight percentages
    #See comments in burnman/composition.py for references to partition coefficent calculation

    if False:
        weight_percents = {'Mg':0.213, 'Fe': 0.08, 'Si':0.27, 'Ca':0., 'Al':0.}
        phase_fractions,relative_molar_percent = burnman.calculate_phase_percents(weight_percents)
        iron_content = lambda p,t: burnman.calculate_partition_coefficient(p,t,relative_molar_percent)

        rock = burnman.composite ( ( (minerals.mg_fe_perovskite_pt_dependent(iron_content,0), phase_fractions['pv'] ),
                                     (minerals.ferropericlase_pt_dependent(iron_content,1), phase_fractions['fp'] ) ) )
        
    #Example 4: three materials
    if False:
        rock = burnman.composite ( ( (minerals.Murakami_fe_perovskite(), 0.7),
                                     (minerals.ferropericlase(0.5), 0.2) ,
                                     (minerals.stishovite(), 0.1 ) ) )
    
    
    #seismic model for comparison:
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700, 2800, number_of_points)
    #alternatively, we could use the values where prem is defined:
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in seis_p]
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(seis_p, temperature, rock)    
    
    [rho_err,vphi_err,vs_err]=burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
        
    
    # PLOTTING
    
    # plot vs
    plt.subplot(2,2,1)
    plt.plot(seis_p/1.e9,mat_vs,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4,label='computation')
    plt.plot(seis_p/1.e9,seis_vs,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4,label='reference')
    plt.title("Vs (km/s)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.ylim(5.1,7.6)
    plt.legend(loc='lower right')
    plt.text(40,7.3,"misfit= %3.3f" % vs_err)
    
    # plot Vphi
    plt.subplot(2,2,2)
    plt.plot(seis_p/1.e9,mat_vphi,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,seis_vphi,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
    plt.title("Vphi (km/s)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.ylim(7,12)
    plt.text(40,11.5,"misfit= %3.3f" % vphi_err)
    
    # plot density
    plt.subplot(2,2,3)
    plt.plot(seis_p/1.e9,mat_rho,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
    plt.plot(seis_p/1.e9,seis_rho,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
    plt.title("density (kg/m^3)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.text(40,4.3,"misfit= %3.3f" % rho_err)
    plt.xlabel("Pressure (GPa)")
    
    
    # plot geotherm
    plt.subplot(2,2,4)
    plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
    plt.title("Geotherm (K)")
    plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
    plt.xlabel("Pressure (GPa)")
    
    plt.savefig("example_composition.png")
    plt.show()
