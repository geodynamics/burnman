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
import matplotlib.image as mpimg
import numpy

if __name__ == "__main__":    
    
    ###Input Model 1
    
    #INPUT for method_1
    """ choose 'slb' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)"""
    
    method = 'mgd2' 
    
    
    # Input for model M_pyrolite as defined in Matas et al. 2007 table 3. Molar proportions are converted to atomic fractions    

    #weight_percents = {'Mg':0.297882, 'Fe': 0.0489699, 'Si':0.1819, 'Ca':0.0228576, 'Al':0.0116446}
    #phase_fractions,relative_molar_percent = burnman.calculate_phase_percents(weight_percents)
    rock = burnman.composite( ( (minerals.Matas_etal_2007.mg_perovskite(),.620 ),
                            (minerals.Matas_etal_2007.fe_perovskite(), .078 ),
                            (minerals.Matas_etal_2007.periclase(), .268 ),
                            (minerals.Matas_etal_2007.wuestite(), .034 )))
    rock2 = burnman.composite( ( (minerals.Matas_etal_2007.mg_perovskite(),.62 ),
                            (minerals.Matas_etal_2007.fe_perovskite(), .078/1.46 ),
                            (minerals.Matas_etal_2007.periclase(), .268 ),
                            (minerals.Matas_etal_2007.wuestite(), .034/1.46 ))) 
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
    
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_G_1 = burnman.velocities_from_rock(rock,seis_p_1, temperature_bs)    
    mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_G_2 = burnman.velocities_from_rock(rock2,seis_p_1, temperature_bs) 

    #Next, we calculate the velocites with 3rd order Birch-Murnaghan
    method='mgd2'
    rock.set_method(method)
    rock2.set_method(method)
    mat_rho_1_3, mat_vp_1_3, mat_vs_1_3, mat_vphi_1_3, mat_K_1_3, mat_G_1_3 = burnman.velocities_from_rock(rock,seis_p_1, temperature_bs)
    mat_rho_2_3, mat_vp_2_3, mat_vs_2_3, mat_vphi_2_3, mat_K_2_3, mat_G_2_3 = burnman.velocities_from_rock(rock2,seis_p_1, temperature_bs)

    
    # seismic velocities for comparison
    class ak135_table(burnman.seismic.radiustable):
        def __init__(self):
            burnman.seismic.radiustable.__init__(self) 
            table = burnman.tools.read_table("input_seismic/ak135_lowermantle.txt") # radius, pressure, density, v_p, v_s
            table = np.array(table)
            self.table_radius = table[:,0]
            self.table_pressure = table[:,1]
            self.table_density = table[:,2]
            self.table_vp = table[:,3]
            self.table_vs = table[:,4]


    ak=ak135_table()
    #seismic model for comparison:
    depths = map(ak.depth, seis_p_1)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = ak.evaluate_all_at(depths)
    prem=burnman.seismic.prem()
    prem_p, prem_rho, prem_vp, prem_vs, prem_vphi = prem.evaluate_all_at(depths) 

    ##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)

    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,seis_vs/1.e3,color='g',linestyle='-',label='ak135')
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='b',linestyle='-',label='rock')
    plt.plot(seis_p_1/1.e9,mat_vs_2/1.e3,color='r',linestyle='-',label='rock2')
    #plt.plot(depths,prem_vs/1.e3,color='r')
    plt.title("Vs (km/s)")
    plt.ylim(5,7.6) 
#    plt.xlim(29,131)
    plt.legend(loc='upper left')
 
    # plot Vp
    plt.subplot(2,2,3)
    plt.title("Vp (km/s)")
    plt.plot(seis_p_1/1.e9,seis_vp/1.e3,color='g',linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_vp_1/1.e3,color='b', linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_vp_2/1.e3,color='r', linestyle='-')
    plt.ylim(6.25,14.0)    
    plt.xlim(29,131)

    # plot density
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,seis_rho/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_rho_1/1.e3,color='g', linestyle='-')
    plt.plot(seis_p_1/1.e9,mat_rho_2/1.e3,color='r', linestyle='-')
    plt.title("Density (kg/m^3)")
    plt.xlim(29,131)    
    
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,temperature_bs,color='k',linestyle='-',label='brown_shankland')
    plt.ylim(1600,3100)
    plt.xlim(29,131)
    plt.title("temperature")
    plt.legend(loc='upper left')
    
    plt.savefig("output_figures/reproduce_matas.png") 
    #plt.show()
    #plt.close()

    plt.subplot(1,1,1)
    plt.ylim(-2,5) 
    plt.xlim(800,3000)
    fig1 = mpimg.imread('input_figures/matas_vs_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vs_1-seis_vs)/seis_vs,color='b',linestyle='--')
    plt.plot(d,100.*(mat_vs_2-seis_vs)/seis_vs,color='r',linestyle='--')
    plt.title("Comparison M-pyrolite model Vs, 2nd order")
    plt.savefig("output_figures/reproduce_matas_cmp_vs.png")
    plt.show()
    plt.close()

    plt.subplot(1,1,1)
    plt.ylim(-2,5)
    plt.xlim(800,3000)
    fig1 = mpimg.imread('input_figures/matas_vp_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vp_1-seis_vp)/seis_vp,color='b',linestyle='--')
    plt.plot(d,100.*(mat_vp_2-seis_vp)/seis_vp,color='r',linestyle='--')
    plt.title("Comparison M-pyrolite model Vp, 2nd order")
    plt.savefig("output_figures/reproduce_matas_cmp_vp.png")
    plt.show()

    plt.subplot(1,1,1)
    plt.ylim(-2,5)
    plt.xlim(800,3000)
    fig1 = mpimg.imread('input_figures/matas_vs_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vs_1_3-seis_vs)/seis_vs,color='b',linestyle='--')
    plt.plot(d,100.*(mat_vs_2_3-seis_vs)/seis_vs,color='r',linestyle='--')
    plt.title("Comparison M-pyrolite model Vs, 3rd order")
    plt.savefig("output_figures/reproduce_matas_cmp_vs.png")
    plt.show()
    plt.close()

    plt.subplot(1,1,1)
    plt.ylim(-2,5)
    plt.xlim(800,3000)
    fig1 = mpimg.imread('input_figures/matas_vp_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vp_1_3-seis_vp)/seis_vp,color='b',linestyle='--')
    plt.plot(d,100.*(mat_vp_2_3-seis_vp)/seis_vp,color='r',linestyle='--')
    plt.title("Comparison M-pyrolite model Vp, 3rd order")
    plt.savefig("output_figures/reproduce_matas_cmp_vp.png")
    plt.show()
