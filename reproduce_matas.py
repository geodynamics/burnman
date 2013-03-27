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
    
    method = 'mgd' 
    
    
    # Input for model M_pyrolite as defined in Matas et al. 2007 table 3. Molar proportions are converted to atomic fractions    

    #weight_percents = {'Mg':0.245, 'Fe': 0.033, 'Si':0.223, 'Ca':0., 'Al':0.}
    #phase_fractions,relative_molar_percent = burnman.calculate_phase_percents(weight_percents)
    rock = burnman.composite( ( (minerals.Matas_mg_perovskite(),.642 ),
                            (minerals.Matas_fe_perovskite(), .056 ),
                            (minerals.Matas_periclase(), .246 ),
                            (minerals.Matas_wuestite(), .056 )))
    #                        (minerals.Murakami_fe_periclase(), phase_fractions['fp']) ) )                                  
    #input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
    rock.set_method(method) 
    seis_p_1 = np.arange(28e9, 128e9, 4.8e9)
    geotherm = burnman.geotherm.brown_shankland
    temperature_bs = [geotherm(p) for p in seis_p_1]
 
 
    #Now we'll calculate the models. 
    
    mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.calculate_velocities(seis_p_1, temperature_bs, rock)    
    


    # seismic velocities for comparison
    class ak135_table(burnman.seismic.radiustable):
        def __init__(self):
            burnman.seismic.radiustable.__init__(self) 
            table = burnman.tools.read_table("data/ak135_lowermantle.txt") # radius, pressure, density, v_p, v_s
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
 
    ##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)
    
    plt.subplot(2,2,2)
    plt.plot(seis_p_1/1.e9,mat_vs_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_1/1.e9,seis_vs/1.e3,color='k',linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.title("Vs (km/s)")
    plt.ylim(5,7.6) 
    plt.xlim(29,131)
 
    # plot Vp
    plt.subplot(2,2,3)
    plt.title("Vp (km/s)")
    plt.plot(seis_p_1/1.e9,mat_vp_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_1/1.e9,seis_vp/1.e3,color='k', linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.ylim(9.25,14.0)    
    plt.xlim(29,131)

    # plot density
    plt.subplot(2,2,4)
    plt.plot(seis_p_1/1.e9,mat_rho_1/1.e3,color='b',linestyle='-')
    plt.plot(seis_p_1/1.e9,seis_rho/1.e3,color='k', linestyle='',marker='o',markerfacecolor='w',markersize=4)
    plt.title("Density (kg/m^3)")
    plt.xlim(29,131)    
    
    
    plt.subplot(2,2,1)
    plt.plot(seis_p_1/1.e9,temperature_bs,color='k',linestyle='-',label='brown_shankland')
    plt.ylim(1600,3100)
    plt.xlim(29,131)
    plt.title("temperature")
    plt.legend(loc='upper left')
    
    plt.savefig("matas.png") 
    #plt.show()


    plt.subplot(1,1,1)
    plt.ylim(-2,5) 
    plt.xlim(800,3000)
    fig1 = mpimg.imread('data/matas_vs_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vs_1-seis_vs)/seis_vs,color='b',linestyle='--')
    plt.savefig("matas_cmp_vs.png")
    plt.show()
    plt.close()

    plt.subplot(1,1,1)
    plt.ylim(-2,5)
    plt.xlim(800,3000)
    fig1 = mpimg.imread('data/matas_vp_forcomparison.png')
    plt.imshow(fig1, extent=[800,3000,-2,5], aspect='auto')
    d=np.array(depths)/1.e3
    plt.plot(d,100.*(mat_vp_1-seis_vp)/seis_vp,color='b',linestyle='--')
    plt.savefig("matas_cmp_vp.png")
    plt.show()

