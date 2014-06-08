# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- compare errors between models
- loops over models

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    

    
    ### input variables ###
    #######################
    
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
    method = 'slb2' 
    
    seismic_model = burnman.seismic.PREM() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    number_of_points = 20 #set on how many depth slices the computations should be done
    depths = np.linspace(700e3,2800e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
    temperature = burnman.geotherm.brown_shankland(seis_p)
    
    def material_error(amount_perovskite):
        rock = burnman.Composite ( [amount_perovskite, 1.0-amount_perovskite], \
                                       [minerals.Murakami_etal_2012.fe_perovskite(), \
                                            minerals.Murakami_etal_2012.fe_periclase()] )
    
        rock.set_method(method)
    
        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
            burnman.velocities_from_rock(rock, seis_p, temperature, burnman.averaging_schemes.VoigtReussHill())
    
        print "Calculations are done for:"
        rock.debug_print()
    
       #[rho_err,vphi_err,vs_err]=burnman.compare_chifactor([mat_vs,mat_vphi,mat_rho],[seis_vs,seis_vphi,seis_rho])
        [rho_err,vphi_err,vs_err]=burnman.compare_l2(depths,[mat_vs,mat_vphi,mat_rho],[seis_vs,seis_vphi,seis_rho])
    
        return vs_err, vphi_err
    
    xx=np.linspace(0.0, 1.0, 40)
    errs=np.array([material_error(x) for x in xx])
    yy_vs=errs[:,0]
    yy_vphi=errs[:,1]
    plt.plot (xx,yy_vs,"r-x",label=("vs error"))
    plt.plot (xx,yy_vphi,"b-x",label=("vphi error"))
    plt.yscale('log')
    plt.xlabel('% Perovskite')
    plt.ylabel('Error')
    plt.legend()
    plt.show()
