# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Compute properties of minerals and creates a table of outputs in a text format
that could be used with other programs.

requires:
- creating minerals
- compute seismic velocities
- geotherms

teaches:
- output computed seismic data to file 

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":    
    
    #See bottom for output tool
    
    
    ### input variables ###
    #######################
    
    #INPUT for method
    """ choose 'slb2' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus, matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shearl modulus, matas et al. 2007)
    or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))"""
    method = 'mgd3' 
    
    #specify material
    amount_perovskite = 0.95
    rock = burnman.composite( ((minerals.mg_fe_perovskite(0.7), amount_perovskite), 
                               (minerals.ferropericlase(0.5), 1.0-amount_perovskite) ) )
    
    #define some pressure range
    pressures = np.arange(25e9,130e9,5e9)
    
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in pressures]
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(pressures, temperature, rock)    
        
    #write to file:
    output_filename = "example_woutput.txt" 
    f = open(output_filename, 'wb')
    f.write("#Pressure\tTemperature\tmat_rho\tmat_vs\tmat_vp\tmat_vphi\tmat_K\tmat_mu\n")
    
    data = zip(pressures,temperature,mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu)
    np.savetxt(f, data, fmt='%.10e', delimiter='\t')
    
        
    print "\nYour data has been saved as: ",output_filename
    
