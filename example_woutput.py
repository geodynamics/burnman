# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example explains how to perform the basic i/o of BurnMan. A method of
calculation is chosen, a composite mineral/material (see
example_composition.py for explanation of this process) is created in the
class "rock," finally a geotherm is created and seismic velocities calculated.

Post-calculation, the results are written to a simple text file to
plot/manipulate at the user's whim.

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
    method = 'mgd3' 
    
    #specify material
    amount_perovskite = 0.95
    rock = burnman.composite( [(minerals.SLB_2005.mg_fe_perovskite(0.7), amount_perovskite), 
                               (minerals.SLB_2005.ferropericlase(0.5), 1.0-amount_perovskite) ] )
    
    #define some pressure range
    pressures = np.arange(25e9,130e9,5e9)
    
    temperature = burnman.geotherm.brown_shankland(pressures)
    
    rock.set_method(method) #append method of calculation to suite of minerals chosen
    
    #Begin calculating velocities and density as depth
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
        burnman.velocities_from_rock(rock, pressures, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())
        
    #write to file:
    output_filename = "example_woutput.txt" 
    f = open(output_filename, 'wb')
    f.write("#Pressure\tTemperature\tmat_rho\tmat_vs\tmat_vp\tmat_vphi\tmat_K\tmat_G\n")
    
    data = zip(pressures,temperature,mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_G)
    np.savetxt(f, data, fmt='%.10e', delimiter='\t')
    
        
    print "\nYour data has been saved as: ",output_filename
    
