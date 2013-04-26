# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
The main focus of this example is to show the mineral physical input constants
necessary for BurnMan to calculate seismic velocity profiles. Furht

Shows user how to input a mineral of his/her choice and which physical values
need to be input for BurnMan to calculate Vs, Vp, Vphi and density at depth.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- how to create your own minerals

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

# A note about units: all the material parameters are expected to be in plain SI units.
# This means that the elastic moduli should be in Pascals and NOT Gigapascals,
# and the Debye temperature should be in K not C.  Additionally, the reference volume
# should be in m^3/(mol molecule) and not in unit cell volume and 'n' should be
# the number of atoms per molecule.  Frequently in the literature the reference volume
# is given in Angstrom^3 per unit cell.  To convert this to m^3/(mol of molecule) 
#you should multiply by 10^(-30) * N_a / Z, where N_a is Avogadro's number and Z is the number of
# atoms per unit cell.  You can look up Z in many places, including www.mindat.org

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
    method = 'slb3'
    
    #in form name_of_mineral (burnman.material <- creates list with parameters)
    class own_material (burnman.material): 
            def __init__(self):
                    burnman.material.__init__(self)
                    self.params = {
                            'ref_V': 10.844e-6, #Molar volume [m^3/(mole molecules)] 
                            	#at room pressure/temperature
                            'ref_K': 135.19e9, #Reference bulk modulus [Pa] 
                            	#at room pressure/temperature
                            'K_prime': 6.04, #pressure derivative of bulk modulus
                            'ref_mu': 175.0e9, #reference shear modulus 
                           	 	#at room pressure/temperature
                            'mu_prime': 1.7, #pressure derivative of shear modulus
                            'molar_mass': .055845, #molar mass in units of [kg/mol]
                            'n': 1, #number of atoms per formula unit
                            'ref_Debye': 998.85, #Debye temperature for material. 
                            	#See Stixrude & Lithgow-Bertelloni, 2005 for values 
                            'ref_grueneisen': 1.368, #Gruneisen parameter for material. 
                            	#See Stixrude & Lithgow-Bertelloni, 2005 for values
                            'q0': 0.917, #isotropic strain derivative of gruneisen
                            	#parameter. Values in Stixrude & Lithgow-Bertelloni, 2005 
       						 'eta_0s': 3.0} #full strain derivative of gruneisen parameter
               					#parameter. Values in Stixrude & Lithgow-Bertelloni, 2005
               					
    
    
    rock = burnman.composite( [(own_material(), 1.0)] )
    
    #seismic model for comparison: (see burnman/seismic.py)
    seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() 
    number_of_points = 20 #set on how many depth slices the computations should be done
    depths = np.linspace(700e3,2800e3, number_of_points)
    #depths = seismic_model.internal_depth_list()
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
    
            
    geotherm = burnman.geotherm.brown_shankland
    temperature = [geotherm(p) for p in seis_p]
    
    rock.set_method(method)
    
    print "Calculations are done for:"
    for ph in rock.phases:
        print ph.fraction, " of phase", ph.mineral.to_string()
    
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = \
        burnman.velocities_from_rock(rock, seis_p, temperature, \
        burnman.averaging_schemes.voigt_reuss_hill())    
    
    [rho_err,vphi_err,vs_err]= \
    burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
