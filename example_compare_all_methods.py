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

if __name__ == "__main__":	
	
	###Input Model 1
	
	#INPUT for method_1
	""" choose 'slb' (finite-strain 2nd order sheer modulus, stixrude and lithgow-bertelloni, 2005)
	or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)
	or 'bm' (birch-murnaghan, if you choose to ignore temperature (your choice in geotherm will not matter in this case))
	or 'slb3 (finite-strain 3rd order shear modulus, stixrude and lithgow-bertelloni, 2005)"""
	
	method_1 = 'mgd' 

        ##input second method
        method_2 = 'slb'

        ##input third method
        method_3 = 'slb3'
        
        ##input fourth method
        method_4 = 'bm'	
	
	#Input composition of model 1. See example_composition for potential choices. We'll just choose something simple here
		
	amount_perovskite = 0.95
	rock = composite( ( (minerals.Murakami_fe_perovskite(), amount_perovskite), 
                            (minerals.Murakami_fe_periclase_LS(), 1.0-amount_perovskite) ) )
	
	#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array
	
	seis_p = np.arange(25e9, 125e9, 5e9)
	
	#input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.
	
	geotherm = burnman.geotherm.brown_shankland
	temperature = [geotherm(p) for p in seis_p]
	
	
	#Now we'll calculate the models. 
	
	rock.set_method(method_1)
	
	print "Calculations are done for:"
	for ph in rock.phases:
		print ph.fraction, " of phase", ph.mineral.to_string()
	
	mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_mu_1 = burnman.calculate_velocities(seis_p, temperature, rock)	
	
	rock.set_method(method_2)
	
	print "Calculations are done for:"
	for ph in rock.phases:
		print ph.fraction, " of phase", ph.mineral.to_string()
	
	mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_mu_2 = burnman.calculate_velocities(seis_p, temperature, rock)	

        rock.set_method(method_3)

        print "Calculations are done for:"
	for ph in rock.phases:
		print ph.fraction, " of phase", ph.mineral.to_string()

        mat_rho_3, mat_vp_3, mat_vs_3, mat_vphi_3, mat_K_3, mat_mu_3 = burnman.calculate_velocities(seis_p, temperature, rock)        

        rock.set_method(method_4)

        print "Calculations are done for:"
	for ph in rock.phases:
		print ph.fraction, " of phase", ph.mineral.to_string()

        mat_rho_4, mat_vp_4, mat_vs_4, mat_vphi_4, mat_K_4, mat_mu_4 = burnman.calculate_velocities(seis_p, temperature, rock)	
	
	##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)
	
	plt.subplot(2,2,1)
	plt.plot(seis_p/1.e9,mat_vs_1,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,mat_vs_2,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
        plt.plot(seis_p/1.e9,mat_vs_3,color='b',linestyle='-',marker='x',markerfacecolor='b',markersize=4)
        plt.plot(seis_p/1.e9,mat_vs_4,color='g',linestyle='-',marker='o',markerfacecolor='g',markersize=4)
	plt.title("Vs (km/s)")
	
	
	# plot Vphi
	plt.subplot(2,2,2)
	plt.plot(seis_p/1.e9,mat_vphi_1,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
	plt.plot(seis_p/1.e9,mat_vphi_2,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
        plt.plot(seis_p/1.e9,mat_vphi_3,color='b',linestyle='-',marker='x',markerfacecolor='b',markersize=4)
        plt.plot(seis_p/1.e9,mat_vphi_4,color='g',linestyle='-',marker='o',markerfacecolor='g',markersize=4)
	plt.title("Vphi (km/s)")
	
	# plot density
	plt.subplot(2,2,3)
	plt.plot(seis_p/1.e9,mat_rho_1,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='mgd')
	plt.plot(seis_p/1.e9,mat_rho_2,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='slb')
        plt.plot(seis_p/1.e9,mat_rho_3,color='b',linestyle='-',marker='x',markerfacecolor='b',markersize=4,label='slb3')
        plt.plot(seis_p/1.e9,mat_rho_4,color='g',linestyle='-',marker='o',markerfacecolor='g',markersize=4,label='bm')
	plt.title("density (kg/m^3)")
	plt.legend(loc='upper left')
	
	
        plt.savefig("compare_all_methods.png")	
	plt.show()
