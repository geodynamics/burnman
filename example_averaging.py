# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
example_averaging
-----------------

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:

1. Voight-Reuss-Hill
2. Voight averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See Watt et al., 1976 Journal of Geophysics and Space Physics for explanations
of each averaging scheme.

requires:
- geotherms
- compute seismic velocities

teaches:
- averaging

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

if __name__ == "__main__":
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
	
	amount_perovskite = 0.6
	method = 'slb3'

	rock = burnman.Composite( [amount_perovskite, 1.0-amount_perovskite], \
					  [minerals.SLB_2005.mg_perovskite(), \
						   minerals.SLB_2005.periclase()] )
	rock.set_method(method)

	perovskitite = minerals.SLB_2005.mg_perovskite()
	perovskitite.set_method(method)

	periclasite = minerals.SLB_2005.periclase()
	periclasite.set_method(method)
		   
	#seismic model for comparison:
	# pick from .prem() .slow() .fast() (see burnman/seismic.py)
	seismic_model = burnman.seismic.PREM()
	#set on how many depth slices the computations should be done
	number_of_points = 20
	# we will do our computation and comparison at the following depth values:
	depths = np.linspace(700e3, 2800e3, number_of_points)
	#alternatively, we could use the values where prem is defined:
	#depths = seismic_model.internal_depth_list()
	pressures, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

	temperatures = burnman.geotherm.brown_shankland(pressures)
	

	print "Calculations are done for:"
	rock.debug_print()

        #calculate the seismic velocities of the rock using a whole battery of averaging schemes:

        # do the end members, here averaging scheme does not matter (though it defaults to Voigt-Reuss-Hill)
	rho_pv, vp_pv, vs_pv, vphi_pv, K_pv, G_pv = \
            burnman.velocities_from_rock(perovskitite, pressures, temperatures)
	rho_fp, vp_fp, vs_fp, vphi_fp, K_fp, G_fp = \
            burnman.velocities_from_rock(periclasite, pressures, temperatures)

        #Voigt Reuss Hill averaging
	rho_vrh, vp_vrh, vs_vrh, vphi_vrh, K_vrh, G_vrh = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.VoigtReussHill())

        #Voigt averaging
	rho_v, vp_v, vs_v, vphi_v, K_v, G_v = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.Voigt())

        #Reuss averaging
	rho_r, vp_r, vs_r, vphi_r, K_r, G_r = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.Reuss())

        #Upper bound for Hashin-Shtrikman averaging
	rho_hsu, vp_hsu, vs_hsu, vphi_hsu, K_hsu, G_hsu = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.HashinShtrikmanUpper())

        #Lower bound for Hashin-Shtrikman averaging
	rho_hsl, vp_hsl, vs_hsl, vphi_hsl, K_hsl, G_hsl = \
            burnman.velocities_from_rock(rock, pressures, temperatures, averaging_scheme=burnman.averaging_schemes.HashinShtrikmanLower())

	
	# PLOTTING
	
	# plot vs
	fig=plt.figure()

	plt.plot(pressures/1.e9,vs_v/1.e3,color='c',linestyle='-',marker='^',\
	    markersize=4,label='Voigt')
	plt.plot(pressures/1.e9,vs_r/1.e3,color='k',linestyle='-',marker='v',\
	    markersize=4,label='Reuss')
	plt.plot(pressures/1.e9,vs_vrh/1.e3,color='b',linestyle='-',marker='x',\
	    markersize=4,label='Voigt-Reuss-Hill')
	plt.plot(pressures/1.e9,vs_hsu/1.e3,color='r',linestyle='-',marker='x',\
	    markersize=4,label='Hashin-Shtrikman')
	plt.plot(pressures/1.e9,vs_hsl/1.e3,color='r',linestyle='-',marker='x',\
	    markersize=4)
	plt.plot(pressures/1.e9,vs_pv/1.e3,color='y',linestyle='-',marker='x',\
	    markersize=4,label='Mg Perovskite')
	plt.plot(pressures/1.e9,vs_fp/1.e3,color='g',linestyle='-',marker='x',\
	    markersize=4,label='Periclase')
	plt.xlim(min(pressures)/1.e9,max(pressures)/1.e9)
	plt.legend(loc='upper left',prop={'size':11},frameon=False)
	plt.xlabel('pressure (GPa)')
	plt.ylabel('Vs (km/s)')

	vs_pv_norm=(vs_fp-vs_fp)/(vs_pv-vs_fp)
	vs_fp_norm=(vs_fp-vs_fp)/(vs_pv-vs_fp)
	vs_vrh_norm=(vs_vrh-vs_fp)/(vs_pv-vs_fp)
	vs_v_norm=(vs_v-vs_fp)/(vs_pv-vs_fp)
	vs_r_norm=(vs_r-vs_fp)/(vs_pv-vs_fp)
	vs_hsu_norm=(vs_hsu-vs_fp)/(vs_pv-vs_fp)
	vs_hsl_norm=(vs_hsl-vs_fp)/(vs_pv-vs_fp)

	ax=fig.add_axes([0.58, 0.18, 0.3, 0.3])
	plt.plot(pressures/1.e9,vs_v_norm,color='c',linestyle='-',marker='^',\
	markersize=4,label='Voigt')
	plt.plot(pressures/1.e9,vs_r_norm,color='k',linestyle='-',marker='v',\
	markersize=4,label='Reuss')
	plt.plot(pressures/1.e9,vs_vrh_norm,color='b',linestyle='-',marker='x',\
	markersize=4,label='Voigt-Reuss-Hill')
	plt.plot(pressures/1.e9,vs_hsl_norm,color='r',linestyle='-',marker='x',\
	markersize=4,label='Hashin-Shtrikman')
	plt.plot(pressures/1.e9,vs_hsu_norm,color='r',linestyle='-',marker='x',\
	markersize=4)
	plt.plot(pressures/1.e9,vs_pv_norm,color='y',linestyle='-',marker='x',\
	markersize=4,label='Mg Perovskite')
	plt.plot(pressures/1.e9,vs_fp_norm,color='g',linestyle='-',marker='x',\
	markersize=4,label='Periclase')
	ax.tick_params(labelsize=10)
	plt.title("normalized by mixture endmembers",fontsize=10)
	plt.xlim(min(pressures)/1.e9,max(pressures)/1.e9)
	plt.ylim(-0.005,1.005)
	plt.xlabel('pressure (GPa)',fontsize=10)
	plt.ylabel('normalized Vs',fontsize=10)
	#plt.legend(loc='lower right')

	plt.savefig("output_figures/example_averaging_normalized.png")
	plt.show()
