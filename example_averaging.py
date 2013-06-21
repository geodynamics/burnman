# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:
1. Voight-Reuss-Hill (volumetric averaging)
2. linear (mol fraction averaging)
3. Voight averaging
4. Reuss averaging

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

	rock = burnman.composite( [ (minerals.SLB2005.mg_perovskite(), amount_perovskite), 
				    (minerals.SLB2005.periclase(), 1.0-amount_perovskite) ] )

		   
	#seismic model for comparison:
	# pick from .prem() .slow() .fast() (see burnman/seismic.py)
	seismic_model = burnman.seismic.prem() 
	#set on how many depth slices the computations should be done
	number_of_points = 20
	# we will do our computation and comparison at the following depth values:
	depths = np.linspace(700e3, 2800e3, number_of_points)
	#alternatively, we could use the values where prem is defined:
	#depths = seismic_model.internal_depth_list()
	seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

	temperature = burnman.geotherm.brown_shankland(seis_p)
	
	rock.set_method(method)

	print "Calculations are done for:"
	for ph in rock.phases:
		print ph.fraction, " of phase", ph.mineral.to_string()
	
	#Moduli are calculated at along a P,T profile for each mineral
	moduli_list = burnman.calculate_moduli(rock, seis_p, temperature)

	#Final moduli are then calculated for each of the averaging scheme
	moduli = burnman.average_moduli(moduli_list, \
	burnman.averaging_schemes.voigt_reuss_hill())
	mat_vs1, mat_vp1, mat_vphi1 = burnman.compute_velocities(moduli)
	mat_K1, mat_mu1, mat_rho1 = moduli.K, moduli.mu, moduli.rho

	moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.linear())
	mat_vs2, mat_vp2, mat_vphi2 = burnman.compute_velocities(moduli)
	mat_K2, mat_mu2, mat_rho2 = moduli.K, moduli.mu, moduli.rho

	moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.voigt())
	mat_vs3, mat_vp3, mat_vphi3 = burnman.compute_velocities(moduli)
	mat_K3, mat_mu3, mat_rho3 = moduli.K, moduli.mu, moduli.rho

	moduli = burnman.average_moduli(moduli_list, burnman.averaging_schemes.reuss())
	mat_vs4, mat_vp4, mat_vphi4 = burnman.compute_velocities(moduli)
	mat_K4, mat_mu4, mat_rho4 = moduli.K, moduli.mu, moduli.rho
	
	moduli = moduli_list[0]
	mat_vsa, mat_vpa, mat_vphia = burnman.compute_velocities(moduli)
	mat_Ka, mat_mua, mat_rhoa = moduli.K, moduli.mu, moduli.rho

	moduli = moduli_list[1]
	mat_vsb, mat_vpb, mat_vphib = burnman.compute_velocities(moduli)
	mat_Kb, mat_mub, mat_rhob = moduli.K, moduli.mu, moduli.rho
		
	
	# PLOTTING
	
	# plot vs
	fig=plt.figure()
	plt.plot(seis_p/1.e9,mat_vs3/1.e3,color='c',linestyle='-',marker='^',\
	markersize=4,label='voigt')
	plt.plot(seis_p/1.e9,mat_vs4/1.e3,color='k',linestyle='-',marker='v',\
	markersize=4,label='reuss')
	plt.plot(seis_p/1.e9,mat_vs1/1.e3,color='b',linestyle='-',marker='x',\
	markersize=4,label='VRH')
	plt.plot(seis_p/1.e9,mat_vs2/1.e3,color='r',linestyle='--',marker='x',\
	markersize=4,label='linear')
	plt.plot(seis_p/1.e9,mat_vsa/1.e3,color='y',linestyle='-',marker='x',\
	markersize=4,label='mg perovskite')
	plt.plot(seis_p/1.e9,mat_vsb/1.e3,color='g',linestyle='-',marker='x',\
	markersize=4,label='periclase')
	#plt.title("Vs (km/s), mix: 60% mg perovskite")
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	#plt.ylim(5.2,7.0)
	plt.legend(loc='upper left',prop={'size':11},frameon=False)
	plt.xlabel('pressure (GPa)')
	plt.ylabel('Vs (km/s)')
	#plt.savefig("output_figures/example_averaging.png")
	#plt.show()

	mat_vsa_norm=(mat_vsa-mat_vsb)/(mat_vsa-mat_vsb)
	mat_vsb_norm=(mat_vsb-mat_vsb)/(mat_vsa-mat_vsb)
	mat_vs1_norm=(mat_vs1-mat_vsb)/(mat_vsa-mat_vsb)
	mat_vs2_norm=(mat_vs2-mat_vsb)/(mat_vsa-mat_vsb)
	mat_vs3_norm=(mat_vs3-mat_vsb)/(mat_vsa-mat_vsb)
	mat_vs4_norm=(mat_vs4-mat_vsb)/(mat_vsa-mat_vsb)
	ax=fig.add_axes([0.58, 0.18, 0.3, 0.3])
	plt.plot(seis_p/1.e9,mat_vs3_norm,color='c',linestyle='-',marker='^',\
	markersize=4,label='Voigt')
	plt.plot(seis_p/1.e9,mat_vs4_norm,color='k',linestyle='-',marker='v',\
	markersize=4,label='Reuss')
	plt.plot(seis_p/1.e9,mat_vs1_norm,color='b',linestyle='-',marker='x',\
	markersize=4,label='VRH')
	plt.plot(seis_p/1.e9,mat_vs2_norm,color='r',linestyle='--',marker='',\
	markersize=4,label='linear')
	plt.plot(seis_p/1.e9,mat_vsa_norm,color='y',linestyle='-',marker='x',\
	markersize=4,label='perovskite')
	plt.plot(seis_p/1.e9,mat_vsb_norm,color='g',linestyle='-',marker='x',\
	markersize=4,label='periclase')
	ax.tick_params(labelsize=10)
	plt.title("normalized by mixture endmembers",fontsize=10)
	plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	plt.ylim(-0.005,1.005)
	plt.xlabel('pressure (GPa)',fontsize=10)
	plt.ylabel('normalized Vs',fontsize=10)
	#plt.legend(loc='lower right')

	plt.savefig("output_figures/example_averaging_normalized.png")
	plt.show()
