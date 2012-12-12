"""

This example shows you how to create two materials from wt% determines the
optimum mixing between the two to match the seismic model of your choice.

requires:
- geotherms
- seismic models
- compute seismic velocities
- creating minerals

teaches:
- weight percent materials

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals


###Input Model 1

##input method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)


#Input composition of model 1. See example_composition for potential choices. We'll just choose something simple here
	
weight_percents_pyro = {'Mg':0.228, 'Fe': 0.0626, 'Si':0.21, 'Ca':0., 'Al':0.} #From Mcdonough 2003
phase_fractions_pyro,relative_molar_percent_pyro = burnman.calculate_phase_percents(weight_percents_pyro)
iron_content = lambda p,t: burnman.calculate_partition_coefficient(p,t,relative_molar_percent_pyro)
phases_pyro = [minerals.mg_fe_perovskite_pt_dependent(iron_content,0), \
	minerals.ferropericlase_pt_dependent(iron_content,1)]
molar_abundances_pyro = [phase_fractions_pyro['pv'], phase_fractions_pyro['fp']]

#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

seis_p_1 = np.arange(25e9, 125e9, 5e9)

#input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

geotherm = burnman.geotherm.brown_shankland
temperature_1 = [geotherm(p) for p in seis_p_1]

##Now onto the second model parameters

##input second method

#Input composition of model enstatite chondrites. See example_composition for potential choices.
	
weight_percents_enst = {'Mg':0.213, 'Fe': 0.0721, 'Si':0.242, 'Ca':0., 'Al':0.} #Javoy 2009 Table 6 PLoM
phase_fractions_enst,relative_molar_percent_enst = burnman.calculate_phase_percents(weight_percents_enst)
iron_content = lambda p,t: burnman.calculate_partition_coefficient(p,t,relative_molar_percent_enst)
phases_enst = [minerals.mg_fe_perovskite_pt_dependent(iron_content,0), \
	minerals.ferropericlase_pt_dependent(iron_content,1)]
molar_abundances_enst = [phase_fractions_enst['pv'], phase_fractions_enst['fp']]

#input pressure range for first model. This could be from a seismic model or something you create. For this example we will create an array

seis_p_2 = np.arange(25e9, 125e9, 5e9)

#input your geotherm. Either choose one (See example_geotherms.py) or create one.We'll use Brown and Shankland.

geotherm = burnman.geotherm.brown_shankland
temperature_2 = [geotherm(p) for p in seis_p_2]


#Now we'll calculate the models. 

for ph in phases_pyro:
	ph.set_method(method)



mat_rho_pyro, mat_vs_pyro, mat_vp_pyro, mat_vphi_pyro, mat_K_pyro, mat_mu_pyro = burnman.calculate_velocities(seis_p_1, temperature_1, phases_pyro, molar_abundances_pyro)	

print "Calculations are done for:"
for i in range(len(phases_pyro)):
	print molar_abundances_pyro[i], " of phase", phases_pyro[i].to_string()

for ph in phases_enst:
	ph.set_method(method)

print "Calculations are done for:"
for i in range(len(phases_enst)):
	print molar_abundances_enst[i], " of phase", phases_enst[i].to_string()

mat_rho_enst, mat_vs_enst, mat_vp_enst, mat_vphi_enst, mat_K_enst, mat_mu_enst = burnman.calculate_velocities(seis_p_2, temperature_2, phases_enst, molar_abundances_enst)	


##let's create PREM for reference
s=burnman.seismic.prem()
depths = map(s.depth, seis_p_1) 
pressures, rho_prem, vp_prem, vs_prem, v_phi_prem = s.evaluate_all_at(depths)


##Now let's plot the comparison. You can conversely just output to a data file (see example_woutput.py)

plt.subplot(2,2,1)
plt.plot(seis_p_1/1.e9,mat_vs_pyro,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
plt.plot(seis_p_2/1.e9,mat_vs_enst,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
plt.plot(seis_p_1/1.e9,vs_prem,color='k',linestyle='-',marker='x',markerfacecolor='k',markersize=4)
plt.title("Vs (km/s)")

# plot Vphi
plt.subplot(2,2,2)
plt.plot(seis_p_1/1.e9,mat_vphi_pyro,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
plt.plot(seis_p_2/1.e9,mat_vphi_enst,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
plt.plot(seis_p_1/1.e9,v_phi_prem,color='k',linestyle='-',marker='x',markerfacecolor='k',markersize=4)
plt.title("Vphi (km/s)")

# plot density
plt.subplot(2,2,3)
plt.plot(seis_p_1/1.e9,mat_rho_pyro,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4,label="C-chondrite")
plt.plot(seis_p_2/1.e9,mat_rho_enst,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4,label="enstatite")
plt.plot(seis_p_1/1.e9,rho_prem,color='k',linestyle='-',marker='x',markerfacecolor='k',markersize=4,label="PREM")
plt.title("density (kg/m^3)")
plt.legend(loc='lower right')
plt.xlabel("Pressure (GPa)")


#plot percent differences
mat_vs_enst_interp = np.interp(seis_p_1, seis_p_2,mat_vs_enst)
mat_vphi_enst_interp = np.interp(seis_p_1, seis_p_2,mat_vphi_enst)
mat_rho_enst_interp = np.interp(seis_p_1, seis_p_2,mat_rho_enst)

per_diff_vs = 100*(mat_vs_pyro - mat_vs_enst_interp)/mat_vs_pyro
per_diff_vphi = 100*(mat_vphi_pyro - mat_vphi_enst_interp)/mat_rho_pyro
per_diff_rho = 100*(mat_rho_pyro - mat_rho_enst_interp)/mat_rho_pyro

plt.subplot(2,2,4)
plt.plot(seis_p_1/1.e9,per_diff_vs,color='g',linestyle='-',marker='o',markerfacecolor='k',markersize=4,label="vs")
plt.plot(seis_p_1/1.e9,per_diff_vphi,color='c',linestyle='-',marker='+',markerfacecolor='k',markersize=4,label="vphi")
plt.plot(seis_p_1/1.e9,per_diff_rho,color='m',linestyle='-',marker='x',markerfacecolor='k',markersize=4,label="density")

plt.title("percent difference")
plt.legend(loc='center right')
plt.xlabel("Pressure (GPa)")

plt.show()
