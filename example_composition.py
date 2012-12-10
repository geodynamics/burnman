"""

This example shows how to create different minerals,
how to compute seismic velocities,
and how to compare them to a seismic reference model.

requires:
- geotherms
- seismic models
- compute seismic velocities

teaches:
- creating minerals
- seismic comparison

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

# To compute seismic velocities and other properties, we need to supply
# burnman with a list of minerals (phaes) and their molar abundances. Minerals
# are classes found in burnman.minerals and are derived from
# burnman.minerals.material.
# Here are a few ways to define phases and molar_abundances:

#Example 1: two simple fixed minerals
if False:
	phases = [minerals.Murakami_perovskite(), minerals.Murakami_fp_LS()]
	amount_perovskite = 0.95
	molar_abundances = [amount_perovskite, 1.0-amount_perovskite]

#Example 2: specify fixed iron content
if False:
	phases = [minerals.mg_fe_perovskite(0.7), minerals.ferropericlase(0.5)]
	amount_perovskite = 0.95
	molar_abundances = [amount_perovskite, 1.0-amount_perovskite]

#Example 3: input weight percentages
#See comments in code/composition.py for references to partition coefficent calculation
if True:
	weight_percents = {'Mg':0.213, 'Fe': 0.08, 'Si':0.27, 'Ca':0., 'Al':0.}
	phase_fractions,relative_molar_percent = burnman.calculate_phase_percents(weight_percents)
	iron_content = lambda p,t: burnman.calculate_partition_coefficient(p,t,relative_molar_percent)
	phases = [minerals.mg_fe_perovskite_pt_dependent(iron_content,0), \
			  minerals.ferropericlase_pt_dependent(iron_content,1)]
	molar_abundances = [phase_fractions['pv'],phase_fractions['fp']]

#Example 4: three materials
if False:
	phases = [minerals.Murakami_perovskite(), minerals.ferropericlase(0.5), minerals.stishovite()]
	molar_abundances = [0.7, 0.2, 0.1]


#seismic model for comparison:
seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see code/seismic.py)
number_of_points = 20 #set on how many depth slices the computations should be done
# we will do our computation and comparison at the following depth values:
depths = np.linspace(700, 2800, number_of_points)
#alternatively, we could use the values where prem is defined:
#depths = seismic_model.internal_depth_list()
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

        
geotherm = burnman.geotherm.brown_shankland
temperature = [geotherm(p) for p in seis_p]

for ph in phases:
	ph.set_method(method)

print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i].to_string()

mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = burnman.calculate_velocities(seis_p, temperature, phases, molar_abundances)	

[rho_err,vphi_err,vs_err]=burnman.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)
	

# PLOTTING

plt.subplot(2,2,1)
p1,=plt.plot(seis_p/1.e9,mat_vs,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4,label='computation')
p2,=plt.plot(seis_p/1.e9,seis_vs,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4,label='reference')
plt.title("Vs (km/s)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.ylim(5.1,7.6)
plt.legend(loc='lower right')
plt.text(40,7.3,"misfit= %3.3f" % vs_err)

# plot Vphi
plt.subplot(2,2,2)
p1,=plt.plot(seis_p/1.e9,mat_vphi,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
p2,=plt.plot(seis_p/1.e9,seis_vphi,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("Vphi (km/s)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.ylim(7,12)
#	plt.legend([p1,p2],["Murakami (0.93 Pv, 0.07 fp)", "seismic model (PREM)"], loc=4)
plt.text(40,11.5,"misfit= %3.3f" % vphi_err)

# plot density
plt.subplot(2,2,3)
p1,=plt.plot(seis_p/1.e9,mat_rho,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
p2,=plt.plot(seis_p/1.e9,seis_rho,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("density (kg/m^3)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.text(40,4.3,"misfit= %3.3f" % rho_err)
plt.xlabel("Pressure (GPa)")


# plot geotherm
plt.subplot(2,2,4)
p1,=plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
plt.title("Geotherm (K)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.xlabel("Pressure (GPa)")

	


plt.savefig("example_composition.png")
plt.show()
