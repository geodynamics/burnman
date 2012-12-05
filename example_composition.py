import os, sys, numpy as np
import matplotlib.pyplot as plt

lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import minerals 
from code import main as main

### input variables ###
#######################

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

#Example 1: simple fixed minerals
if False:
	phases = (minerals.Murakami_perovskite(), minerals.Murakami_fp_LS())
	amount_perovskite = 0.95
	molar_abundances = ( amount_perovskite, 1.0-amount_perovskite)

#Example 2: specify iron content
if True:
	phases = (minerals.mg_fe_perovskite(0.7), minerals.ferropericlase(0.5))
	amount_perovskite = 0.95
	molar_abundances = ( amount_perovskite, 1.0-amount_perovskite)

#Example 3: input weight percentages
#See comments in code/composition.py for references to partition coefficent calculation
if False:
	weight_percents = {'Mg':0.213, 'Fe': 0.0626, 'Si':0.242, 'Ca':0., 'Al':0.}
	phase_fractions,relative_molar_percent = part.conv_inputs(weight_percents)
	iron_content = lambda p,t: part.calculate_partition_coefficient(p,t,relative_molar_percent)
	phases = (minerals.mg_fe_perovskite_pt_dependent(iron_content,0), \
			  minerals.ferropericlase_pt_dependent(iron_content,1))
	molar_abundances = (phase_fractions['pv'],phase_fractions['fp'])

#Example 4: more complicated, three materials
if False:
	phases = (minerals.Murakami_perovskite(), minerals.ferropericlase(0.5), minerals.stishovite())
	molar_abundances = ( 0.7, 0.2, 0.1)


#INPUT for seismic models
name = 'prem'  					# choose from 'prem'(at 1 second, Dziewonski & Anderson, 1981), 'ref_fast' (Lekic et al. 2012), 'ref_slow' (Lekic et al. 2012)
attenuation_correction= 'off'   		# correction for attuation (not required for PREM at 1s, see Matas et al. 2007, page 4) choose from 'on', 'off'
depth_min = 0.0   				# minimum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_max = 0.0   				# maximum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_step = 100.0 				# steps at which the seismic velocity and calculated velocity are compared (in between is interpolated)

seis_p, seis_r, seis_vp, seis_vphi, seis_vs, seis_rho = main.seismo_load(name, attenuation_correction, depth_min, depth_max,depth_step)
 
        
geotherm = main.get_geotherm("brown_shankland")
temperature = [geotherm(p) for p in seis_p]

for ph in phases:
	ph.set_method(method)

print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i].to_string()

mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = main.calculate_velocities(seis_p, temperature, phases, molar_abundances)	

[rho_err,vphi_err,vs_err]=main.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho/1e3) #check units on rho
	

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
p2,=plt.plot(seis_p/1.e9,seis_rho/1.e3,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("density (kg/m^3)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.text(40,5.3,"misfit= %3.3f" % rho_err)

# plot geotherm
plt.subplot(2,2,4)
p1,=plt.plot(seis_p/1e9,temperature,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=4)
plt.title("Geotherm (K)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
	


plt.savefig("example_composition.png")
plt.show()
