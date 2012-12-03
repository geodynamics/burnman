import os, sys, numpy as np
import matplotlib.pyplot as plt


lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import minerals 
from code import main as main

#See bottom for output tool


### input variables ###
#######################

#INPUT for method
method = 'slb' # choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)

#INPUT for geotherm
geotherm = 'geotherm_brown_shankland'


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
 
#input the output filename
output_filename = "output.txt"
 
        
temperature = main.build_geotherm(geotherm, seis_p)

for ph in phases:
	ph.set_method(method)

print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i].to_string()

mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = main.calculate_velocities(seis_p, temperature, phases, molar_abundances)	

[rho_err,vphi_err,vs_err]=main.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho/1e3) #check units on rho
	
#write to file

file = open(output_filename, 'wb')
	file.write("Pressure	Temperature	mat_rho		mat_vs		mat_vp		mat_vphi	mat_K		mat_mu\n")
	writer = csv.writer(file, delimiter="\t")
	data = np.array([seis_p,temperature,mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu]).transpose()
	writer.writerows(data)
	
	print "\nYour data has been saved it: ",output_filename

