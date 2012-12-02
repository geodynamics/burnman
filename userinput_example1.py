import os, sys, numpy as np
import matplotlib.pyplot as plt

lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import birch_murnaghan
from code import comparison as comp
from code import composition as part
from code import geotherm as gt
from code import material
from code import mie_grueneisen_debye
from code import minerals 
from code import prem
from code import seismo_in as seis
from code import slb_finitestrain
from code import tools
from code import voigt_reuss_hill as vrh
from code import main as main

### input variables ###
#######################


#INPUT for method
method = 'mgd'					# choose 'slb' (finite-strain, stixrude and lithgow-bertelloni, 2005, ONLY ONE WORKING) or 'mgd' (mie-gruneisen-debeye, matas et al. 2007)


#INPUT for geotherm
geotherm = 'geotherm_brown_shankland'



#INPUT for composition
# for material constants see minerals.py
pv = 'minerals.Murakami_perovskite'            	# choose 'mg_fe_perovskite' or 'Murakami_perovskite'
fp = 'minerals.Murakami_fp_LS'                   	# choose 'ferropericlase', or 'Murakami_fp_LS'

# for material constants see minerals.py
composition_input = '2phase_fractions'		# choose 'weight_percents', '2phase_fractions', or 'nphase_fractions'

# if '2phase_fractions' or 'nphase_fractions'
phase_names=('pv','fp')
phases = (pv,fp)
phase_fractions = {'pv':0.95, 'fp':0.05} 	# should add up to 1.0

# if '2phase_fractions'
calculate_partitioning = 'off'          	# sets if partioning coeefficients are used or calculated, 'on'/'off'/'auto', only for 2 phase fractions
partitioning = {'pv':0.94, 'fp':0.79}   	# ignored if calculate_partioning = 'on'artitioning = {'pv':0.94, 'fp':0.79}   # ignored if calculate_partioning = 'on'

#INPUT for seismic models
name = 'prem'  					# choose from 'prem'(at 1 second, Dziewonski & Anderson, 1981), 'ref_fast' (Lekic et al. 2012), 'ref_slow' (Lekic et al. 2012)
attenuation_correction= 'off'   		# correction for attuation (not required for PREM at 1s, see Matas et al. 2007, page 4) choose from 'on', 'off'
depth_min = 0.0   				# minimum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_max = 0.0   				# maximum depth considered in km, choose 0.0 to use the limits given by the seismic model
depth_step = 100.0 				# steps at which the seismic velocity and calculated velocity are compared (in between is interpolated)

seis_p, seis_r, seis_vp, seis_vphi, seis_vs, seis_rho = main.seismo_load(name, attenuation_correction, depth_min, depth_max,depth_step)
 
        
temperature = main.build_geotherm(geotherm, seis_p)


if (composition_input=='weight_percents'):
	weight_percents = weight_percents 
	phase_fractions,relative_molar_percent= part.conv_inputs(weight_percents)
	iron_content = lambda p,t: part.calculate_partition_coefficient(p,t,relative_molar_percent)
	pv_d = mg_fe_perovskite_pt_dependent(iron_content)
	fp_d = ferropericlase_pt_dependent(iron_content)
        phases = (pv_d, fp_d)
	iron_content=("(P,T) dependent","(P,T) dependent")	
	molar_abundances = (phase_fractions['pv'],phase_fractions['fp'])
elif (composition_input=='2phase_fractions'):
	if(calculate_partitioning=='auto'):
		iron_content = lambda p,t: part.calculate_partition_coefficient(p,t,relative_molar_percent)
		pv_d = mg_fe_perovskite_pt_dependent(iron_content)
		fp_d = ferropericlase_pt_dependent(iron_content)
        	pv = eval(pv)(pv_d)
        	fp = eval(fp)(fp_d)
        elif (calculate_partitioning=='on'):
		pv_d=partitioning['pv']
		fp_d=partitioning['fp']
		iron_content=(pv_d,fp_d)
        	pv = eval(pv)(pv_d)
        	fp = eval(fp)(fp_d)
	else:
		pv = eval(pv)()
  		fp = eval(fp)()
		iron_content=('not used', 'not used')
        phases = (pv, fp)
	molar_abundances = (phase_fractions['pv'],	phase_fractions['fp'])
elif (composition_input=='nphase_fractions'):
        list=phase_names
        molar_abundances=np.empty(len(list))
        iron_content=np.empty(len(list))
        phases=[]
	for i in range(len(list)):
		try:
			molar_abundances[i]=float(phase_fractions[list[i]]) 
		except:
			molar_abundances[i]=0.
		try:	
			iron_content[i]=partitioning[list[i]]
		except:
			iron_content[i]=None
		print iron_content[i]
		if np.isnan(iron_content[i]):
			try: 
			        phase=eval(str(phases[i]))()
			        print phase
			except:
				raise Exception("Iron partitioning needs to be defined for", userinput.phases[i])
		else:
			try:
                		phase=eval(str(phases[i]))(iron_content[i])
                	except:
				raise Exception(phases[i], "does not handle an iron_partioning #")		
								
		phases.append(phase)
else:
	raise("Choose method to determine composition: 'weight_percents' or 'phase_fractions'")

#set method
for ph in phases:
	ph.set_method(method)


print "Calculations are done for:"
for i in range(len(phases)):
	print molar_abundances[i], " of phase", phases[i], " with partitioning fraction", iron_content[i]

mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu = main.calculate_velocities(seis_p, temperature, phases, molar_abundances)
	

[rho_err,vphi_err,vs_err]=main.compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho/1e3) #check units on rho
	


# read in composition


plt.subplot(2,2,1)
p1,=plt.plot(seis_p/1.e9,mat_vs,color='b',linestyle='-',marker='o',markerfacecolor='b',markersize=4)
p2,=plt.plot(seis_p/1.e9,seis_vs,color='k',linestyle='-',marker='o',markerfacecolor='k',markersize=4)
plt.title("Vs (km/s)")
plt.xlim(min(seis_p)/1.e9,max(seis_p)/1.e9)
plt.ylim(5.1,7.6)
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
	


plt.savefig("output_plots.png")
plt.show()
