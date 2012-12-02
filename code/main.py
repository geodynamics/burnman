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

def calculate_velocities (pressure, temperature, phases, molar_abundances):
	mat_vs = np.empty_like(pressure)
	mat_vp = np.empty_like(pressure)
	mat_vphi = np.empty_like(pressure)
	mat_rho = np.empty_like(pressure)
	mat_K = np.empty_like(pressure)
	mat_mu = np.empty_like(pressure)
	print "Calculating elastic properties for phase assemblage \n"
	print "seismic p (GPa)	T (K)	density(kg/m^3)	K(Gpa) G(GPa)	Vs (km/s)	Vp(km/s)	Vphi (km/s)"
	for i in range(len(pressure)):
		rho,vp,vs,vphi,K,mu = \
		vrh.voigt_reuss_hill(pressure[i], temperature[i], phases, molar_abundances)
		print pressure[i]/1.e9,"	", temperature[i],"	",rho,"	", K,"	", mu,"	", vs/1.e3,"	", vp/1.e3,"	", vphi/1.e3
		#print pressure[i]/1.e9,"	",rho,"	", mu
		mat_rho[i] = rho/1.e3
		mat_vs[i] = vs/1.e3
		mat_vp[i] = vp/1.e3
    		mat_vphi[i] = vphi/1.e3
   		mat_K[i] = K
		mat_mu[i] = mu

	return mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu

def compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho):


	rho_err_tot = chi_factor(mat_rho,seis_rho)
	vphi_err_tot = chi_factor(mat_vphi,seis_vphi)
	vs_err_tot = chi_factor(mat_vs,seis_vs)
    	err_tot=rho_err_tot+vphi_err_tot+vs_err_tot

	return rho_err_tot, vphi_err_tot, vs_err_tot

def chi_factor(calc,obs):
	#assuming 1% a priori uncertainty on the seismic model

	err=np.empty_like(calc)
	for i in range(len(calc)):
		err[i]=pow((calc[i]-obs[i])/(0.01*np.mean(obs)),2.)

	err_tot=np.sum(err)/len(err)

	return err_tot

def seismo_load(name,attenuation_correction,depth_min,depth_max,depth_step):
	s=seis.seismo_in()
	s.params['name']=name
	s.params['attenuation_correction']=attenuation_correction
	s.params['depth_min']=depth_min
	s.params['depth_max']=depth_max
	s.params['depth_step']=depth_step

	print "Looking up seismic model:",s.params['name']
	seis_p = np.array(s.pressure())*1.e9     
	seis_r = np.array(s.radius())/1.e3
	seis_vp = np.array(s.v_p())/1.e3#[seis.seis_V(y/1.e9)[0]/1.e3 for y in seis_p]
	seis_vphi = np.array(s.v_phi())/1.e3#
	seis_vs = np.array(s.v_s())/1.e3#[seis.seis_V(y/1.e9)[1]/1.e3 for y in seis_p]
	seis_rho = np.array(s.density())#[seis.seis_density(y/1.e9) for y in seis_p]
	
	return seis_p, seis_r, seis_vp, seis_vphi, seis_vs, seis_rho


def build_geotherm(name,pressure):
	if (name=='geotherm_brown_shankland'):
		geotherm = lambda p: gt.geotherm_brown_shankland(p)
		print "geotherm is from Brown & Shankland (1981)"
	else:
		print "define geothermal gradient"
        
	temperature = [geotherm(p) for p in pressure]
	
	return temperature

def phase_builder(kind,

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
