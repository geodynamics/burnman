import os, sys, numpy as np
import matplotlib.pyplot as plt

import birch_murnaghan
import comparison as comp
import composition as part
import geotherm as gt
import mie_grueneisen_debye
import minerals 
import seismic
import slb_finitestrain
import tools
import voigt_reuss_hill as vrh
import scipy.integrate as integrate


def calculate_velocities (pressure, temperature, phases, molar_abundances):

	mat_vs = np.empty_like(pressure)
	mat_vp = np.empty_like(pressure)
	mat_vphi = np.empty_like(pressure)
	mat_rho = np.empty_like(pressure)
	mat_K = np.empty_like(pressure)
	mat_mu = np.empty_like(pressure)
	#print "Calculating elastic properties for phase assemblage \n"
	#print "seismic p (GPa)	T (K)	density(kg/m^3)	K(Gpa) G(GPa)	Vs (km/s)	Vp(km/s)	Vphi (km/s)"
	for i in range(len(pressure)):
		rho,vp,vs,vphi,K,mu = \
		vrh.voigt_reuss_hill(pressure[i], temperature[i], phases, molar_abundances)
		#print pressure[i]/1.e9,"	", temperature[i],"	",rho,"	", K,"	", mu,"	", vs/1.e3,"	", vp/1.e3,"	", vphi/1.e3
		#print pressure[i]/1.e9,"	",rho,"	", mu
		mat_rho[i] = rho/1.e3 # convert from kg/m^3 to g/cc
		mat_vs[i] = vs/1.e3
		mat_vp[i] = vp/1.e3
    		mat_vphi[i] = vphi/1.e3
   		mat_K[i] = K
		mat_mu[i] = mu

	return mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_mu

# apply attenuation correction to a list of v_p, v_s, v_phi
# takes a list of pressures and the seismic model
def apply_attenuation_correction(seismic_model,pressure,v_p,v_s,v_phi):
	length = len(pressure)
	ret_v_p = np.zeros(length)
	ret_v_s = np.zeros(length)
	ret_v_phi = np.zeros(length)
	for i in range(length):
		ret_v_p[i],ret_v_s[i],ret_v_phi[i] = \
		    vrh.attenuation_correction(seismic_model, pressure, v_p[i], v_s[i], v_phi[i])
	
	return ret_v_p, ret_v_s, ret_v_phi


# returns the differences of the two velocity profiles sampled at the given depth
# it computes the L2 norm of the two functions (assumed to be linear between points)
def compare_two(depth,mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho):

	rho_err_tot = l2(depth,mat_rho,seis_rho)
	vphi_err_tot = l2(depth,mat_vphi,seis_vphi)
	vs_err_tot = l2(depth,mat_vs,seis_vs)
    	err_tot=rho_err_tot+vphi_err_tot+vs_err_tot

	return rho_err_tot, vphi_err_tot, vs_err_tot

# weighted comparison
def compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho):

	rho_err_tot = chi_factor(mat_rho,seis_rho)
	vphi_err_tot = chi_factor(mat_vphi,seis_vphi)
	vs_err_tot = chi_factor(mat_vs,seis_vs)
    	err_tot=rho_err_tot+vphi_err_tot+vs_err_tot

	return rho_err_tot, vphi_err_tot, vs_err_tot

def l2(x,funca,funcb):
	diff=np.array(funca-funcb)
	diff=diff*diff
	length=x[-1]-x[0]
	assert(length>0)
	return integrate.trapz(diff,x) / length
	

def chi_factor(calc,obs):
	#assuming 1% a priori uncertainty on the seismic model

	err=np.empty_like(calc)
	for i in range(len(calc)):
		err[i]=pow((calc[i]-obs[i])/(0.01*np.mean(obs)),2.)

	err_tot=np.sum(err)/len(err)

	return err_tot


def get_geotherm(name):
	if (name=='brown_shankland'):
		return lambda p: gt.geotherm_brown_shankland(p)
	elif (name=='watson_baxter'):
		return lambda p: gt.geotherm_watson_baxter(p)
	else:
		raise ValueError, "Error, unknown geotherm '" + name + "'"
	return []
