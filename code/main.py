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
from code import seismic
from code import slb_finitestrain
from code import tools
from code import voigt_reuss_hill as vrh



def calculate_velocities (pressure, temperature, phases, molar_abundances):
	attenuation = False # TODO, make parameter?

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
		vrh.voigt_reuss_hill(pressure[i], temperature[i], phases, molar_abundances, attenuation)
		#print pressure[i]/1.e9,"	", temperature[i],"	",rho,"	", K,"	", mu,"	", vs/1.e3,"	", vp/1.e3,"	", vphi/1.e3
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


def get_geotherm(name):
	if (name=='brown_shankland'):
		return lambda p: gt.geotherm_brown_shankland(p)
	elif (name=='watson_baxter'):
		return lambda p: gt.geotherm_watson_baxter(p)
	else:
		raise ValueError, "Error, unknown geotherm '" + name + "'"
	return []
