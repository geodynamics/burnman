#standard numpy, scipy imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

#eos imports
import seismo_in as seis
import mie_grueneisen_debye as mgd
import comparison as comp
import geotherm as gt
import voigt_reuss_hill as vrh

def compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho):


	rho_err_tot = chi_factor(mat_rho,seis_rho)
	vphi_err_tot = chi_factor(mat_vphi,seis_vphi)
	vs_err_tot = chi_factor(mat_vs,seis_vs)
    	err_tot=rho_err_tot+vphi_err_tot+vs_err_tot

	print 'density misfit=',rho_err_tot, 'vphi misfit=',vphi_err_tot,'vs misfit=', vs_err_tot, 'total=',err_tot

	return rho_err_tot, vphi_err_tot, vs_err_tot


def madeup_misfit(calc,obs):

	
	err = np.empty_like(calc)
	for i in range(len(calc)):
		err[i] = (calc[i]-obs[i])
		err[i] = pow(err[i],2.)/pow(obs[i],2.)
	
	err_tot=100.*integrate.trapz(err)

	return err_tot

def chi_factor(calc,obs):
	#assuming 1% a priori uncertainty on the seismic model

	err=np.empty_like(calc)
	for i in range(len(calc)):
		err[i]=pow((calc[i]-obs[i])/(0.01*np.mean(obs)),2.)

	err_tot=np.sum(err)/len(err)

	return err_tot
