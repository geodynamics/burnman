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

def compare_with_seismic_model(vs_err,vphi_err,rho_err,seis_vs,seis_vphi,seis_rho):

	
	for i in range(len(seis_vs)):
		rho_err = (rho_err-seis_rho)
		vp_err = (vp_err/1.e3- seis_vp)
   		vs_err = (vs_err/1.e3- seis_vs)
    		rho_err = pow(rho_err,2.)/pow(seis_rho,2.)
    		vp_err = pow(vp_err,2.)/pow(seis_vp,2.)
    		vs_err = pow(vs_err,2.)/pow(seis_vs,2.)

	rho_err_tot = 100.*integrate.trapz(rho_err)
	vp_err_tot = 100.*integrate.trapz(vp_err)
	vs_err_tot = 100.*integrate.trapz(vs_err)
    	err_tot=rho_err_tot+vp_err_tot+vs_err_tot

	print 'density misfit=',rho_err_tot, 'vp misfit=',vp_err_tot,'vs misfit=', vs_err_tot, 'total=',err_tot

	return rho_err_tot, vp_err_tot, vs_err_tot
