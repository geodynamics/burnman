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

def compare_with_seismic_model(phases,molar_abundances,geoth):

        s=seis.seismo_in()
    	seis_p = np.array(s.pressure())*1.e9           #np.arange(26.e9,135.0e9,1e9)
    	seis_r = np.array(s.radius())/1.e3
    	seis_vp = np.array(s.v_p())/1.e3#[seis.seis_V(y/1.e9)[0]/1.e3 for y in seis_p]
    	seis_vs = np.array(s.v_s())/1.e3#[seis.seis_V(y/1.e9)[1]/1.e3 for y in seis_p]
    	seis_density = np.array(s.density())#[seis.seis_density(y/1.e9) for y in seis_p]
        #	geoth = lambda p: sum(molar_abundances[i]*phases[i].geotherm(p) for i in range(len(phases)))
	# Pearson's chi squared test	
	rho_err = np.empty_like(seis_p)
	vs_err = np.empty_like(seis_p)
	vp_err = np.empty_like(seis_p)
	
	for i in range(len(seis_p)):
		rho_err[i], vp_err[i], vs_err[i],v_phi,K,mu = \
			vrh.voigt_reuss_hill(seis_p[i], geoth(seis_p[i]), phases, molar_abundances)
	rho_err = (rho_err-seis_density)
	vp_err = (vp_err/1.e3- seis_vp)
   	vs_err = (vs_err/1.e3- seis_vs)
    	rho_err = pow(rho_err,2.)/pow(seis_density,2.)
    	vp_err = pow(vp_err,2.)/pow(seis_vp,2.)
    	vs_err = pow(vs_err,2.)/pow(seis_vs,2.)

	rho_err_tot = 100.*integrate.trapz(rho_err)
	vp_err_tot = 100.*integrate.trapz(vp_err)
	vs_err_tot = 100.*integrate.trapz(vs_err)
    	err_tot=rho_err_tot+vp_err_tot+vs_err_tot

	print 'density misfit=',rho_err_tot, 'vp misfit=',vp_err_tot,'vs misfit=', vs_err_tot, 'total=',err_tot

	return rho_err_tot, vp_err_tot, vs_err_tot
