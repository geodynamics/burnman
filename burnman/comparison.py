"""
    BurnMan- a lower mantle toolkit
    Copyright (C) 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def compare_with_seismic_model(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho):


	rho_err_tot = chi_factor(mat_rho,seis_rho)
	vphi_err_tot = chi_factor(mat_vphi,seis_vphi)
	vs_err_tot = chi_factor(mat_vs,seis_vs)
    	err_tot=rho_err_tot+vphi_err_tot+vs_err_tot

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
