# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import os, sys, numpy as np
import matplotlib.pyplot as plt

import composition as part
import mie_grueneisen_debye
import minerals 
import seismic
import tools
import averaging_schemes as ave 
import scipy.integrate as integrate


def calculate_velocities (pressure, temperature, rock):

    mat_vs = np.empty_like(pressure)
    mat_vp = np.empty_like(pressure)
    mat_vphi = np.empty_like(pressure)
    mat_rho = np.empty_like(pressure)
    mat_K = np.empty_like(pressure)
    mat_mu = np.empty_like(pressure)
    #print "Calculating elastic properties for phase assemblage \n"
    #print "seismic p (GPa)    T (K)    density(kg/m^3)    K(Gpa) G(GPa)    Vs (km/s)    Vp(km/s)    Vphi (km/s)"
    for i in range(len(pressure)):
        #determine the partial molar volumes of the phases
        rock.set_state(pressure[i], temperature[i])
        V_ph = [(ph.fraction*ph.mineral.molar_volume()) for ph in rock.phases]
        V_tot= sum(V_ph)
        #calculate the density of the phase assemblage
        mat_rho[i] = (1./V_tot) * sum( ph.fraction* ph.mineral.molar_mass() for ph in rock.phases)
        K_ph = [ph.mineral.adiabatic_bulk_modulus() for ph in rock.phases]
        mat_K[i] = ave.voigt_reuss_hill(V_ph,K_ph)
        mu_ph = [ph.mineral.shear_modulus() for ph in rock.phases]
        mat_mu[i] = ave.voigt_reuss_hill(V_ph,mu_ph)
        #compute seismic velocities
        mat_vs[i] = np.sqrt( mat_mu[i] / mat_rho[i])
        mat_vp[i] = np.sqrt( (mat_K[i] + 4./3.*mat_mu[i]) / mat_rho[i])
        mat_vphi[i] = np.sqrt( (mat_K[i]) / mat_rho[i])

    return mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu

# apply attenuation correction to a list of v_p, v_s, v_phi
# takes a list of pressures and the seismic model
def apply_attenuation_correction(v_p,v_s,v_phi,Qs,Qphi):
    length = len(v_p)
    ret_v_p = np.zeros(length)
    ret_v_s = np.zeros(length)
    ret_v_phi = np.zeros(length)
    for i in range(length):
        ret_v_p[i],ret_v_s[i],ret_v_phi[i] = \
            seismic.attenuation_correction(v_p[i], v_s[i], v_phi[i],Qs,Qphi)
    
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
