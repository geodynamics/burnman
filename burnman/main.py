# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import os, sys, numpy as np
import matplotlib.pyplot as plt

import partitioning as part
import mie_grueneisen_debye
import minerals 
import seismic
import tools
import averaging_schemes
import scipy.integrate as integrate
from composite import composite

#phase = namedtuple('phase', ['mineral', 'fraction'])

class elastic_properties:
    """
    class that contains volume V, density rho, bulk_modulus K and 
    shear_modulus G for a list of pressures.
    """
    def __init__(self, V=None, rho=None, K=None, G=None, fraction=None):
        self.V = V
        self.rho = rho
        self.K = K
        self.G = G
        self.fraction = fraction

    def set_size(self, size):
        self.V = np.ndarray(size)
        self.rho = np.ndarray(size)
        self.K = np.ndarray(size)
        self.G = np.ndarray(size)
        self.fraction = np.ndarray(size)


def calculate_moduli(rock, pressures, temperatures):
    moduli = [ elastic_properties() for ph in rock.phases ]
    for m in moduli:
        m.set_size(len(pressures))
        
    for idx in range(len(pressures)):
        rock.set_state(pressures[idx], temperatures[idx])
        
        for midx in range(len(moduli)):
            m = moduli[midx]
            ph = rock.phases[midx]
            m.V[idx] = ph.fraction * ph.mineral.molar_volume()
            m.K[idx] = ph.mineral.adiabatic_bulk_modulus()
            m.G[idx] = ph.mineral.shear_modulus()
            m.rho[idx] = ph.mineral.molar_mass() / ph.mineral.molar_volume()
            m.fraction[idx] = ph.fraction
        
    return moduli

def average_moduli(moduli_list, averaging_scheme=averaging_schemes.voigt_reuss_hill):
    n_pressures = len(moduli_list[0].V)
    result = elastic_properties()
    result.set_size(n_pressures)
    
    for idx in range(n_pressures):
        fractions = [m.fraction[idx] for m in moduli_list]
        V_ph = [m.V[idx] for m in moduli_list]
        K_ph = [m.K[idx] for m in moduli_list]
        G_ph = [m.G[idx] for m in moduli_list]
        massfraction_ph = [m.rho[idx]*m.V[idx] for m in moduli_list]
               
        result.V[idx] = sum(V_ph)
        
        result.K[idx] = averaging_scheme.average(fractions, V_ph, K_ph)
        result.G[idx] = averaging_scheme.average(fractions, V_ph, G_ph)
        result.rho[idx] = sum(massfraction_ph) / result.V[idx]
        result.fraction[idx] = 1.0
    return result

def compute_velocities(moduli):
    mat_vs = np.ndarray(len(moduli.V))
    mat_vp = np.ndarray(len(moduli.V))
    mat_vphi = np.ndarray(len(moduli.V))
    
    for i in range(len(moduli.V)):
        mat_vs[i] = np.sqrt( moduli.G[i] / moduli.rho[i])
        mat_vp[i] = np.sqrt( (moduli.K[i] + 4./3.*moduli.G[i]) / moduli.rho[i])
        mat_vphi[i] = np.sqrt( moduli.K[i] / moduli.rho[i])
    
    return mat_vp, mat_vs, mat_vphi
 
 
def velocities_from_rock(rock, pressures, temperatures, averaging_scheme=averaging_schemes.voigt_reuss_hill()):
    
    moduli_list = calculate_moduli(rock, pressures, temperatures)
    moduli = average_moduli(moduli_list, averaging_scheme)
    mat_vp, mat_vs, mat_vphi = compute_velocities(moduli)
    return moduli.rho, mat_vp, mat_vs, mat_vphi, moduli.K, moduli.G

# def whole_darn_thing_one_mat(mineral, pressures, temperatures):
# 
#     rock = composite( ( (mineral, 1.0) ) )
#     moduli_list = calculate_moduli(rock, pressures, temperatures)
#     mat_vp, mat_vs, mat_vphi = compute_velocities(moduli_list[0])
#     return moduli_list[0].rho, mat_vp, mat_vs, mat_vphi, moduli_list[0].K, moduli_list[0].G 






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
