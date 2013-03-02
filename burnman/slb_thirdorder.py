# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

#system libs:
import numpy
import scipy.optimize as opt
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt

#own libs:
from tools import *
import birch_murnaghan as bm
import debye

# TODO: add up weight percent and check <100 and tell them how much
## Based on Stixrude & Lithgow-Bertelloni (2005), all equation numbers refer to this paper. 


#compute the grueneisen parameter with depth, according
#to q0.  Takes x=ref_V/V.
def grueneisen_parameter(x, params):
    f = 1./2. * (pow(x, 2./3.) - 1.)
    a1_ii = 6. * params['ref_grueneisen'] # EQ 47
    a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
    nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # EQ 41
    gr = 1./6./nu_o_nu0_sq * (2*f+1.) * ( a1_ii + a2_iikk*f )
    return gr

def isotropic_eta_s(x, params):
    f = 1./2. * (pow(x, 2./3.) - 1.)
    a2_s = -2.*params['ref_grueneisen'] - 2.*params['eta_0s'] # EQ 47 
    a1_ii = 6. * params['ref_grueneisen'] # EQ 47
    a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
    nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # EQ 41
    gr = 1./6./nu_o_nu0_sq * (2*f+1.) * ( a1_ii + a2_iikk*f )
    eta_s = - gr - (1./2. * pow(nu_o_nu0_sq,-1.) * pow((2.*f)+1.,2.)*a2_s) # EQ 46 NOTE the typo from Stixrude 2005
    return eta_s

def debye_temperature(x,params):
    return params['ref_Debye']


def volume(p,T,params):
    debye_T = lambda x : debye_temperature(params['ref_V']/x, params)
    gr = lambda x : grueneisen_parameter(params['ref_V']/x, params)
    E_th =  lambda x : debye.thermal_energy(T,debye_T(x), params['n']) #thermal energy at temperature T
    E_th_ref = lambda x : debye.thermal_energy(300.,debye_T(x), params['n']) #thermal energy at reference temperature
      
    b_iikk= 9.*params['ref_K'] # EQ 28
    b_iikkmm= 27.*params['ref_K']*(params['K_prime']-4.) # EQ 29

    f = lambda x: 0.5*(pow(params['ref_V']/x,2./3.)-1.) # EQ 24
    func = lambda x: (1./3.)*(pow(1.+2.*f(x),5./2.))*((b_iikk*f(x)) \
            +(0.5*b_iikkmm*pow(f(x),2.))) + gr(x)*(E_th(x) - E_th_ref(x))/x - p #EQ 21 
    
    V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V']) 
    return V

def shear_modulus(T,V,params):
    debye_T = debye_temperature(params['ref_V']/V, params)
    gr = grueneisen_parameter(params['ref_V']/V, params)
    eta_s = isotropic_eta_s(params['ref_V']/V, params)

    E_th = debye.thermal_energy(T,debye_T, params['n'])
    E_th_ref = debye.thermal_energy(300.,debye_T, params['n'])

    f =.5*(pow(params['ref_V']/V,2./3.)-1.) # EQ 24
  
    G_twoterms = bm.shear_modulus_second_order(V, params) - eta_s * (E_th-E_th_ref) / V
    G_threeterms = bm.shear_modulus_third_order(V, params) - eta_s * (E_th-E_th_ref) / V

    return G_threeterms
    


def bulk_modulus(T,V,params):

    debye_T = debye_temperature(params['ref_V']/V, params)
    gr = grueneisen_parameter(params['ref_V']/V, params)

    E_th = debye.thermal_energy(T,debye_T, params['n']) #thermal energy at temperature T
    E_th_ref = debye.thermal_energy(300.,debye_T, params['n']) #thermal energy at reference temperature
    
    C_v = debye.heat_capacity_v(T,debye_T, params['n']) #heat capacity at temperature T
    C_v_ref = debye.heat_capacity_v(300.,debye_T, params['n']) #heat capacity at reference temperature

    f =.5*(pow(params['ref_V']/V,2./3.)-1) # EQ 24
    
    K = bm.bulk_modulus(V, params) \
             + (gr + 1.-params['q0'])* ( gr / V ) * (E_th - E_th_ref) \
             - ( pow(gr , 2.) / V )*(C_v*T - C_v_ref*300.)   
    
    return K

#heat capacity at constant volume
def heat_capacity_v(T, V, params):
    debye_T = debye_temperature(params['ref_V']/V, params)
    return debye.heat_capacity_v(T, debye_T,params['n'])

def thermal_expansivity(T,V,params):
    C_v = heat_capacity_v(T,V,params)
    gr = grueneisen_parameter(params['ref_V']/V, params)
    K = bulk_modulus(T,V,params)
    alpha = gr * C_v / K / V
    return alpha

#heat capacity at constant pressure
def heat_capacity_p(T,V,params):
    alpha = thermal_expansivity(T,V,params)
    gr = grueneisen_parameter(params['ref_V']/V, params)
    C_v = heat_capacity_v(T,V,params)
    C_p = C_v*(1. + gr * alpha * T)
    return C_p

#calculate the adiabatic bulk modulus (K_S) as a function of P, T, and V
# alpha is basically 1e-5
def bulk_modulus_adiabatic(T,V,params):
    K_T=bulk_modulus(T,V,params)
    alpha = thermal_expansivity(T,V,params)
    gr = grueneisen_parameter(params['ref_V']/V, params)
    K_S = K_T*(1. + gr * alpha * T)
    return K_S


