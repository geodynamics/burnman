# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt
import birch_murnaghan as bm
import debye


#calculate isotropic thermal pressure, see
# Matas et. al. (2007) eq B4
def thermal_pressure(T,V, params):
    Debye_T = debye_temperature(params['ref_V']/V, params) 
    gr = grueneisen_parameter(params['ref_V']/V, params)
    P_th = gr * debye.thermal_energy(T,Debye_T, params['n'])/V
    return P_th

#compute the Debye temperature in K.  Takes the
#parameter x, which is ref_V/V (molar volumes).
#Depends on the reference grueneisen parameter,
#the reference Debye temperature, and the factor
#q0, see Matas eq B6
def debye_temperature(x, params):
    return params['ref_Debye']*np.exp((params['ref_grueneisen']- \
        grueneisen_parameter(x, params))/params['q0'])

#compute the grueneisen parameter with depth, according
#to q0.  Takes x=ref_V/V. See Matas eq B6
def grueneisen_parameter(x, params):
    return params['ref_grueneisen']*pow(1./x, params['q0'])

#invert the mie-grueneisen-debye eqn of state with thermal corrections
#for the volume
def volume(pressure,T,params):
    func = lambda x: bm.birch_murnaghan(params['ref_V']/x, params) + \
            thermal_pressure(T, x, params) - \
            thermal_pressure(300, x, params) - pressure
    V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V'])
    return V

#do the forward problem of the mie-grueneisen-debye equation of state, i.e.
#get pressure from volume
def pressure(T, V, params):
    return bm.birch_murnaghan(params['ref_V']/V, params) + \
            thermal_pressure(T,V, params) - \
            thermal_pressure(300,V, params)
    
#calculate the thermal correction for the mgd
#bulk modulus (see matas et al, 2007)
def thermal_bulk_modulus(T,V,params):
    gr = grueneisen_parameter(params['ref_V']/V, params)
    Debye_T = debye_temperature(params['ref_V']/V, params) 
    K_th = 3.*params['n']*debye.R*T/V * gr * \
        ((1. - params['q0'] - 3.*gr)*debye.debye_fn(Debye_T/T)+3.*gr*(Debye_T/T)/(np.exp(Debye_T/T) - 1.)) # EQ B5
    return K_th

#calculate the mgd bulk modulus (K_T) as a function of P, T, and V
def bulk_modulus(T,V, params):
    K_T = bm.bulk_modulus(V, params) + \
        thermal_bulk_modulus(T,V, params) - \
        thermal_bulk_modulus(300.,V, params)  #EQB13
    return K_T


#calculate the thermal correction to the shear modulus as a function of V, T
def thermal_shear_modulus(T, V, params):
    gr = grueneisen_parameter(params['ref_V']/V, params)
    Debye_T = debye_temperature(params['ref_V']/V, params) 
    mu_th= 3./5. * ( thermal_bulk_modulus(T,V,params) - \
            6*debye.R*T*params['n']/V * gr * debye.debye_fn(Debye_T/T) ) # EQ B10
    return mu_th

#calculate the mgd shear modulus as a function of P, V, and T
def shear_modulus(T,V, params):
    mu = bm.shear_modulus(V,params) + \
        thermal_shear_modulus(T,V, params) - \
        thermal_shear_modulus(300.,V, params) # EQ B11
    return mu

#heat capacity at constant volume
def heat_capacity_v(T, V, params):
    Debye_T = debye_temperature(params['ref_V']/V, params)
    C_v = debye.heat_capacity_v(T, Debye_T, params['n'])
    return C_v

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

