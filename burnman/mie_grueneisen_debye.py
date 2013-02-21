# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt
import birch_murnaghan as bm


#Evaluate the Debye function.  Takes the parameter
#xi = Debye_T/T
def debye_fn(x):
	sol = integrate.quad( lambda xi: pow(xi,3.)/(np.exp(xi)-1.) , 0.0, x) # EQ B3
	return 3.*sol[0]/pow(x,3.)

#calculate isotropic thermal pressure, see
# Matas et. al. (2007) eq B4
def mgd_thermal_pressure(V, T, params):
	R = 8.314462175
	Debye_T = debye_temperature(params['ref_V']/V, params) 
	gr = grueneisen_parameter(params['ref_V']/V, params)
	xi = Debye_T/T
	P_th = 3.*params['n']*gr*R*T/V * debye_fn(xi)
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
	func = lambda x: bm.birch_murnaghan(params['ref_V']/x, params)*1.e9 + \
			mgd_thermal_pressure(x, T, params) - \
			mgd_thermal_pressure(x, 300, params) - pressure
	V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V'])
	return V

#do the forward problem of the mie-grueneisen-debye equation of state, i.e.
#get pressure from volume
def pressure(volume, T, params):
	return bm.birch_murnaghan(params['ref_V']/volume, params)*1.e9 + \
			mgd_thermal_pressure(volume, T, params) - \
			mgd_thermal_pressure(volume, 300, params)
	
#calculate the thermal correction for the mgd
#bulk modulus (see matas et al, 2007)
def thermal_bulk_modulus(V, T, params):
	R = 8.314462175
	gr = grueneisen_parameter(params['ref_V']/V, params)
	Debye_T = debye_temperature(params['ref_V']/V, params) 
	K_th = 3.*params['n']*R*T/V * gr * \
		((1. - params['q0'] - 3.*gr)*debye_fn(Debye_T/T)+3.*gr*(Debye_T/T)/(np.exp(Debye_T/T) - 1.)) # EQ B5
	return K_th

#calculate the mgd bulk modulus (K_T) as a function of P, T, and V
def bulk_modulus(T,V, params):
	K_T = bm.bulk_modulus(V, params) + \
		thermal_bulk_modulus(V,T, params)/1.e9 - \
		thermal_bulk_modulus(V,300., params)/1.e9  #EQB13
	return K_T

#calculate the adiabatic bulk modulus (K_S) as a function of P, T, and V
# alpha is basically 1e-5
def bulk_modulus_adiabatic(T,V,params):
	K_T=bulk_modulus(T,V,params)
	alpha=1./K_T*((mgd_thermal_pressure(V,T+1.,params)-mgd_thermal_pressure(V,T-1.,params))/2.)/1.e9
	return K_T*(1.+alpha*grueneisen_parameter(params['ref_V']/V,params)*T)


#calculate the thermal correction to the shear modulus as a function of V, T
def thermal_shear_modulus(V, T,  params):
	R = 8.314462175
	gr = grueneisen_parameter(params['ref_V']/V, params)
	Debye_T = debye_temperature(params['ref_V']/V, params) 
        mu_th= 3./5. * ( thermal_bulk_modulus(V,T,params) - \
			6*R*T*params['n']/V * gr * debye_fn(Debye_T/T) ) # EQ B10
	return mu_th

#calculate the mgd shear modulus as a function of P, V, and T
def shear_modulus(T,V, params):
	mu = bm.shear_modulus(V,params) + \
		thermal_shear_modulus(V,T, params)/1.e9 - \
		thermal_shear_modulus(V,300, params)/1.e9 # EQ B11
	return mu


