# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt

"""
Functions for the Debye model.  Note that this is not Mie-Grueneisen-Debye, 
just Debye, so is pretty limited.  Combine this with Mie-Grueneisen and
Birch-Murnaghan to get a full EOS
"""

R = 8.314462175

#Evaluate the Debye function.  Takes the parameter
#xi = Debye_T/T
def debye_fn(x):
    sol = integrate.quad( lambda xi: pow(xi,3.)/(np.exp(xi)-1.) , 0.0, x) # EQ B3
    return 3.*sol[0]/pow(x,3.)

#calculate the thermal energy of a substance.  Takes the temperature,
#the Debye temperature, and n, the number of atoms per molecule
def thermal_energy(T, debye_T, n):
    if T == 0:
      return 0
    E_th = 3.*n*R*T * debye_fn(debye_T/T)
    return E_th

#heat capacity at constant volume
def heat_capacity_v(T, debye_T, n):
    if T ==0:
      return 0
    deb = integrate.quad( lambda x : pow(x,4.)*np.exp(x)/pow((np.exp(x)-1.),2.), 0.0, debye_T/T)
    C_v = 9.*n*R*deb[0]/pow(debye_T/T,3.)
    return C_v

