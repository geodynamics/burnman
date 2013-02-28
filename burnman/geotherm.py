# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import bisect
import matplotlib.pyplot as pyplot
from tools import *
import seismic
import scipy.integrate as integrate
from minerals import *

import mie_grueneisen_debye as mgd

# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in GPa
# return: temperature in K
def watson_baxter(pressure):
    if (pressure <= 15e9):
        return 1900-1420*pow(0.8,pressure/1e9)
    else:
        return 1680+11.1*pressure/1e9





# geotherm from Brown and Shankland 81
def brown_shankland(pressure):
    depth = seismic.prem_model.depth(pressure)
    return lookup_and_interpolate(table_brown_depth, table_brown_temperature, depth)    

#This integrates dT/dP = gr * T / K_s
def self_consistent(pressures, T0, phases, molar_abundances):
    temperatures = integrate.odeint(dTdP, T0, pressures, args=(phases,molar_abundances))
    return temperatures

#ODE to integrate temperature with depth for a composite material
#Assumes that the minerals exist at a common pressure (Reuss bound, should be good for
#slow deformations at high temperature), as well as an adiabatic process.  This
#corresponds to conservation of enthalpy.  Is this formula correct?  It works for
#single phases, and seems right for composites -- the tricky thing is working out
#the heat exchange between phases as they heat up at different rates due to differing
#grueneisen parameters and bulk moduli
def dTdP(temperature, pressure, phases, molar_abundances):
    top = 0
    bottom = 0
    for i in range(len(phases)):
        phases[i].set_state(pressure, temperature)

        gr = phases[i].grueneisen_parameter()
        K_s = phases[i].adiabatic_bulk_modulus()*1.e9
        C_p = phases[i].heat_capacity_p()

        top += molar_abundances[i]*gr*C_p/K_s
        bottom += molar_abundances[i]*C_p
    
    return temperature*top/bottom


table_brown = read_table("data/brown_81.txt")
table_brown_depth = np.array(table_brown)[:,0]
table_brown_temperature = np.array(table_brown)[:,1]

# test geotherm
if __name__ == "__main__":
    p = np.arange(1.0e9,128.0e9,3e9)
    t1 = [watson_baxter(y) for y in p]
    t2 = [brown_shankland(y) for y in p]
  
    test_mineral = material()
    test_mineral.params ={'name':'test',
                          'ref_V': 6.844e-6,
                          'ref_K': 256,
                          'K_prime':0.0,
                          'ref_mu': 175.,
                          'mu_prime': 1.7,
                          'molar_mass': 0.0,
                          'n': 1.,
                          'ref_Debye': 1.,
                          'ref_grueneisen': 1.5,
                          'q0': .0001}
 
    pressure = np.linspace(0., 140.e9, 100)
    fp = ferropericlase(0.4)
    pv = mg_fe_perovskite(0.7)
    phases = [pv,fp]
    for ph in phases:
        ph.set_method('mgd')
    molar_abundances = [7., .3]
    t3 = self_consistent(p, 1600, phases, molar_abundances)

    p1,=pyplot.plot(p,t1,'x--r')
    p2,=pyplot.plot(p,t2,'*-g')
    p3,=pyplot.plot(p,t3,'*-b')
    #pyplot.xlim(25,135)
    #pyplot.ylim(1600e9,3100e9)
    pyplot.legend([p1,p2,p3],[ "watson", "brown", "self consistent"], loc=4)

    pyplot.show()
