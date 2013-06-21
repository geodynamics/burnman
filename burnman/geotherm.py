# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import bisect
import matplotlib.pyplot as pyplot
from tools import *
import seismic
import scipy.integrate as integrate
import minerals
from composite import *


# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in Pa
# return: temperature in K
def watson_baxter(pressure):
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
      if (pressure[i] <= 15e9):
          temperature[i] = 1900-1420*pow(0.8,pressure[i]/1e9)
      else:
          temperature[i] = 1680+11.1*pressure[i]/1e9
    return temperature





# geotherm from Brown and Shankland 1981
def brown_shankland(pressure):
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
      depth = seismic.prem_model.depth(pressure[i])/1.e3
      temperature[i] = lookup_and_interpolate(table_brown_depth, table_brown_temperature, depth)    
    return temperature

# geotherm from Anderson 1982
def anderson(pressure):
    temperature = np.empty_like(pressure)
    for i in range(len(pressure)):
      depth = seismic.prem_model.depth(pressure[i])/1.e3
      temperature[i] = lookup_and_interpolate(table_anderson_depth, table_anderson_temperature, depth)    
    return temperature

#This integrates dT/dP = gr * T / K_s
def adiabatic(pressures, T0, rock):
    temperatures = integrate.odeint(lambda t,p : dTdP(t,p,rock), T0, pressures)
    return temperatures.ravel()

#ODE to integrate temperature with depth for a composite material
#Assumes that the minerals exist at a common pressure (Reuss bound, should be good for
#slow deformations at high temperature), as well as an adiabatic process.  This
#corresponds to conservation of enthalpy.  Is this formula correct?  It works for
#single phases, and seems right for composites -- the tricky thing is working out
#the heat exchange between phases as they heat up at different rates due to differing
#grueneisen parameters and bulk moduli
def dTdP(temperature, pressure, rock):
    top = 0
    bottom = 0
    for ph in rock.phases:
        ph.mineral.set_state(pressure, temperature)

        gr = ph.mineral.grueneisen_parameter()
        K_s = ph.mineral.adiabatic_bulk_modulus()
        C_p = ph.mineral.heat_capacity_p()

        top += ph.fraction*gr*C_p/K_s
        bottom += ph.fraction*C_p
    
    return temperature*top/bottom


table_brown = read_table("input_geotherm/brown_81.txt")
table_brown_depth = np.array(table_brown)[:,0]
table_brown_temperature = np.array(table_brown)[:,1]

table_anderson = read_table("input_geotherm/anderson_82.txt")
table_anderson_depth = np.array(table_anderson)[:,0]
table_anderson_temperature = np.array(table_anderson)[:,1]

# test geotherm
if __name__ == "__main__":
    p = np.arange(1.0e9,128.0e9,3e9)
  
    pyrolite = composite( [ (minerals.SLB2011.mg_fe_perovskite(0.2), 0.8), (minerals.SLB2011.ferropericlase(0.4), 0.2) ] )
    pyrolite.set_method('slb3')
    pyrolite.set_state(40.e9, 2000)

    t1 = watson_baxter(p)
    t2 = brown_shankland(p)
    t3 = adiabatic(p, 1600, pyrolite)

    p1,=pyplot.plot(p,t1,'x--r')
    p2,=pyplot.plot(p,t2,'*-g')
    p3,=pyplot.plot(p,t3,'*-b')
    pyplot.legend([p1,p2,p3],[ "watson", "brown", "self consistent"], loc=4)

    pyplot.show()
