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


import numpy
import bisect
import matplotlib.pyplot as pyplot
from tools import *
import seismic

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
def self_consistent(pressure, T0, params):
        minP = 0
        #upper mantle potential temperature, in K
	# integrate the adiabatic gradient equation
	lnT = integrate.quad( lambda x: (params['ref_grueneisen']/(1.e9*bm_bulk_modulus(x, params))), minP, pressure)
	T = T0*np.exp(lnT[0])
	return T
    
    


table_brown = read_table("data/brown_81.txt")
table_brown_depth = numpy.array(table_brown)[:,0]
table_brown_temperature = numpy.array(table_brown)[:,1]

# test geotherm
if __name__ == "__main__":
    p = numpy.arange(1.0,128.0,3)
    t1 = [geotherm_watsona_baxter(y) for y in p]
    t2 = [geotherm_brown_shankland(y) for y in p]
    p2,=pyplot.plot(p,t2,'x--r')
    p3,=pyplot.plot(p,t3,'*-g')
    pyplot.xlim(25,135)
    pyplot.ylim(1600,3100)
    pyplot.legend([p1,p2],[ "watson", "brown"], loc=4)

    pyplot.show()
