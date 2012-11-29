import numpy
import bisect
import matplotlib.pyplot as pyplot
import seismic_models as seis
from tools import *


# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in GPa
# return: temperature in K
def geotherm_watson_baxter(pressure):
    if (pressure <= 15):
        return 1900-1420*pow(0.8,pressure)
    else:
        return 1680+11.1*pressure





# geotherm from Brown and Shankland 81
def geotherm_brown_shankland(pressure):
    depth = seis.prem_depth(pressure/1e9)/1.e3
    idx = bisect.bisect_left(table_brown_depth, depth) - 1
    if (idx < 0):
        return table_brown[0][1]
    elif (idx < len(table_brown)-1):
        return linear_interpol(depth, table_brown_depth[idx], table_brown_depth[idx+1], table_brown[idx][1], table_brown[idx+1][1])
    else:
        return table_brown[idx][1]


#This integrates dT/dP = gr * T / K_s
def geotherm_self_consistent(pressure, T0, params):
        minP = 0
        #upper mantle potential temperature, in K
	# integrate the adiabatic gradient equation
	lnT = integrate.quad( lambda x: (params['ref_grueneisen']/(1.e9*bm_bulk_modulus(x, params))), minP, pressure)
	T = T0*np.exp(lnT[0])
	return T
    
    


table_brown = read_table("data/brown_81.txt")
table_brown_depth = numpy.array(table_brown)[:,0]

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
