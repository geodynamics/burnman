import numpy
import bisect
import matplotlib.pyplot as pyplot
from tools import *

# this loads the PREM seismic velocities from prem_table.txt

# pressure: in GPa
# return: V_p, V_s in km/s
def prem_V(pressure):
    return lookup_(pressure, 3), lookup_(pressure, 4)

def prem_radius(pressure):
    return lookup_(pressure, 0)

def prem_density(pressure):
    return lookup_(pressure, 2)

#return linear interpolated column colidx
def lookup_(pressure, colidx):
    idx = bisect.bisect_left(table_p, pressure) - 1

    if (idx < 0):
        return table[0][colidx]
    elif (idx < len(table_p)-1):
        return linear_interpol(pressure, table_p[idx], table_p[idx+1], table[idx][colidx], table[idx+1][colidx])
    else:
        return table[idx][colidx]

#radius pressure density V_p V_s
table=[] 

# TODO: use tools.read_table
table = read_table("data/prem_table.txt")
table = sort_table(table, 1)
table = numpy.array(table)
table[:,1] = table[:,1] * 0.1 # convert kbar to GPa
table_p = table[:,1]


# test
if __name__ == "__main__":
    p = numpy.arange(1.0,360.0,3)
    vp = [prem_V(y)[0] for y in p]
    vs = [prem_V(y)[1] for y in p]
    pyplot.plot(p,vp,'+-r')
    pyplot.plot(p,vs,'+-')
    pyplot.show()
