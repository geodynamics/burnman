import numpy
import bisect
import matplotlib.pyplot as pyplot
import prem
from tools import *

# loads a simple geotherm from geotherm.txt


# geotherm from geotherm.txt, from Cayman
# pressure: in GPa
# return: temperature in K
def geotherm(pressure):
    idx = bisect.bisect_left(table_p, pressure) - 1
    if (idx < 0):
        return table_T[0]
    elif (idx < len(table_p)-1):
        return linear_interpol(pressure, table_p[idx], table_p[idx+1], table_T[idx], table_T[idx+1])
    else:
        return table_T[idx]


# polynomial fit from Watson, Baxter, EPSL, 2007
# pressure: in GPa
# return: temperature in K
def geotherm_formula(pressure):
    if (pressure <= 15):
        return 1900-1420*pow(0.8,pressure)
    else:
        return 1680+11.1*pressure





# geotherm from Brown81
def geotherm_brown(pressure):
    depth = 6371. - prem.prem_radius(pressure)
    idx = bisect.bisect_left(table_brown_depth, depth) - 1
    if (idx < 0):
        return table_brown[0][1]
    elif (idx < len(table_brown)-1):
        return linear_interpol(depth, table_brown_depth[idx], table_brown_depth[idx+1], table_brown[idx][1], table_brown[idx+1][1])
    else:
        return table_brown[idx][1]
    
    


table_brown = read_table("data/brown_81.txt")
table_brown_depth = numpy.array(table_brown)[:,0]

    

geotherm_table = read_table("data/geotherm.txt")

geotherm_table = sort_table(geotherm_table, 0)
table_p=numpy.array(geotherm_table)[:,0]
table_T=numpy.array(geotherm_table)[:,1]







# test geotherm
if __name__ == "__main__":
    p = numpy.arange(1.0,128.0,3)
    t = [geotherm(y) for y in p]
    t2 = [geotherm_formula(y) for y in p]
    t3 = [geotherm_brown(y) for y in p]
    p1,=pyplot.plot(p,t,'+-')
    p2,=pyplot.plot(p,t2,'x--r')
    p3,=pyplot.plot(p,t3,'*-g')
    pyplot.xlim(25,135)
    pyplot.ylim(1600,3100)
    pyplot.legend([p1,p2,p3],["cayman", "watson", "brown"], loc=4)

    pyplot.show()
