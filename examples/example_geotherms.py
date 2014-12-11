# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
example_geotherms
-----------------

This example shows each of the geotherms currently possible with BurnMan.
These are:

1. Brown and Shankland, 1981
2. Anderson, 1982
3. Watson and Baxter, 2007
4. linear extrapolation
5. Read in from file from user
6. Adiabatic from potential temperature and choice of mineral (pyrolite in this example)

*Uses:*

* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.geotherm.anderson`
* input geotherm file *input_geotherm/example_geotherm.txt* (optional)
* :class:`burnman.composite.Composite` for adiabat

*Demonstrates:*

* the available geotherms

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":

# we want to evaluate several geotherms at these values
    pressures = np.arange(1e9,128e9,3e9)

    #load two builtin geotherms and evaluate the temperatures at all pressures
    temperature1 = burnman.geotherm.brown_shankland(pressures)
    temperature2 = burnman.geotherm.anderson(pressures)

    #a geotherm is actually just a function that returns a list of temperatures given pressures in Pa
    #so we can just write our own function
    my_geotherm_function = lambda p:  [ 1500+(2500-1500)*x/128e9 for x in p]
    temperature3 = my_geotherm_function(pressures)

    #what about a geotherm defined from datapoints given in a file (our inline)?
    table = [[1e9,1600],[30e9,1700],[130e9,2700]]
    #this could also be loaded from a file, just uncomment this
    #table = burnman.tools.read_table("input_geotherm/example_geotherm.txt")

    table_pressure = np.array(table)[:,0]
    table_temperature = np.array(table)[:,1]

    my_geotherm_interpolate = lambda p:  [ np.interp(x, table_pressure, \
                                                     table_temperature) for x in p]

    temperature4 = my_geotherm_interpolate(pressures)


    #finally, we can also calculate a self consistent
    #geotherm for an assemblage of minerals
    #based on self compression of the composite rock.
    #First we need to define an assemblage
    pyrolite = burnman.Composite([0.7, 0.3],
                                 [minerals.SLB_2005.mg_fe_perovskite(0.1),
                                  minerals.SLB_2005.ferropericlase(0.4)])
    #next, define an anchor temperature at which we are starting.
    #Perhaps 1500 K for the upper mantle
    T0 = 1500.
    #then generate temperature values using the self consistent function.
    # This takes more time than the above methods
    temperature5 = burnman.geotherm.adiabatic(pressures, T0, pyrolite)

    #you can also look at burnman/geotherm.py to see how the geotherms are implemented


    plt.plot(pressures/1e9,temperature1,'-r',label="Brown, Shankland")
    plt.plot(pressures/1e9,temperature2,'-c',label="Anderson")
    plt.plot(pressures/1e9,temperature3,'-b',label="handwritten linear")
    plt.plot(pressures/1e9,temperature4,'-k',label="handwritten from table")
    plt.plot(pressures/1e9,temperature5,'-m',label="Adiabat with pv (70%) and fp(30%)")

    plt.legend(loc='lower right')
    plt.xlim([0, 130])
    plt.xlabel('Pressure/GPa')
    plt.ylabel('Temperature')
    plt.savefig("output_figures/example_geotherm.png")
    plt.show()

