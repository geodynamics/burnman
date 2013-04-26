# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows each of the geotherms currently possible with BurnMan.
These are:

1. Brown and Shankland, 1981
2. Watson and Baxter, 2007
3. linear extrapolation
4. Read in from file from user
5. Adiabatic from potential temperature and choice of mineral (pyrolite in this example)

requires:

- input geotherm files
- mineral definitions for adiabat
- burnman/composite.py for adiabat (if necessary)

teaches:
- geotherms

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
    geotherm1 = burnman.geotherm.brown_shankland
    temperature1 = [geotherm1(p) for p in pressures]
    geotherm2 = burnman.geotherm.watson_baxter
    temperature2 = [geotherm2(p) for p in pressures]
    
    #a geotherm is actually just a function that returns a temperature given pressure in Pa
    #so we can just write our own function
    geotherm3 = lambda p:  1500+(2500-1500)*p/128e9
    temperature3 = [geotherm3(p) for p in pressures]
    
    #what about a geotherm defined from datapoints given in a file (our inline)?
    table = [[1e9,1600],[30e9,1700],[130e9,2700]]
    #this could also be loaded from a file, just uncomment this
    #table = tools.read_table("input_geotherm/example_geotherm.txt")

    table_pressure = np.array(table)[:,0]
    table_temperature = np.array(table)[:,1]
    
    my_geotherm = lambda p:  burnman.tools.lookup_and_interpolate\
    (table_pressure, table_temperature, p)
    temperature4 = [my_geotherm(p) for p in pressures]


    #finally, we can also calculate a self consistent 
    #geotherm for an assemblage of minerals
    #based on self compression of the composite rock.  
    #First we need to define an assemblage
    pyrolite = burnman.composite( [ (minerals.mg_fe_perovskite(0.1), 0.7), 
                                   (minerals.ferropericlase(0.4),   0.3) ] )
    pyrolite.set_method("mgd3")
    #next, define an anchor temperature at which we are starting.  
    #Perhaps 1500 K for the upper mantle
    T0 = 1500.
    #then generate temperature values using the self consistent function. 
    # This takes more time than the above methods
    temperature5 = burnman.geotherm.self_consistent(pressures, T0, pyrolite)
    
    #you can also look at burnman/geotherm.py to see how the geotherms are implemented
    
    
    plt.plot(pressures/1e9,temperature1,'-r',label="Brown, Shankland")
    plt.plot(pressures/1e9,temperature2,'-g',label="Watson, Baxter")
    plt.plot(pressures/1e9,temperature3,'-b',label="handwritten linear")
    plt.plot(pressures/1e9,temperature4,'-k',label="handwritten from table")
    plt.plot(pressures/1e9,temperature5,'-m',label="Adiabat with pv (70%) and fp(30%)")
    
    plt.legend(loc='lower right')
    plt.xlim([0, 130])
    plt.xlabel('Pressure/GPa')
    plt.ylabel('Temperature')
    plt.show()
    
