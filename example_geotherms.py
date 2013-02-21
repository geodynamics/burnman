# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# under GPL v2 or later.

"""
Shows the various ways to input geotherms: Built-in geotherms (geotherm1 and 2), basic linear (geotherm3),
loaded in from a data file (geotherm4) of your choice. Geotherm 1 is from Brown & Shankland (1981) and 
geotherm2 from Watson & Baxter (2007).

requires:

teaches:
- geotherms

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

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
#table = tools.read_table("data/example_geotherm.txt")

table_pressure = np.array(table)[:,0]
table_temperature = np.array(table)[:,1]

my_geotherm = lambda p:  burnman.tools.lookup_and_interpolate(table_pressure, table_temperature, p)
temperature4 = [my_geotherm(p) for p in pressures]

#you can also look at burnman/geotherm.py to see how the geotherms are implemented


plt.plot(pressures/1e9,temperature1,'-r',label="Brown, Shankland")
plt.plot(pressures/1e9,temperature2,'-g',label="Watson, Baxter")
plt.plot(pressures/1e9,temperature3,'-b',label="handwritten linear")
plt.plot(pressures/1e9,temperature4,'-k',label="handwritten from table")

plt.legend(loc='lower right')
plt.xlim([0, 130])
plt.xlabel('Pressure/GPa')
plt.ylabel('Temperature')
plt.show()

