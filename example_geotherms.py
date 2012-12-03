import os, sys, numpy as np
import matplotlib.pyplot as plt

lib_path = os.path.abspath('code/')
sys.path.append(lib_path)

from code import minerals 
from code import main as main
from code import tools

# we want to evaluate several geotherms at these values
pressures = np.arange(1e9,128e9,3e9)

#load two builtin geotherms and evaluate the temperatures at all pressures
geotherm1 = main.get_geotherm('brown_shankland')
temperature1 = [geotherm1(p) for p in pressures]
geotherm2 = main.get_geotherm('watson_baxter')
temperature2 = [geotherm2(p) for p in pressures]

#a geotherm is actually just a function that returns a pressure given pressure in Pa
#so we can just write our own function
geotherm3 = lambda p:  1500+(2500-1500)*p/128e9
temperature3 = [geotherm3(p) for p in pressures]

#what about a geotherm defined from datapoints given in a file (our inline)?
table = [[1e9,1600],[30e9,1700],[130e9,2700]]
#this could also be loaded from a file, just uncomment this
#table = tools.read_table("data/example_geotherm.txt")

table_pressure = np.array(table)[:,0]
table_temperature = np.array(table)[:,1]

my_geotherm = lambda p:  tools.lookup_and_interpolate(table_pressure, table_temperature, p)
temperature4 = [my_geotherm(p) for p in pressures]

#you can also look at code/geotherm.py to see how the geotherms are implemented


plt.plot(pressures,temperature1,'-r',label="Brown Shankland")
plt.plot(pressures,temperature2,'-g',label="Watson Baxter")
plt.plot(pressures,temperature3,'-b',label="handwritten linear")
plt.plot(pressures,temperature4,'x-k',label="handwritten from table")

plt.legend(loc=4)
plt.xlabel('Pressure/Pa')
plt.ylabel('Temperature')
plt.show()

