# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals


garnet=minerals.HP_2011_ds62.garnet()
composition=np.array([ 0.5, 0.2, 0.1, 0.2 ])
garnet.set_composition(composition)

print garnet.molar_fraction
print garnet.alpha
print 'Ideal activities'
print garnet.ideal_activity
print ''

# Excess volumes for the pyrope-grossular join
n=100
pyrope_proportion= np.empty(shape=(n+1))
garnet_excess_volume= np.empty(shape=(n+1))
for i in range(n+1):
    pyrope_proportion[i]=float(i)/n
    composition=([ pyrope_proportion[i], 0.0, 1.-pyrope_proportion[i], 0.0 ])
    garnet.set_composition(composition)
    garnet_excess_volume[i]=garnet.V_excess

pressure=1.e9 # Pa
temperature=573.15 # K
composition=[0.9, 0.0, 0.1, 0.0]
garnet.set_state(pressure, temperature, composition)

# Excess gibbs for the pyrope-grossular join
n=100
pyrope_proportion= np.empty(shape=(n+1))
garnet_excess_gibbs= np.empty(shape=(n+1))
for i in range(n+1):
    pyrope_proportion[i]=float(i)/n
    composition=([ pyrope_proportion[i], 0.0, 1.-pyrope_proportion[i], 0.0 ])
    garnet.set_state(pressure, temperature, composition)
    garnet_excess_gibbs[i]=garnet.gibbs_excess


import matplotlib.pyplot as plt
plt.subplot(1,2,1)
plt.plot(pyrope_proportion,garnet_excess_volume,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=0)
plt.xlim(min(pyrope_proportion),max(pyrope_proportion))
plt.xlabel("p(pyrope)")
plt.title("V excess (m^3/mol) \nfor pyrope-grossular garnets")

plt.subplot(1,2,2)
plt.plot(pyrope_proportion,garnet_excess_gibbs,color='r',linestyle='-',marker='o',markerfacecolor='r',markersize=0)
plt.xlim(min(pyrope_proportion),max(pyrope_proportion))
plt.xlabel("p(pyrope)")
plt.title("Excess Gibbs (J/mol) for pyrope-grossular garnets\n"+str(pressure/1.e9)+" GPa, "+str(temperature)+" K")
plt.show()

