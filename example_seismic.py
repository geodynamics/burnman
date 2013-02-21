# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# under GPL v2 or later.

"""
Shows the various ways to input seismic models (Vs, Vp, Vphi, Density) as a
function of depth (or P) as well as different velocity models available:
PREM (Dziewonski & Anderson, 1981)
reference model for fast regionsi (outside the LLSVP's) in the lower mantle (Lekic et al. 2012)
reference model for slow regions (LLSVP's) in the lower mantle (Lekic et la. 2012)

requires:
- burnman.seismic

teaches:
- seismic models

"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
	sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

#create a seismic dataset from prem:
s=burnman.seismic.prem()

# specify where we want to evaluate, here we map from pressure to depth, because we can
p = np.arange(1.0e9,360.0e9,5e9)
depths = map(s.depth, p) 
#we could also just specify some depth levels directly like this:
#depths = np.arange(35,5600,100)
#we could also use the data points where the seismic model is specified:
depths = s.internal_depth_list()

#now evaluate everything at the given depths levels (using interpolation)
pressures, density, v_p, v_s, v_phi = s.evaluate_all_at(depths)

# plot vs and vp and v_phi (note that v_phi is computed!)
plt.subplot(2,2,1)
plt.title('prem')
plt.plot(depths,v_p,'+-r', label='v_p')
plt.plot(depths,v_s,'+-b', label='v_s')
plt.plot(depths,v_phi,'--g', label='v_phi')
plt.legend(loc='lower left')
plt.xlabel('depth in km')
plt.ylabel('km/s')

# plot pressure,density vs depth from prem:
plt.subplot(2,2,2)
plt.title('prem')
plt.plot(depths,pressures/1e9,'-r', label='pressure')
plt.ylabel('GPa')
plt.xlabel('depth in km')
plt.legend(loc='upper left')
plt.twinx()
plt.ylabel('g/cc')
plt.plot(depths,density,'-b', label='density')
plt.legend(loc='lower right')


#now load a different seismic model:
sslow = burnman.seismic.slow()
depths2 = sslow.internal_depth_list()
pressures2, density2, v_p2, v_s2, v_phi2 = sslow.evaluate_all_at(depths2)

sfast = burnman.seismic.fast()
depths3 = sfast.internal_depth_list()
pressures3, density3, v_p3, v_s3, v_phi3 = sfast.evaluate_all_at(depths3)


plt.subplot(2,2,3)
plt.plot(pressures/1e9,v_p,'-k', label='v_p prem')
plt.plot(pressures2/1e9,v_p2,'-r', label='v_p slow')
plt.plot(pressures3/1e9,v_p3,'-b', label='v_p fast')

plt.legend(loc='lower right')
plt.xlim([30,136])
plt.ylim([11,14])
plt.xlabel('pressure')
plt.ylabel('km/s')

plt.subplot(2,2,4)
plt.plot(pressures/1e9,v_s,'-k', label='v_s prem')
plt.plot(pressures2/1e9,v_s2,'-r', label='v_s slow')
plt.plot(pressures3/1e9,v_s3,'-b', label='v_s fast')

plt.legend(loc='upper left')
plt.xlim([30,136])
plt.ylim([6,8])
plt.xlabel('pressure')
plt.ylabel('km/s')

plt.show()
