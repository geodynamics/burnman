# BurnMan - a lower mantle toolkit
# Copyright (C) 2014, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
    
example_uncertainties
----------------

This example demonstrates how to perturb minerals with given errors to obtain the spread of possible seimsic velocities.

*Uses:*

* :func:`burnman.mineral.perturb`


*Demonstrates:*

* Including uncertainties on mineral parameters.


"""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import normal
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap

if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))


import burnman
from burnman import minerals




# Again, we set up the seismic model and ask for its properties
# in the lower mantle.  Basically the same thing as in steps 1 and 2.

seismic_model = burnman.seismic.PREM()
min_depth = 850.e3
max_depth = 2800.e3
n_depths = 10
depths = np.linspace(min_depth, max_depth, n_depths)
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
pressures_sampled = np.linspace(pressure[0], pressure[-1], 20*len(pressure))

temperature = burnman.geotherm.brown_shankland(pressure)


"""
Finally, we make 1000 realizations of the rock and save the seismic profiles
for later visualization.  
"""

n_realizations = 1000
outfile = open('uncertainty.dat', 'w')

rock=burnman.Composite([0.8,0.2],[minerals.SLB_2011.mg_fe_perovskite(.1),minerals.SLB_2011.ferropericlase(.2)])
rock.set_method('slb3')

for i in range(n_realizations):

    print "realization", i+1
    try:
      #setup a rock

      # Perturb the rock. The default treats the given error as one sigma.
      rock.perturb()
      # Calculate the wavespeed profiles, just as before.
      rho, vp, vs, vphi, K, G = \
        burnman.velocities_from_rock(rock, pressure, temperature, burnman.averaging_schemes.VoigtReussHill())

      # This interpolates the resulting densities and wavespeeds\
      # to a higher resolution line to make the plot look nicer.
      func_rho = interpolate.interp1d(pressure, rho)
      func_vs = interpolate.interp1d(pressure, vs)
      func_vphi = interpolate.interp1d(pressure, vphi)

      pressure_list = pressures_sampled
      density_list = func_rho(pressures_sampled)
      vs_list = func_vs(pressures_sampled)
      vphi_list = func_vphi(pressures_sampled)

      # Save the output to a file
      data = zip(pressure_list, vs_list, vphi_list, density_list)
      np.savetxt(outfile,data,fmt='%.10e',delimiter='\t')

    # It is possible for the Birch-Murnaghan equation of state to go unstable for
    # some values of the parameters, which can make it fail.  If that happens, we
    # simply disregard this realization of the rock.
    except ValueError:
      print "failed, skipping"


"""
The rest of this script is concerned with plotting the results so that 
you may visualize the uncertainty space that we have sampled with our
mineral model.  In different colors we see 2D histograms, where the color 
intensity corresponds to how many of the models went through that portion of
the velocity/density space.  The dashed lines show the PREM values.  Tighter
histograms indicate that the perturbations to the mineral physical values 
make less of a difference, more spread out ones show larger effects.

In general, this plotting is more complex than the ones in step 1 and step 2
and you may safely ignore it.  
"""

# Read in the data from uncertainty.dat
outfile.close()
data = np.loadtxt('uncertainty.dat')
pressure_list = data[:,0]
vs_list = data[:,1]
vphi_list = data[:,2]
density_list = data[:,3]


# Create 2D histograms showing the spread of the data
density_hist,rho_xedge,rho_yedge = np.histogram2d(pressure_list, density_list, bins=len(pressures_sampled), normed = True)
vs_hist,vs_xedge,vs_yedge = np.histogram2d(pressure_list, vs_list, bins=len(pressures_sampled), normed = True)
vphi_hist,vphi_xedge,vphi_yedge = np.histogram2d(pressure_list, vphi_list, bins=len(pressures_sampled), normed = True)


vs_xedge = vs_xedge/1.e9
vphi_xedge = vphi_xedge/1.e9
rho_xedge = rho_xedge/1.e9
vs_yedge = vs_yedge/1.e3
vphi_yedge = vphi_yedge/1.e3
rho_yedge = rho_yedge/1.e3

left_edge = min(vs_xedge[0], vphi_xedge[0], rho_xedge[0])
right_edge = max(vs_xedge[-1], vphi_xedge[-1], rho_xedge[-1])
bottom_edge= 4.3
top_edge=11.3
aspect_ratio = (right_edge-left_edge)/(top_edge-bottom_edge)
gamma = 0.5  #Mess with this to change intensity of colormaps near the edges

plt.subplot(111, aspect='equal')
plt.xlim(left_edge, right_edge)
plt.ylim(bottom_edge, top_edge)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Wave Speed (km/s)')

#plot density
density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#edf8e9'), (0.4, '#bae4b3'), (0.6, '#74c476'), (0.8, '#31a354'), (1.0, '#006d2c') ] , gamma=0.1)
c.set_bad('w', alpha=1.0)
plt.imshow(density_hist.transpose(), origin='low', cmap=c, interpolation = 'gaussian', alpha=.7,\
       aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="--",color='g',linewidth=2.0,label='Density')


#plot v_s
vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#eff3ff'), (0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c') ] , gamma=0.5)
c.set_bad('w', alpha=1.0)
plt.imshow(vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', alpha=.7,\
       aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="--",color='b',linewidth=2.0,label='Vs')

#plot v_phi
vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#fee5d9'), (0.4, '#fcae91'), (0.6, '#fb6a4a'), (0.8, '#de2d26'), (1.0, '#a50f15') ] , gamma=0.5)
c.set_bad('w', alpha=1.0)
plt.imshow(vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7, \
       aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
plt.plot(pressure/1.e9,seis_vphi/1.e3,linestyle="--",color='r',linewidth=2.0,label='Vphi')

plt.legend(loc = 'lower right')

plt.show()
