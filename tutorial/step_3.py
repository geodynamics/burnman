# BurnMan - a lower mantle toolkit
# Copyright (C) 2014, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
from numpy.random import normal
from scipy import interpolate
from matplotlib.colors import LinearSegmentedColormap

import burnman
from burnman import minerals
from burnman.mineral_helpers import HelperSolidSolution


def realize_mineral( mineral ):
    K_prime_std_dev = 0.2  #<----------------- One sigma uncertainty in K prime
    G_prime_std_dev = 0.2  #<----------------- One sigma uncertainty in G prime

    mineral.params['Kprime_0'] = mineral.params['Kprime_0'] + normal(scale=K_prime_std_dev)
    mineral.params['Gprime_0'] = mineral.params['Gprime_0'] + normal(scale=G_prime_std_dev)
    return mineral

def realize_rock():

    #approximate four component pyrolite model
    x_pv = 0.7    #<---------------------------- Fraction of perovskite in preferred model
    x_fp = 1.0-x_pv   
    pv_fe_num = 0.05   #<------------------------- Fraction of iron in perovskite
    fp_fe_num = 0.3    #<------------------------- Franction of iron in ferropericlase

    mg_perovskite = minerals.SLB_2011.mg_perovskite(); realize_mineral(mg_perovskite)
    fe_perovskite = minerals.SLB_2011.fe_perovskite(); realize_mineral(fe_perovskite)
    wuestite = minerals.SLB_2011.wuestite(); realize_mineral(wuestite)
    periclase = minerals.SLB_2011.periclase(); realize_mineral(periclase)

    perovskite = HelperSolidSolution( [ mg_perovskite, fe_perovskite], [1.0-pv_fe_num, pv_fe_num])
    ferropericlase = HelperSolidSolution( [ periclase, wuestite], [1.0-fp_fe_num, fp_fe_num])

    mantle_rock = burnman.Composite( [x_pv, x_fp], [perovskite, ferropericlase] )
    mantle_rock.set_method('slb3')

    return mantle_rock

#set up the seismic model
seismic_model = burnman.seismic.PREM()
min_depth = 850.e3
max_depth = 2800.e3
n_depths = 10
depths = np.linspace(min_depth, max_depth, n_depths)
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)
pressures_sampled = np.linspace(pressure[0], pressure[-1], 20*len(pressure))



n_realizations = 30
outfile = open('uncertainty.dat', 'w')

for i in range(n_realizations):

    print "realization", i+1
    try:
      #create the ith model
      rock = realize_rock()
      temperature = burnman.geotherm.brown_shankland(pressure)

      #calculate the seismic observables
      rho, vp, vs, vphi, K, G = \
        burnman.velocities_from_rock(rock, pressure, temperature, burnman.averaging_schemes.VoigtReussHill())

      #interpolate to a higher resolution line
      func_rho = interpolate.interp1d(pressure, rho)
      func_vs = interpolate.interp1d(pressure, vs)
      func_vphi = interpolate.interp1d(pressure, vphi)

      pressure_list = pressures_sampled
      density_list = func_rho(pressures_sampled)
      vs_list = func_vs(pressures_sampled)
      vphi_list = func_vphi(pressures_sampled)


      data = zip(pressure_list, vs_list, vphi_list, density_list)
      np.savetxt(outfile,data,fmt='%.10e',delimiter='\t')

    except ValueError:
      print "failed, skipping"

outfile.close()
data = np.loadtxt('uncertainty.dat')
pressure_list = data[:,0]
density_list = data[:,3]
vs_list = data[:,1]
vphi_list = data[:,2]

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
gamma = 0.8  #Mess with this to change intensity of colormaps near the edges

#do some setup for the figure
plt.subplots_adjust(wspace=0.3)

plt.subplot(111, aspect='equal')
plt.xlim(left_edge, right_edge)
plt.ylim(bottom_edge, top_edge)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Wave Speed (km/s)')


#plot v_s
vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#eff3ff'), (0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c') ] , gamma=gamma)
c.set_bad('w', alpha=1.0)
plt.imshow(vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', alpha=.7,\
       aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')

#plot v_phi
vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#fee5d9'), (0.4, '#fcae91'), (0.6, '#fb6a4a'), (0.8, '#de2d26'), (1.0, '#a50f15') ] , gamma=gamma)
c.set_bad('w', alpha=1.0)
plt.imshow(vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7, \
       aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
plt.plot(pressure/1.e9,seis_vphi/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')

#plot density
density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
c = LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#edf8e9'), (0.4, '#bae4b3'), (0.6, '#74c476'), (0.8, '#31a354'), (1.0, '#006d2c') ] , gamma=0.2)
c.set_bad('w', alpha=1.0)
plt.imshow(density_hist.transpose(), origin='low', cmap=c, interpolation = 'gaussian', alpha=.7,\
       aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')

plt.show()
