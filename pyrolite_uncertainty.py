# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


import numpy as np, matplotlib.pyplot as plt
import numpy.random
import burnman
from burnman import minerals
from burnman.minerals_base import helper_solid_solution
from matplotlib.colors import colorConverter
from scipy import interpolate
import matplotlib.colors

def normal(loc = 0.0, scale = 1.0):
    if scale <= 0.0:
        return 0.0
    else:
        return numpy.random.normal(loc, scale)

def realize_mineral( mineral ):
    #some of the minerals are missing an uncertainty for this.  Assign a large(ish) nominal value for those
    if mineral.uncertainties['err_Gprime_0'] == 0.0:
      mineral.uncertainties['err_Gprime_0'] = 0.2

    # sample the uncertainties for all the relevant parameters.  Assume that molar mass, V0, and n are well known
    mineral.params['K_0'] = mineral.params['K_0'] + normal(scale=mineral.uncertainties['err_K_0'])
    mineral.params['Kprime_0'] = mineral.params['Kprime_0'] + normal(scale=mineral.uncertainties['err_Kprime_0'])
    mineral.params['G_0'] = mineral.params['G_0'] + normal(scale=mineral.uncertainties['err_G_0'])
    mineral.params['Gprime_0'] = mineral.params['Gprime_0'] + normal(scale=mineral.uncertainties['err_Gprime_0'])
    mineral.params['Debye_0'] = mineral.params['Debye_0'] + normal(scale=mineral.uncertainties['err_Debye_0'])
    mineral.params['grueneisen_0'] = mineral.params['grueneisen_0'] + normal(scale=mineral.uncertainties['err_grueneisen_0'])
    mineral.params['q_0'] = mineral.params['q_0'] + normal(scale=mineral.uncertainties['err_q_0'])
    mineral.params['eta_s_0'] = mineral.params['eta_s_0'] + normal(scale=mineral.uncertainties['err_eta_s_0'])
    return mineral

def realize_pyrolite():

    x_pv = 0.67
    x_fp = 0.33
    pv_fe_num = 0.1
    fp_fe_num = 0.3 

    mg_perovskite = minerals.SLB_2011_ZSB_2013.mg_perovskite(); realize_mineral(mg_perovskite)
    fe_perovskite = minerals.SLB_2011_ZSB_2013.fe_perovskite(); realize_mineral(fe_perovskite)
    wuestite = minerals.SLB_2011_ZSB_2013.wuestite(); realize_mineral(wuestite)
    periclase = minerals.SLB_2011_ZSB_2013.periclase(); realize_mineral(periclase)

    perovskite = helper_solid_solution( [ mg_perovskite, fe_perovskite], [1.0-pv_fe_num, pv_fe_num])
    ferropericlase = helper_solid_solution( [ periclase, wuestite], [1.0-fp_fe_num, fp_fe_num])

    pyrolite = burnman.composite( [ (perovskite, x_pv), (ferropericlase, x_fp) ] )
    pyrolite.set_method('slb3')

    anchor_temperature = normal(loc = 1935.0, scale = 150.0)

    return pyrolite, anchor_temperature
    

#set up the seismic model
seismic_model = burnman.seismic.prem()
npts = 10
depths = np.linspace(850e3,2700e3, npts)
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

n_realizations = 1000

pressures_sampled = np.linspace(pressure[0], pressure[-1], 20*len(pressure)) 
density_list = np.empty(0)
pressure_list = np.empty(0)
vp_list = np.empty(0)
vs_list = np.empty(0)

for i in range(n_realizations):

  print "realization", i+1
  try:
    pyrolite, anchor_temperature= realize_pyrolite()
    temperature = burnman.geotherm.adiabatic(pressure, anchor_temperature, pyrolite)
    
    rho, vp, vs, _, _, _ = \
      burnman.velocities_from_rock(pyrolite, pressure, temperature, burnman.averaging_schemes.hashin_shtrikman_average())

    pressure_list = np.concatenate((pressure_list, pressures_sampled))
    fr = interpolate.interp1d(pressure, rho) 
    density_list = np.concatenate((density_list, fr(pressures_sampled)))
    fs = interpolate.interp1d(pressure, vs) 
    vs_list = np.concatenate((vs_list, fs(pressures_sampled)))
    fp = interpolate.interp1d(pressure, vp) 
    vp_list = np.concatenate((vp_list, fp(pressures_sampled)))
  except ValueError:
    print "failed, skipping"
  
f=open('output_pyrolite_uncertainty.txt','wb')
f.write("#pressure\t Vs \t Vp \t rho \n")
data=zip(pressure_list, vs_list, vp_list, density_list)
np.savetxt(f,data,fmt='%.10e',delimiter='\t')

density_hist,rho_xedge,rho_yedge = np.histogram2d(pressure_list, density_list, bins=len(pressures_sampled), normed = True)
vs_hist,vs_xedge,vs_yedge = np.histogram2d(pressure_list, vs_list, bins=len(pressures_sampled), normed = True)
vp_hist,vp_xedge,vp_yedge = np.histogram2d(pressure_list, vp_list, bins=len(pressures_sampled), normed = True)

vs_xedge = vs_xedge/1.e9
vp_xedge = vp_xedge/1.e9
rho_xedge = rho_xedge/1.e9
vs_yedge = vs_yedge/1.e3
vp_yedge = vp_yedge/1.e3
rho_yedge = rho_yedge/1.e3

#do some setup for the figure
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
plt.subplots_adjust(wspace=0.3)

#plot v_s
plt.subplot(131, aspect='equal')
plt.xlim(vs_xedge[0], vs_xedge[-1])
plt.ylim(vs_yedge[0], vs_yedge[-1])
aspect_ratio = (vs_xedge[-1]-vs_xedge[0])/(vs_yedge[-1]-vs_yedge[0])
plt.imshow(vs_hist.transpose(), origin='low', cmap='Blues',  interpolation='gaussian', \
           aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('S Wave Speed (km/s)')

#plot v_p
plt.subplot(132)
plt.xlim(vp_xedge[0], vp_xedge[-1])
plt.ylim(vp_yedge[0], vp_yedge[-1])
aspect_ratio = (vp_xedge[-1]-vp_xedge[0])/(vp_yedge[-1]-vp_yedge[0])
plt.imshow(vp_hist.transpose(), origin='low', cmap='Blues', interpolation='gaussian',\
           aspect=aspect_ratio, extent=[vp_xedge[0], vp_xedge[-1], vp_yedge[0], vp_yedge[-1]])
plt.plot(pressure/1.e9,seis_vp/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('P Wave Speed (km/s)')

#plot density
plt.subplot(133)
plt.xlim(rho_xedge[0], rho_xedge[-1])
plt.ylim(rho_yedge[0], rho_yedge[-1])
aspect_ratio = (rho_xedge[-1]-rho_xedge[0])/(rho_yedge[-1]-rho_yedge[0])
plt.imshow(density_hist.transpose(), origin='low', cmap='Blues', interpolation = 'gaussian',\
           aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Density (kg/m$^3$)')

#save and show the image
fig = plt.gcf()
fig.set_size_inches(12.0, 4.0)
fig.savefig("pyrolite_uncertainty.pdf",bbox_inches='tight', dpi=100)
'''
y_samples = np.linspace(3.e3,15.e3,200)
density_hist,rho_xedge,rho_yedge = np.histogram2d(pressure_list, density_list, bins=[pressures_sampled,y_samples], normed = True)
vs_hist,vs_xedge,vs_yedge = np.histogram2d(pressure_list, vs_list, bins=[pressures_sampled,y_samples], normed = True)
vp_hist,vp_xedge,vp_yedge = np.histogram2d(pressure_list, vp_list, bins=[pressures_sampled,y_samples], normed = True)

region=[pressures_sampled[0],pressures_sampled[-1], y_samples[0], y_samples[-1]]
aspect_ratio = (region[1]-region[0])/(region[3]-region[2])

#do some setup for the figure
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
plt.subplots_adjust(wspace=0.3)

cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap',['white','green'],256)
cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list('my_cmap2',['white','red'],256)
cmap1.set_under(color='w', alpha=1.0)
cmap2.set_under(color='w', alpha=1.0)

plt.subplot(111)
plt.xlim(region[0], region[1])
plt.ylim(region[2], region[3])
plt.imshow(vs_hist.transpose(), cmap='Purples', alpha=0.9,  origin='low', interpolation='gaussian', \
           aspect=aspect_ratio, extent=region)

plt.imshow(vp_hist.transpose(), cmap='Reds', alpha=0.9, origin='low', interpolation='gaussian',\
           aspect=aspect_ratio, extent=region)
plt.plot(pressure/1.e9,seis_vp/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')

plt.imshow(density_hist.transpose(), origin='low', cmap='Blues', alpha=0.9, interpolation = 'gaussian',\
           aspect=aspect_ratio, extent=region)
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')


plt.xlabel('Pressure (GPa)')
plt.ylabel('S Wave Speed (km/s)')

#save and show the image
fig = plt.gcf()
fig.set_size_inches(10.0, 10.0)
fig.savefig("pyrolite_uncertainty.pdf",bbox_inches='tight', dpi=100)
'''
