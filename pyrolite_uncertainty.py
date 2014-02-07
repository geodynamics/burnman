# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, 2014, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


import os.path
import numpy as np, matplotlib.pyplot as plt
import numpy.ma as ma
import numpy.random
import burnman
from burnman import minerals
from burnman.minerals_base import helper_solid_solution
import matplotlib.cm
import matplotlib.colors
from scipy import interpolate

def normal(loc = 0.0, scale = 1.0):
    if scale <= 0.0:
        return 0.0
    else:
        return numpy.random.normal(loc, scale)

def realize_mineral( mineral ):
    #some of the minerals are missing an uncertainty for this.  Assign a characteristic nominal value for those
    if mineral.uncertainties['err_Gprime_0'] == 0.0:
      mineral.uncertainties['err_Gprime_0'] = 0.1

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
 
    #approximate four component pyrolite model
    x_pv = 0.67
    x_fp = 0.33
    pv_fe_num = 0.07
    fp_fe_num = 0.2 

    mg_perovskite = minerals.SLB_2011_ZSB_2013.mg_perovskite(); realize_mineral(mg_perovskite)
    fe_perovskite = minerals.SLB_2011_ZSB_2013.fe_perovskite(); realize_mineral(fe_perovskite)
    wuestite = minerals.SLB_2011_ZSB_2013.wuestite(); realize_mineral(wuestite)
    periclase = minerals.SLB_2011_ZSB_2013.periclase(); realize_mineral(periclase)

    perovskite = helper_solid_solution( [ mg_perovskite, fe_perovskite], [1.0-pv_fe_num, pv_fe_num])
    ferropericlase = helper_solid_solution( [ periclase, wuestite], [1.0-fp_fe_num, fp_fe_num])

    pyrolite = burnman.composite( [ (perovskite, x_pv), (ferropericlase, x_fp) ] )
    pyrolite.set_method('slb3')

    anchor_temperature = normal(loc = 1935.0, scale = 200.0)

    return pyrolite, anchor_temperature


def output_rock( rock, file_handle ):
  for ph in rock.staticphases:
    if( isinstance(ph.mineral, burnman.minerals_base.helper_solid_solution) ):
      for min in ph.mineral.base_materials:
        file_handle.write( '\t' + min.to_string() + '\n')
        for key in min.params:
          file_handle.write('\t\t' + key + ': ' + str(min.params[key]) + '\n')
    else:
      file_handle.write( '\t' + ph.mineral.to_string() + '\n' )
      for key in ph.mineral.params:
        file_handle.write('\t\t' + key + ': ' + str(ph.mineral.params[key]) + '\n')
    

#set up the seismic model
seismic_model = burnman.seismic.prem()
npts = 10
depths = np.linspace(850e3,2700e3, npts)
pressure, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

n_realizations = 10000
min_error = np.inf

pressures_sampled = np.linspace(pressure[0], pressure[-1], 20*len(pressure)) 
fname = 'output_pyrolite_uncertainty.txt'

if os.path.isfile(fname) == False:
  outfile = open(fname, 'w')
  outfile.write("#pressure\t Vs \t Vp \t rho \n")
  best_fit_file = open('output_pyrolite_closest_fit.txt', 'w') 

  for i in range(n_realizations):

    print "realization", i+1
    try:
      #create the ith model
      pyrolite, anchor_temperature= realize_pyrolite()
      temperature = burnman.geotherm.adiabatic(pressure, anchor_temperature, pyrolite)
     
      #calculate the seismic observables
      rho, vp, vs, vphi, K, G = \
        burnman.velocities_from_rock(pyrolite, pressure, temperature, burnman.averaging_schemes.hashin_shtrikman_average())

      #estimate the misfit with the seismic model 
      err_rho, err_vphi, err_vs = burnman.compare_l2(depths/np.mean(depths), vs/np.mean(seis_vs), vphi/np.mean(seis_vphi), \
        rho/np.mean(seis_rho), seis_vs/np.mean(seis_vs), seis_vphi/np.mean(seis_vphi), seis_rho/np.mean(seis_rho))
      error = np.sum([err_rho, err_vphi, err_vs])
      if error < min_error:
        min_error = error
        print error
        best_fit_file.write('Current best fit : '+str(error) + '\n' )
        output_rock(pyrolite, best_fit_file)
      

      #interpolate to a higher resolution line
      frho = interpolate.interp1d(pressure, rho) 
      fs = interpolate.interp1d(pressure, vs) 
      fphi = interpolate.interp1d(pressure, vphi) 

      pressure_list = pressures_sampled
      density_list = frho(pressures_sampled)
      vs_list = fs(pressures_sampled)
      vphi_list = fphi(pressures_sampled)

      

      data=zip(pressure_list, vs_list, vphi_list, density_list)
      np.savetxt(outfile,data,fmt='%.10e',delimiter='\t')

    except ValueError:
      print "failed, skipping"
  outfile.close()
  best_fit_file.close()
  
infile=open(fname,'r')
data = np.loadtxt(fname, skiprows=1)
pressure_list = data[:,0]
density_list = data[:,3]
vs_list = data[:,1]
vphi_list = data[:,2]
infile.close()

density_hist,rho_xedge,rho_yedge = np.histogram2d(pressure_list, density_list, bins=len(pressures_sampled), normed = True)
vs_hist,vs_xedge,vs_yedge = np.histogram2d(pressure_list, vs_list, bins=len(pressures_sampled), normed = True)
vphi_hist,vphi_xedge,vphi_yedge = np.histogram2d(pressure_list, vphi_list, bins=len(pressures_sampled), normed = True)


vs_xedge = vs_xedge/1.e9
vphi_xedge = vphi_xedge/1.e9
rho_xedge = rho_xedge/1.e9
vs_yedge = vs_yedge/1.e3
vphi_yedge = vphi_yedge/1.e3
rho_yedge = rho_yedge/1.e3

'''
#do some setup for the figure
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
plt.subplots_adjust(wspace=0.3)
c = matplotlib.colors.LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#eff3ff'), (0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c') ])
c.set_bad('w', alpha=1.0)

#plot v_s
plt.subplot(131, aspect='equal')
plt.xlim(vs_xedge[0], vs_xedge[-1])
plt.ylim(vs_yedge[0], vs_yedge[-1])
aspect_ratio = (vs_xedge[-1]-vs_xedge[0])/(vs_yedge[-1]-vs_yedge[0])
vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
plt.imshow(vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', \
           aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('S Wave Speed (km/s)')

#plot v_phi
plt.subplot(132)
plt.xlim(vphi_xedge[0], vphi_xedge[-1])
plt.ylim(vphi_yedge[0], vphi_yedge[-1])
aspect_ratio = (vphi_xedge[-1]-vphi_xedge[0])/(vphi_yedge[-1]-vphi_yedge[0])
vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
plt.imshow(vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', \
           aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
plt.plot(pressure/1.e9,seis_vphi/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Bulk Sound Speed (km/s)')

#plot density
plt.subplot(133)
plt.xlim(rho_xedge[0], rho_xedge[-1])
plt.ylim(rho_yedge[0], rho_yedge[-1])
aspect_ratio = (rho_xedge[-1]-rho_xedge[0])/(rho_yedge[-1]-rho_yedge[0])
density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
plt.imshow(density_hist.transpose(), origin='low', cmap=c, interpolation = 'gaussian',\
           aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="-",color='k',linewidth=2.0,label='PREM')
plt.xlabel('Pressure (GPa)')
plt.ylabel('Density (kg/m$^3$)')

#save and show the image
fig = plt.gcf()
fig.set_size_inches(12.0, 4.0)
#fig.savefig("pyrolite_uncertainty.pdf",bbox_inches='tight', dpi=100)
plt.show()
'''

left_edge = min(vs_xedge[0], vphi_xedge[0], rho_xedge[0])
right_edge = max(vs_xedge[-1], vphi_xedge[-1], rho_xedge[-1])
bottom_edge= 4.3
top_edge=11.3
aspect_ratio = (right_edge-left_edge)/(top_edge-bottom_edge)
gamma = 0.8  #Mess with this to change intensity of colormaps near the edges

#do some setup for the figure
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'
plt.rc('font', family='sanserif')
plt.subplots_adjust(wspace=0.3)

plt.subplot(111, aspect='equal')
plt.xlim(left_edge, right_edge)
plt.ylim(bottom_edge, top_edge)
plt.xlabel('Pressure (GPa)')
plt.ylabel('Wave Speed (km/s)')


#plot v_s
vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
c = matplotlib.colors.LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#eff3ff'), (0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c') ] , gamma=gamma)
c.set_bad('w', alpha=1.0)
plt.imshow(vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', alpha=.7,\
           aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
plt.plot(pressure/1.e9,seis_vs/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')

#plot v_phi
vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
c = matplotlib.colors.LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#fee5d9'), (0.4, '#fcae91'), (0.6, '#fb6a4a'), (0.8, '#de2d26'), (1.0, '#a50f15') ] , gamma=gamma)
c.set_bad('w', alpha=1.0)
plt.imshow(vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7, \
           aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
plt.plot(pressure/1.e9,seis_vphi/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')

#plot density
density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
c = matplotlib.colors.LinearSegmentedColormap.from_list('vphi', [ (0, '#ffffff'), (0.2, '#edf8e9'), (0.4, '#bae4b3'), (0.6, '#74c476'), (0.8, '#31a354'), (1.0, '#006d2c') ] , gamma=gamma)
c.set_bad('w', alpha=1.0)
plt.imshow(density_hist.transpose(), origin='low', cmap=c, interpolation = 'gaussian', alpha=.7,\
           aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
plt.plot(pressure/1.e9,seis_rho/1.e3,linestyle="--",color='k',linewidth=2.0,label='PREM')


#save and show the image
fig = plt.gcf()
fig.set_size_inches(6.0, 6.0)
fig.savefig("pyrolite_uncertainty.pdf",bbox_inches='tight', dpi=100)
plt.show()
