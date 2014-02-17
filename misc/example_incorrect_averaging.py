# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
Attempt to reproduce Figure 6.12 From chapter 6 of Physics and Chemistry of the Deep Earth, 2013
Book chapter by Motohiko Murakami, editor Shun-ichiro Karato
"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals
import matplotlib.image as mpimg
import burnman.minerals_base as mb
import numpy as np
import colors


plt.figure(dpi=100,figsize=(12,6))
prop={'size':12}
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = '\usepackage{relsize}'

dashstyle2=(7,3)
dashstyle3=(3,2)

method = 'slb2'
    
#define the minerals from table 6.3
mg_perovskite = burnman.material()
mg_perovskite.params = {       'name': 'Mg perovskite',
                    'molar_mass' : 0.1004,
                    'V_0': 24.43e-6,
                    'K_0': 253.0e9,
                    'Kprime_0': 3.9,
                    'G_0' : 172.9e9,
                    'Gprime_0' : 1.56,
                    'n': 5.0,
                    'Debye_0': 1100.,
                    'grueneisen_0': 1.40,
                    'q_0': 1.40, 
                    'eta_s_0' : 2.6}
mg_perovskite.set_method('slb2')

fe_perovskite = burnman.material()
fe_perovskite.params = {       'name': 'Fe perovskite',
                    'molar_mass' : 0.1319,
                    'V_0': 25.49e-6,
                    'K_0': 281.0e9,
                    'Kprime_0': 4.1,
                    'G_0' : 138.0e9,
                    'Gprime_0' : 1.70,
                    'n': 5.0,
                    'Debye_0': 841.,
                    'grueneisen_0': 1.48,
                    'q_0': 1.40, 
                    'eta_s_0' : 2.1}
fe_perovskite.set_method(method)

periclase = burnman.material()
periclase.params = {       'name': 'periclase',
                    'molar_mass' : 0.0403,
                    'V_0': 11.24e-6,
                    'K_0': 161.0e9,
                    'Kprime_0': 3.9,
                    'G_0' : 130.9e9,
                    'Gprime_0' : 1.92,
                    'n': 2.0,
                    'Debye_0': 773.,
                    'grueneisen_0': 1.50,
                    'q_0': 1.50, 
                    'eta_s_0' : 2.3}
periclase.set_method(method)

wustite = burnman.material()
wustite.params = {       'name': 'wustite',
                    'molar_mass' : 0.07184,
                    'V_0': 12.06e-6,
                    'K_0': 152.0e9,
                    'Kprime_0': 4.9,
                    'G_0' : 47.0e9,
                    'Gprime_0' : 0.70,
                    'n': 2.0,
                    'Debye_0': 455.,
                    'grueneisen_0': 1.28,
                    'q_0': 1.50, 
                    'eta_s_0' : 0.8}
wustite.set_method(method)


#in the text for the book chapter a linear relationship in elastic properties
#for the solid solutions is assumed...
class ferropericlase(mb.helper_solid_solution):
    def __init__(self, fe_num):
        base_materials = [periclase, wustite]
        molar_fraction = [1. - fe_num, 0.0 + fe_num]
        mb.helper_solid_solution.__init__(self, base_materials, molar_fraction)



class perovskite(mb.helper_solid_solution):
    def __init__(self, fe_num):
        base_materials = [mg_perovskite, fe_perovskite]
        molar_fraction = [1. - fe_num, 0.0 + fe_num]
        mb.helper_solid_solution.__init__(self, base_materials, molar_fraction)


#define the P-T path
pressure = np.linspace(28.0e9, 129e9, 25.)
temperature_bs = burnman.geotherm.brown_shankland(pressure)
temperature_an = burnman.geotherm.anderson(pressure)

#seismic model for comparison:
seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see burnman/seismic.py)
depths = map(seismic_model.depth, pressure)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)    

#pure perovskite
perovskitite = burnman.composite( ( (perovskite(0.06), 1.0),) )
perovskitite.set_method(method)
 
#pure periclase
periclasite = burnman.composite( ( (ferropericlase(0.21), 1.0),))
periclasite.set_method(method) 

#pyrolite (80% perovskite)
pyrolite = burnman.composite( ( (perovskite(0.06), 0.834),
                              (ferropericlase(0.21), 0.166) ) )
pyrolite.set_method(method) 

#preferred mixture?
amount_perovskite = 0.92
preferred_mixture = burnman.composite( ( (perovskite(0.06), amount_perovskite),
                                         (ferropericlase(0.21), 1.0-amount_perovskite) ) )
preferred_mixture.set_method(method) 
    

mat_rho_1, mat_vp_1, mat_vs_1, mat_vphi_1, mat_K_1, mat_G_1 = burnman.velocities_from_rock(perovskitite,seis_p, temperature_bs)    
mat_rho_2, mat_vp_2, mat_vs_2, mat_vphi_2, mat_K_2, mat_G_2 = burnman.velocities_from_rock(periclasite,seis_p, temperature_bs)    
mat_rho_3, mat_vp_3, mat_vs_3, mat_vphi_3, mat_K_3, mat_G_3 = burnman.velocities_from_rock(pyrolite,seis_p, temperature_bs)
mat_rho_4, mat_vp_4, mat_vs_4, mat_vphi_4, mat_K_4, mat_G_4 = burnman.velocities_from_rock(preferred_mixture,seis_p, temperature_bs)



### HERE IS THE STEP WITH THE INCORRECT MIXING ###
# comment this out to have correct phase averaging, leave it in to have incorrect phase averaging

mat_vs_3_wr = 0.5*((0.834*mat_vs_1 + 0.166*mat_vs_2) + np.ones_like(mat_vs_1)/(0.834/mat_vs_1 + 0.166/mat_vs_2))
mat_vs_4_wr = 0.5*((0.92*mat_vs_1 + 0.08*mat_vs_2) + np.ones_like(mat_vs_1)/(0.92/mat_vs_1 + 0.08/mat_vs_2))

plt.subplot(1,2,2)
plt.ylim(5.2,7.4)
plt.xlim(25,135)
#fig1 = mpimg.imread('input_figures/murakami_book_chapter.png')
#plt.imshow(fig1, extent=[25,135,5.0,7.6], aspect='auto')
plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',linestyle='-',marker='None',markerfacecolor='w',markersize=4,label='PREM',linewidth=3.0,mew=1.5)
plt.plot(seis_p/1.e9,mat_vs_1/1.e3,color=colors.color(3),marker='v',markerfacecolor=colors.color(3), \
         markersize=4, markeredgecolor=colors.color(3),linewidth=1.5,label='perovskite')
plt.plot(seis_p/1.e9,mat_vs_2/1.e3,color=colors.color(1),linestyle='-', \
         linewidth=1.5,marker='^',markerfacecolor=colors.color(1), markersize=4, \
         markeredgecolor=colors.color(1),label='periclase')
plt.plot(seis_p/1.e9,mat_vs_4_wr/1.e3,color=colors.color(4),dashes=dashstyle3, \
         linewidth=1.5,marker='o',markerfacecolor=colors.color(4), markersize=4,  \
         markeredgecolor=colors.color(4),label='92\% pv')
plt.plot(seis_p/1.e9,mat_vs_3_wr/1.e3,color='g',linestyle='-',dashes=dashstyle2, \
         linewidth=1.5,marker='o',markerfacecolor='w', markersize=4, markeredgecolor='g',label='83\% pv')
plt.legend(loc='lower right',prop={'size':12})


plt.title("Phase average on velocities")

plt.xlabel("Pressure (GPa)")

plt.subplot(1,2,1)
plt.ylim(5.2,7.4)
plt.xlim(25,135)
#fig1 = mpimg.imread('input_figures/murakami_book_chapter.png')
#plt.imshow(fig1, extent=[25,135,5.0,7.6], aspect='auto')
plt.plot(seis_p/1.e9,seis_vs/1.e3,color='k',linestyle='-',marker='None',markerfacecolor='w',markersize=4,label='PREM',linewidth=3.0,mew=1.5)
plt.plot(seis_p/1.e9,mat_vs_1/1.e3,color=colors.color(3),marker='v',markerfacecolor=colors.color(3), \
         markersize=4, markeredgecolor=colors.color(3),linewidth=1.5,label='perovskite')
plt.plot(seis_p/1.e9,mat_vs_2/1.e3,color=colors.color(1),linestyle='-', \
         linewidth=1.5,marker='^',markerfacecolor=colors.color(1), markersize=4, \
         markeredgecolor=colors.color(1),label='periclase')
plt.plot(seis_p/1.e9,mat_vs_4/1.e3,color=colors.color(4),dashes=dashstyle3, \
         linewidth=1.5,marker='o',markerfacecolor=colors.color(4), markersize=4,  \
         markeredgecolor=colors.color(4),label='92\% pv')
plt.plot(seis_p/1.e9,mat_vs_3/1.e3,color='g',linestyle='-',dashes=dashstyle2, \
         linewidth=1.5,marker='o',markerfacecolor='w', markersize=4, markeredgecolor='g',label='83\% pv')

plt.title(" V.-R.-H. on moduli")
plt.xlabel("Pressure (GPa)")
plt.ylabel("Shear Velocity Vs (km/s)")
plt.tight_layout()
plt.savefig("example_incorrect_averaging.pdf",bbox_inches='tight')
plt.show()




