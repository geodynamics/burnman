from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.minerals import \
    DKS_2013_liquids, \
    DKS_2013_solids, \
    SLB_2011
from burnman import constants
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.image as mpimg


fig1 = mpimg.imread('figures/MgO-SiO2_enthalpy_mixing.png')
plt.subplot(232)
plt.imshow(fig1, extent=[0, 1, -40000, 10000], aspect='auto')

fig1 = mpimg.imread('figures/MgO-SiO2_entropy_mixing.png')
plt.subplot(233)
plt.imshow(fig1, extent=[0, 1, 0, 24], aspect='auto')

fig1 = mpimg.imread('figures/MgO-SiO2_volume_mixing.png')
plt.subplot(234)
plt.imshow(fig1, extent=[0, 1, -2.0e-6, 0.4e-6], aspect='auto')

phases = [DKS_2013_liquids.SiO2_liquid(),
         DKS_2013_liquids.MgSiO3_liquid(),
         DKS_2013_liquids.MgSi2O5_liquid(),
         DKS_2013_liquids.MgSi3O7_liquid(),
         DKS_2013_liquids.MgSi5O11_liquid(),
         DKS_2013_liquids.Mg2SiO4_liquid(),
         DKS_2013_liquids.Mg3Si2O7_liquid(),
         DKS_2013_liquids.Mg5SiO7_liquid(),
         DKS_2013_liquids.MgO_liquid()
         ]


pressure = 25.e9 # Pa
temperature = 3000. # K

MgO_liq = DKS_2013_liquids.MgO_liquid()
SiO2_liq = DKS_2013_liquids.SiO2_liquid()

pts = [[0, 3000],
       [5, 3000],
       [10, 3000],
       [25, 3000],
       [50, 4000],
       [75, 5000],
       [100, 6000],
       [135, 6000]]



for p, t in pts:
    pressure=p*1.e9
    temperature=t*1.

    MgO_liq.set_state(pressure, temperature)
    SiO2_liq.set_state(pressure, temperature)
    MgO_gibbs = MgO_liq.gibbs
    SiO2_gibbs = SiO2_liq.gibbs

    MgO_H = MgO_liq.H
    SiO2_H = SiO2_liq.H

    MgO_S = MgO_liq.S
    SiO2_S = SiO2_liq.S

    MgO_V = MgO_liq.V
    SiO2_V = SiO2_liq.V

    MgO_K_T = MgO_liq.K_T
    SiO2_K_T = SiO2_liq.K_T

    fSis=[]
    Gexs=[]
    Hexs=[]
    Sexs=[]
    Vexs=[]
    K_Ts=[]
    K_Texs=[]
    for phase in phases:
        
        nSi = phase.params['formula']['Si']
        nMg = phase.params['formula']['Mg']
        
        sum_cations = nSi+nMg
        fSi=nSi/sum_cations
        
        phase.set_state(pressure, temperature)
        Gex = phase.gibbs/sum_cations - (fSi*SiO2_gibbs + (1.-fSi)*MgO_gibbs)       
        Hex = phase.H/sum_cations - (fSi*SiO2_H + (1.-fSi)*MgO_H)

        Sex = phase.S/sum_cations - (fSi*SiO2_S + (1.-fSi)*MgO_S)

        Vex = phase.V/sum_cations - (fSi*SiO2_V + (1.-fSi)*MgO_V)

        K_T = phase.K_T
        K_Tex = (phase.K_T - (fSi*SiO2_K_T + (1.-fSi)*MgO_K_T))/K_T

        fSis.append(fSi)
        Gexs.append(Gex)
        Hexs.append(Hex)
        Sexs.append(Sex)
        Vexs.append(Vex)
        K_Ts.append(K_T)
        K_Texs.append(K_Tex)


    plt.subplot(231)
    plt.title('Excess Gibbs') 
    plt.plot(fSis, Gexs, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')
    plt.subplot(232)
    plt.title('Excess Enthalpies') 
    plt.plot(fSis, Hexs, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')
    plt.subplot(233)
    plt.title('Excess Entropies') 
    plt.plot(fSis, Sexs, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')
    plt.subplot(234)
    plt.title('Excess Volumes') 
    plt.plot(fSis, Vexs, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')
    plt.subplot(235)
    plt.title('Fractional excess K_T') 
    plt.plot(fSis, K_Texs, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')
    plt.subplot(236)
    plt.title('K_T') 
    plt.plot(fSis, K_Ts, marker='o', linestyle='None', label=str(p)+' GPa, '+str(t)+' K')

plt.legend(loc='lower right')
plt.show()

