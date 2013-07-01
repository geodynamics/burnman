# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example is under construction.

requires:

teaches:


"""

import os, sys, numpy as np, matplotlib.pyplot as plt

import scipy.optimize as opt
import burnman
from burnman import minerals

#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
        sys.path.insert(1,os.path.abspath('..')) 


def calculate_forward_problem(frac, pressures):
    print frac
    
    rock = burnman.composite( [ (minerals.SLB2011.mg_perovskite(),frac[2]*frac[0] ),
                            (minerals.SLB2011.fe_perovskite(), frac[2]*(1.0-frac[0]) ),
                            (minerals.SLB2011.periclase(), (1.0-frac[2])) ,
                            (minerals.SLB2011.wuestite(), 0.0*(1.0-frac[2])*(1.0-frac[0]) ) ] )
    rock.set_method('slb3')
    temperature = burnman.geotherm.self_consistent(pressures, frac[1], rock)
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = burnman.velocities_from_rock(rock,pressures, temperature)
    return mat_rho,mat_vs, mat_vphi
   
def error(guess, pressures, obs_rho, obs_vs, obs_vphi):
    [rho,vs,vphi] = calculate_forward_problem(guess,  pressures )

    vs_l2 = [ (vs[i] - obs_vs[i])*(vs[i] - obs_vs[i])/(obs_vs[i]*obs_vs[i]) for i in range(len(obs_vs)) ]
    vphi_l2 = [ (vphi[i] - obs_vphi[i])*(vphi[i] - obs_vphi[i])/(obs_vphi[i]*obs_vphi[i]) for i in range(len(obs_vphi)) ]
    rho_l2=[(rho[i] - obs_rho[i])*(rho[i]-obs_rho[i])/(obs_rho[i]*obs_rho[i]) for i in range(len(obs_rho)) ]
    l2_error= sum(vphi_l2)+sum(vs_l2)+ sum(rho_l2)
    print "current error:", l2_error
    return l2_error

pressures=np.linspace(30e9,120e9,20)


prem = burnman.seismic.prem()
depths = map(prem.depth,pressures) 
seis_p, prem_density, prem_vp, prem_vs, prem_vphi = prem.evaluate_all_at(depths)

#make the mineral to fit
guess = [0.95,2100,0.65]
lowerbounds=[0.5,0, 0,0,1800]
upperbounds=[1.0,0.2,0.5,0.2,2100]

#first, do the second-order fit
func = lambda x : error( x, pressures, prem_density, prem_vs, prem_vphi)
sol = opt.fmin(func, guess, xtol=0.8)
print sol 



