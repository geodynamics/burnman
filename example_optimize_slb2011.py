import os, sys, numpy as np, matplotlib.pyplot as plt

import scipy.optimize as opt
import burnman
from burnman import minerals

#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
        sys.path.insert(1,os.path.abspath('..')) 


def calculate_forward_problem(frac, pressures):
    print frac
    rock = burnman.composite( ( (minerals.SLB2011_mg_perovskite(),frac[0] ),
                            (minerals.SLB2011_fe_perovskite(), frac[1] ),
                            (minerals.SLB2011_periclase(), frac[2] ),
                            (minerals.SLB2011_wuestite(), frac[3] )))
    rock.set_method('slb3')
    temperature = burnman.geotherm.self_consistent(pressures, frac[4], rock)

    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.velocities_from_rock(rock,pressures, temperature)
    return mat_vs, mat_vphi
   
def error(guess, pressures, obs_vs, obs_vphi):
    [vs,vphi] = calculate_forward_problem(guess,  pressures )

    vs_l2 = [ (vs[i] - obs_vs[i])*(vs[i] - obs_vs[i]) for i in range(len(obs_vs)) ]
    vphi_l2 = [ (vphi[i] - obs_vphi[i])*(vphi[i] - obs_vphi[i]) for i in range(len(obs_vphi)) ]
    l2_error= sum(vs_l2)+sum(vphi_l2)
    return l2_error


pressures=np.linspace(30e9,120e9,20)


prem = burnman.seismic.prem()
depths = map(prem.depth,pressures) 
seis_p, prem_density, prem_vp, prem_vs, prem_vphi = prem.evaluate_all_at(depths)

#make the mineral to fit
guess = [0.62,0.078,0.36,0.022,1900]

#first, do the second-order fit
func = lambda x : error( x, pressures, prem_vs, prem_vphi)
sol = opt.fmin(func, guess, xtol=0.1)
print sol 



