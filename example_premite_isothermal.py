# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example is under construction.

requires:

teaches:


"""

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
        sys.path.insert(1,os.path.abspath('..')) 

import burnman
from burnman import minerals

import pymc

seismic_model = burnman.seismic.prem() # pick from .prem() .slow() .fast() (see code/seismic.py)
number_of_points = 20 #set on how many depth slices the computations should be done
depths = np.linspace(750.e3,2890.e3, number_of_points)
seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

print "preparations done"


def calc_velocities(ref_rho, ref_K, K_prime, ref_mu, mu_prime): 

    test = burnman.minerals_base.material()

    test.params['ref_V'] = 10.e-6
    test.params['molar_mass'] = ref_rho*test.params['ref_V']
    test.params['ref_K'] = ref_K
    test.params['K_prime'] = K_prime
    test.params['ref_mu'] = ref_mu
    test.params['mu_prime'] = mu_prime

    rock = burnman.composite( [(test, 1.0 )] )
    rock.set_method('bm3')

    temperature = np.empty_like(seis_p)
    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_mu = burnman.velocities_from_rock(rock,seis_p, temperature)	

    return mat_rho, mat_vphi, mat_vs

def error(ref_rho, ref_K, K_prime, ref_mu, mu_prime):
    rho, vphi, vs = calc_velocities(ref_rho, ref_K, K_prime, ref_mu, mu_prime)

    vphi_chi = burnman.chi_factor(vphi, seis_vphi)
    vs_chi = burnman.chi_factor(vs, seis_vs)
    rho_chi = burnman.chi_factor(rho, seis_rho)

    return rho_chi+vphi_chi+vs_chi


# Priors on unknown parameters:
ref_rho = pymc.Uniform('ref_rho', lower=3300., upper=4500.)
ref_K = pymc.Uniform('ref_K', lower=200.e9, upper=300.e9)
K_prime = pymc.Uniform('K_prime', lower=3., upper=6.)
ref_mu = pymc.Uniform('ref_mu', lower=50.e9, upper=250.e9)
mu_prime = pymc.Uniform('mu_prime', lower=0., upper=3.)


minerr = 1e100
@pymc.deterministic
def theta(p1=ref_rho,p2=ref_K,p3=K_prime,p4=ref_mu,p5=mu_prime):
    global minerr
    if (p1<0 or p2<0 or p3<0 or p4<0 or p5 < 0):
        return 1e30
    try:
        e = error(p1,p2,p3,p4,p5)
        if (e<minerr):
            minerr=e
            print "best fit", e, "values:", p1,p2/1.e9,p3,p4/1.e9,p5
        return e
    except ValueError:
        return 1e20


sig = 1e-4
misfit = pymc.Normal('d',mu=theta,tau=1.0/(sig*sig),value=0,observed=True,trace=True)
model = dict(ref_rho=ref_rho, ref_K=ref_K, K_prime=K_prime, ref_mu=ref_mu, mu_prime=mu_prime, misfit=misfit)
things = ['ref_rho', 'ref_K', 'K_prime', 'ref_mu', 'mu_prime']

S = pymc.MAP(model)
S.fit( method = 'fmin')

rho, vphi, vs = calc_velocities(S.ref_rho.value, S.ref_K.value, S.K_prime.value, S.ref_mu.value, S.mu_prime.value)

plt.subplot(2,2,1)
plt.plot(seis_p/1.e9,vs/1000.,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
plt.plot(seis_p/1.e9,seis_vs/1000.,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
plt.ylim([4, 8])
plt.title("Vs (km/s)")

plt.subplot(2,2,2)
plt.plot(seis_p/1.e9,vphi/1000.,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4)
plt.plot(seis_p/1.e9,seis_vphi/1000.,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4)
plt.ylim([7, 12])
plt.title("Vphi (km/s)")
   
    
# plot density
plt.subplot(2,2,3)
plt.plot(seis_p/1.e9,rho/1000.,color='r',linestyle='-',marker='^',markerfacecolor='r',markersize=4,label='model 1')
plt.plot(seis_p/1.e9,seis_rho/1000.,color='k',linestyle='-',marker='v',markerfacecolor='k',markersize=4,label='ref')
plt.title("density (kg/m^3)")
plt.legend(loc='upper left')
plt.ylim([3, 7 ])
plt.show()

print "done"


