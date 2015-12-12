from __future__ import absolute_import
from __future__ import print_function
# Benchmarks for the solid solution class
import os.path, sys
sys.path.insert(1,os.path.abspath('../..'))

import burnman
from burnman.minerals import SLB_2011
from burnman.minerals import HP_2011_ds62
from burnman.minerals import Sundman_1991
from burnman import constants
import numpy as np


def p(v1, v2):
    return (v2-v1)/v1

###
filemin=[['SLB2011', '../../burnman/data/input_perplex/fo_SLB2011_params.dat', SLB_2011.fo()],['HP2011', '../../burnman/data/input_perplex/fo_HP2011_params.dat', HP_2011_ds62.fo()]]

for database, f, mineral in filemin:
    f = open(f, 'r')
    datalines = [ line.strip() for idx, line in enumerate(f.read().split('\n')) if line.strip() and idx>0 ]
    data = [ list(map(float,"%".join(line.split("%")[:1]).split())) for line in datalines ]
    P, T, H, S, V, C_p, alpha, beta, rho = list(zip(*data))

    variables=['H','S','V','C_p','alpha','beta','rho']
    
    fo = mineral
    percentage_diff=[]
    PT=[]

    print('Benchmarks for', database, 'database with method', fo.params['equation_of_state'])
    print(variables)

    for line in data:
        P, T, H, S, V, C_p, alpha, beta, rho = line
        fo.set_state(P*1.e5,T)
        gibbs=H-T*S
        PT.append([P/1.e4,T])
        diff=[p(fo.gibbs, gibbs), p(fo.H, H), p(fo.S, S), p(fo.V, V/1.e5), p(fo.C_p, C_p), p(fo.alpha, alpha), p(fo.K_T, 1.e5/beta), p(fo.density, rho)]
        percentage_diff.append(diff)

    percentage_diff=np.array(percentage_diff)
    i,j = np.unravel_index(percentage_diff.argmax(), percentage_diff.shape)

    print('Maximum error in', database, 'database:')
    print(variables[j], ':', percentage_diff[i,j], '% at', PT[i][0], 'GPa and', PT[i][1], 'K')
    print('')





variables=['V','beta','rho']
    
fo = HP_2011_ds62.fo()
fo.set_method('mt')
percentage_diff=[]
PT=[]

print('Benchmarks for', database, 'database with method', fo.params['equation_of_state'])
print(variables)


perplex_output=[[1., 4.3660, 0.77818E-06, 3222.4],[50000., 4.2104,  0.67868E-06,   3341.5],[100000., 4.0778, 0.60406E-06,   3450.2]]
T=298.15
for P, V, beta, rho in perplex_output:
    fo.set_state(P*1.e5,T)
    PT.append([P/1.e4,T])
    diff=[p(fo.V, V/1.e5), p(fo.K_T, 1.e5/beta), p(fo.density, rho)]
    print(diff)
    percentage_diff.append(diff)

percentage_diff=np.array(percentage_diff)
i,j = np.unravel_index(percentage_diff.argmax(), percentage_diff.shape)

print('Maximum error in', database, 'database:')
print(variables[j], ':', percentage_diff[i,j], '% at', PT[i][0], 'GPa and', PT[i][1], 'K')
print('')


'''
Sundman (1991) models for fcc and bcc iron
'''


def magnetic_gibbs(T, Tc, beta, p):
    A = (518./1125.) + (11692./15975.)*((1./p) - 1.)
    tau=T/Tc
    if tau < 1: 
        f=1.-(1./A)*(79./(140.*p*tau) + (474./497.)*(1./p - 1.)*(np.power(tau, 3.)/6. + np.power(tau, 9.)/135. + np.power(tau, 15.)/600.))
    else:
        f=-(1./A)*(np.power(tau,-5)/10. + np.power(tau,-15)/315. + np.power(tau, -25)/1500.)
    return constants.gas_constant*T*np.log(beta + 1.)*f

def HSERFe(T):
    if T < 1811:
        gibbs=1224.83 + 124.134*T - 23.5143*T*np.log(T) - 0.00439752*T*T - 5.89269e-8*T*T*T + 77358.3/T
    else:
        gibbs=-25384.451 + 299.31255*T - 46.*T*np.log(T) + 2.2960305e31*np.power(T,-9.)
    return gibbs 

def gibbs_bcc_1bar(T):
    Tc=1043.
    beta=2.22
    p=0.4
    return HSERFe(T) + magnetic_gibbs(T, Tc, beta, p)

def gibbs_fcc_1bar(T):
    Tc=201.
    beta=2.10
    p=0.28
    if T < 1811:
        gibbs=HSERFe(T) - 1462.4 + 8.282*T - 1.15*T*np.log(T) + 0.00064*T*T
    else:
        gibbs= - 27098.266 + 300.25256*T - 46.*T*np.log(T) + 2.78854e31*np.power(T,-9.)
    return gibbs + magnetic_gibbs(T, Tc, beta, p)

bcc=Sundman_1991.bcc_iron()
fcc=Sundman_1991.fcc_iron()


temperatures=np.linspace(298.15, 1798.15, 101)
Pr=1.e5

fcc_gibbs=[]
bcc_gibbs=[]
fcc_gibbs_model=[]
bcc_gibbs_model=[]
for T in temperatures:
    fcc.set_state(Pr, T)
    bcc.set_state(Pr, T)

    fcc_gibbs.append(fcc.gibbs)
    bcc_gibbs.append(bcc.gibbs)

    fcc_gibbs_model.append(gibbs_fcc_1bar(T))
    bcc_gibbs_model.append(gibbs_bcc_1bar(T))

import matplotlib.pyplot as plt
plt.plot( temperatures, np.array(fcc_gibbs)-np.array(bcc_gibbs), 'r-', linewidth=1., label='FCC-BCC')
plt.plot( temperatures, np.array(fcc_gibbs_model)-np.array(bcc_gibbs_model), 'b-', linewidth=1., label='FCC-BCC model')

plt.title("Gibbs free energy difference between bcc and fcc iron")
plt.ylabel("Gibbs free energy (J/mol)")
plt.xlabel("Temperature (K)")
plt.xlim(1100, 1700)
plt.ylim(-10, 0)
plt.legend(loc='lower left')
plt.show()
