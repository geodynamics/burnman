# Benchmarks for the solid solution class
import os.path, sys
sys.path.insert(1,os.path.abspath('../..'))

import burnman
from burnman.minerals import SLB_2011
from burnman.minerals import HP_2011_ds62
import numpy as np


def p(v1, v2):
    return (v2-v1)/v1

###
filemin=[['SLB2011', '../../burnman/data/input_perplex/fo_SLB2011_params.dat', SLB_2011.fo()],['HP2011', '../../burnman/data/input_perplex/fo_HP2011_params.dat', HP_2011_ds62.fo()]]

for database, f, mineral in filemin:
    f = open(f, 'r')
    datalines = [ line.strip() for idx, line in enumerate(f.read().split('\n')) if line.strip() and idx>0 ]
    data = [ map(float,"%".join(line.split("%")[:1]).split()) for line in datalines ]
    P, T, H, S, V, C_p, alpha, beta, rho = zip(*data)

    variables=['H','S','V','C_p','alpha','beta','rho']
    
    fo = mineral
    percentage_diff=[]
    PT=[]

    print 'Benchmarks for', database, 'database with method', fo.params['equation_of_state']
    print variables

    for line in data:
        P, T, H, S, V, C_p, alpha, beta, rho = line
        fo.set_state(P*1.e5,T)
        gibbs=H-T*S
        PT.append([P/1.e4,T])
        diff=[p(fo.gibbs, gibbs), p(fo.H, H), p(fo.S, S), p(fo.V, V/1.e5), p(fo.C_p, C_p), p(fo.alpha, alpha), p(fo.K_T, 1.e5/beta), p(fo.density(), rho)]
        print diff
        percentage_diff.append(diff)

    percentage_diff=np.array(percentage_diff)
    i,j = np.unravel_index(percentage_diff.argmax(), percentage_diff.shape)

    print 'Maximum error in', database, 'database:'
    print variables[j], ':', percentage_diff[i,j], '% at', PT[i][0], 'GPa and', PT[i][1], 'K'
    print ''





variables=['V','beta','rho']
    
fo = HP_2011_ds62.fo()
fo.set_method('mtait')
percentage_diff=[]
PT=[]

print 'Benchmarks for', database, 'database with method', fo.params['equation_of_state']
print variables


perplex_output=[[1., 4.3660, 0.77818E-06, 3222.4],[50000., 4.2104,  0.67868E-06,   3341.5],[100000., 4.0778, 0.60406E-06,   3450.2]]
T=298.15
for P, V, beta, rho in perplex_output:
    fo.set_state(P*1.e5,T)
    PT.append([P/1.e4,T])
    diff=[p(fo.V, V/1.e5), p(fo.K_T, 1.e5/beta), p(fo.density(), rho)]
    print diff
    percentage_diff.append(diff)

percentage_diff=np.array(percentage_diff)
i,j = np.unravel_index(percentage_diff.argmax(), percentage_diff.shape)

print 'Maximum error in', database, 'database:'
print variables[j], ':', percentage_diff[i,j], '% at', PT[i][0], 'GPa and', PT[i][1], 'K'
print ''
