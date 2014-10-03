import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.minerals.SLB_2011 import *
from burnman.gibbsminimization import *
import numpy as np
from scipy.linalg import lstsq
import matplotlib.pyplot as plt

composition = { 'Mg': 1.8/7., 'Fe': 0.2/7., 'O': 4./7., 'Si': 1./7}

minlist = [mg_fe_olivine(), mg_fe_wadsleyite()]

stoic=assemble_stoichiometric_matrix ( minlist )
null=sparsify_basis(compute_nullspace(stoic[0]))

print null

# Create a compositional vector with elements in the same order as the 
# stoiciometric matrix
comp_vector=np.empty( (len(stoic[1])) )
for idx, element in enumerate(stoic[1]):
    comp_vector[idx]=composition[element]

# Check that the bulk composition can be described by a linear set of the
# endmembers
x,resid,rank,s = lstsq(stoic[0],comp_vector)
resid = np.dot(stoic[0],x) - comp_vector
tol=1e-10


if (any(abs(resid[i]) > tol for i in range(len(resid)))):
    print "Oh no! Bulk composition outside the compositional range of the chosen minerals ( Maximum residual:", max(resid),  ")"
    
    tol=1e-3
    if (any(abs(resid[i]) > tol for i in range(len(resid)))):
        print "Residuals too high to recast composition. Program exiting"
        sys.exit()
    else:
        print "Residuals low enough to recast composition. Doing this now..."
        comp_vector=np.dot(stoic[0],x)
else:
    print 'Composition good'




# Create list of unknowns and their starting guesses
gvars=['P (kbar)', 'T (K)']
guesses=[250, 1000, 1.0] # First *mineral* proportion starts at 1.0

gvars.extend(['p(' + mineral.name + ')' for mineral in minlist])



for mineral in minlist:
    if isinstance(mineral, burnman.SolidSolution):
        for i in range(len(mineral.base_material)-1):
            gvars.extend(['x(' + str(mineral.base_material[0][0].params['name']) + ')'])
            #guesses.extend(map(float,zip(*gueslist)[1]))

print gvars


# Fix T, P, proportions or mixing parameters
fvars=[]
print '\nSelect two variables to fix'
for i in range(len(gvars)):
    print i, gvars[i]

for i in range(2):
    varinp=raw_input("Variable index and value (or min,max,nsteps): ")
    var=varinp.split()
    if len(var) == 2:
        fvars.append([int(var[0]), float(var[1]), float(var[1]), 1])
    elif len(var) == 4:
        fvars.append([int(var[0]), float(var[1]), float(var[2]), int(var[3])])
    else:
        print 'Incorrect number of lines'
        sys.exit()

print gvars[fvars[0][0]], fvars[0][1], '<= x <=', fvars[0][2], 'nsteps=', fvars[0][3]
print gvars[fvars[1][0]], fvars[1][1], '<= x <=', fvars[1][2], 'nsteps=', fvars[1][3]

'''
# Set up problem and solve it with fsolve
print "\n", ' '.join(x.rjust(10) for x in gvars)
fixed_vars=[]
for i in linspace(fvars[0][1],fvars[0][2],fvars[0][3]):
    guesses[fvars[0][0]]=i
    if fvars[0][0] == 2:
        guesses[3]=1.0
    for j in linspace(fvars[1][1],fvars[1][2],fvars[1][3]):
        guesses[fvars[1][0]]=j
        fixed_vars=[[fvars[0][0], i],[fvars[1][0], j]]
        if fvars[1][0] == 2:
            guesses[3]=1.0
        elif fvars[0][0] == 2 and fvars[1][0] == 3:
            guesses[4]=1.0
        soln=fsolve(set_eqns,guesses,args=(comp, minlist, mbrlist, mbrcomp, mbr, model, null, lookup, fixed_vars), full_output=1, xtol=1e-10)
        if soln[2]==1:
            print " ".join(str(("%10.5f" % x)) for x in soln[0])
            guesses=soln[0]
        else:
            print " ".join(str(("%10.5f" % x)) for x in soln[0]), soln[3]

'''
