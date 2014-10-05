import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman.minerals.SLB_2011 import *
from burnman.gibbsminimization import *
import numpy as np
from scipy.linalg import lstsq
import scipy.optimize as opt
import matplotlib.pyplot as plt
from burnman.equilibriumequations import *
composition = { 'Mg': 1.8, 'Fe': 0.2, 'O': 4., 'Si': 1.}

mineral_list = [mg_fe_olivine(), mg_fe_wadsleyite()]
for mineral in mineral_list:
    mineral.set_method('slb3')
    if isinstance(mineral, burnman.SolidSolution):
        molar_fractions=np.zeros(len(mineral.base_material))
        molar_fractions[0]=1.0
        mineral.set_composition(molar_fractions) 
    mineral.set_state(0.,0.)

'''
molar_fractions=np.array([0.9,0.1])
ol=mg_fe_olivine()
fo=forsterite()
fa=fayalite()

ol.set_composition(molar_fractions)
ol.set_state(P, T)
fo.set_method('slb3')
fa.set_method('slb3')
fo.set_state(P,T)
fa.set_state(P,T)

# Test excess function
print ol.gibbs
print molar_fractions[0]*fo.gibbs + molar_fractions[1]*fa.gibbs + ol.excess_gibbs

# Test partial gibbs functions
print ol.excess_gibbs
print 0.9*ol.calcpartialgibbsexcesses(P,T,molar_fractions)[0] + 0.1*ol.calcpartialgibbsexcesses(P,T,molar_fractions)[1]
'''

stoic=assemble_stoichiometric_matrix ( mineral_list )
null=sparsify_basis(compute_nullspace(stoic[0]))

# Create a compositional vector with elements in the same order as the 
# stoichiometric matrix
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
gvars=['P (Pa)', 'T (K)']
P=1e9
T=2000.
guesses=[P, T] 

for mineral in mineral_list:
    gvars.extend(['p(' + mineral.name + ')'])
    guesses.extend([1.0])
    if isinstance(mineral, burnman.SolidSolution):
        for i in range(len(mineral.base_material)-1):
            gvars.extend(['x(' + str(mineral.base_material[0][0].params['name']) + ')'])
            if i==0:
                guesses.extend([1.0])
            else:
                guesses.extend([0.0])

print gvars
print guesses


# Fix T, P, proportions or mixing parameters
fvars=[]
fvars.append([1, 1000., 2000., 11])
fvars.append([4, 0., 0., 1.]) # For some reason the solver fails when I set "2, 0., 0., 1", i.e. p(olivine)=0. This seems to be order dependent ... the solver works fine when I change the order of olivine and wadsleyite... :/ 

'''
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
'''

print gvars[fvars[0][0]], fvars[0][1], '<= x <=', fvars[0][2], 'nsteps=', fvars[0][3]
print gvars[fvars[1][0]], fvars[1][1], '<= x <=', fvars[1][2], 'nsteps=', fvars[1][3]

# Set up problem and solve it
print "\n", ' '.join(x.rjust(10) for x in gvars)
fixed_vars=[]
for i in np.linspace(fvars[0][1],fvars[0][2],fvars[0][3]):
    guesses[fvars[0][0]]=i
    if fvars[0][0] == 2:
        guesses[3]=1.0
    for j in np.linspace(fvars[1][1],fvars[1][2],fvars[1][3]):
        guesses[fvars[1][0]]=j
        fixed_vars=[[fvars[0][0], i],[fvars[1][0], j]]
        if fvars[1][0] == 2:
            guesses[3]=1.0
        elif fvars[0][0] == 2 and fvars[1][0] == 3:
            guesses[4]=1.0
        soln=opt.fsolve(set_eqns,guesses,args=(comp_vector, mineral_list, stoic, null, fixed_vars), full_output=1, xtol=1e-10)
        if soln[2]==1:
            print " ".join(str(("%10.5f" % x)) for x in soln[0])
            guesses=soln[0]
        else:
            print " ".join(str(("%10.5f" % x)) for x in soln[0]), soln[3]



# Plot up the last iteration

P=soln[0][0]
T=soln[0][1]

ol=mg_fe_olivine()
wad=mg_fe_wadsleyite()

comp = np.linspace(0, 0.2, 100)
ol_gibbs = np.empty_like(comp)
wad_gibbs = np.empty_like(comp)


for i,c in enumerate(comp):
   ol.set_composition( np.array([1.0-c, c]) )
   wad.set_composition( np.array([1.0-c, c]) )
   ol.set_state( P, T )
   wad.set_state( P, T )
   ol_gibbs[i] = ol.gibbs
   wad_gibbs[i] = wad.gibbs


# Detrend
a=(ol_gibbs[len(comp)-1] - ol_gibbs[0])/(comp[len(comp)-1] - comp[0])
b=ol_gibbs[0] - a*comp[0]

for i,c in enumerate(comp):
    ol_gibbs[i] = ol_gibbs[i] - (a*comp[i] + b)
    wad_gibbs[i] = wad_gibbs[i] - (a*comp[i] + b)


comp_tangent=np.array([1.-soln[0][3], 1.-soln[0][5]])
ol.set_composition( np.array([1.-comp_tangent[0], comp_tangent[0]]) )
wad.set_composition( np.array([1.-comp_tangent[1], comp_tangent[1]]) )
ol.set_state( P, T )
wad.set_state( P, T )

gibbs_tangent=np.array([ol.gibbs - (a*comp_tangent[0] + b), wad.gibbs - (a*comp_tangent[1] + b)])


plt.plot( comp, ol_gibbs, 'b', linewidth=1., label='ol')
plt.plot( comp, wad_gibbs, 'g', linewidth=1., label='wad')
plt.plot( comp_tangent, gibbs_tangent, 'r', linewidth=1., label='ol-wad equilibrium')
plt.title('olivine-wadsleyite equilibrium at '+str(round(P/1.e9,2))+' GPa and '+str(round(T,0))+' K')
plt.xlim(0.0,0.2)
plt.ylabel("Detrended Gibbs free energy (kJ/mol)")
plt.xlabel("x(ol), x(wad)")
plt.legend(loc='upper right')
plt.show()
