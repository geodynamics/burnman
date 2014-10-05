import os, sys, numpy as np, matplotlib.pyplot as plt
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
import numpy as np


def set_eqns(arg, bulk_composition, mineral_list, stoic, null, fixed_vars): 
    P=float(arg[0])
    T=float(arg[1]) # arg[0] and arg[1] are P, T

    mbridx=0
    minidx=0
    # mineral proportions and solid solution compositions
    propscomps=np.array(arg[2:])
    partialgibbs=np.empty(len(propscomps))
    endmember_proportions=np.empty(len(propscomps))

    # The following for loop
    # finds the partial Gibbs free energies of the endmembers
    for mineral in mineral_list:
        mineral_proportion=propscomps[mbridx]
        if isinstance(mineral, burnman.SolidSolution):
            mbrprops=np.zeros(len(mineral.base_material))
            mbrprops[0:len(mineral.base_material)-1]=propscomps[mbridx+1:mbridx+len(mineral.base_material)]
            mbrprops[len(mineral.base_material)-1]=1.0-sum(mbrprops) 
            mineral.set_composition(mbrprops)
            partial_gibbs_excesses=mineral.calcpartialgibbsexcesses(P, T, mbrprops)
            for idx, endmember in enumerate(mineral.base_material):
                endmember_proportions[minidx]=mineral_proportion*mbrprops[idx]
                minidx=minidx+1
                partialgibbs[mbridx]=endmember[0].calcgibbs(P, T) + partial_gibbs_excesses[idx]
                mbridx=mbridx+1
        else:
            endmember_proportions[minidx]=mineral_proportion
            minidx=minidx+1
            partialgibbs[mbridx]=mineral.calcgibbs(P, T)
            mbridx=mbridx+1
            
    #### ADD EQUATIONS TO SOLVE ####
    eqns=[]

    # Add independent endmember reactions
    gibbs_inequalities=np.dot(partialgibbs,null)
    for ieq in gibbs_inequalities:
        eqns.append(ieq)  

    # Add fixed variable constraints
    eqns.append(arg[fixed_vars[0][0]]-fixed_vars[0][1])
    eqns.append(arg[fixed_vars[1][0]]-fixed_vars[1][1])

    # Add bulk compositional constraints
    mineral_bulk_composition=np.dot(endmember_proportions,stoic[0].T)

#    for i in range(len(mineral_bulk_composition)):
    for i in range(2): 
        eqns.append(mineral_bulk_composition[i] - bulk_composition[i])

    return eqns
