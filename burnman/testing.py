# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.



# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import os, sys, numpy as np, matplotlib.pyplot as plt

#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))


# Here we import the relevant modules from BurnMan.  The burnman
# module imports several of the most important functionalities of
# the library, including the ability to make composites, and compute
# thermoelastic properties of them.  The minerals module includes
# the mineral physical parameters for the predefined minerals in
# BurnMan
import burnman
from burnman import minerals


# The following script checks the Gibbs Free energies of *solid* Holland and Powell (2011) endmembers satisfying one of several cases

print 'N.B. Small differences (<100 J/mol) between the burnman and thermocalc outputs are presumably the result of Holland and Powell missing .15 K somewhere in their code, either in the standard state temperature (298.15 K), or in the conversion from Celsius to Kelvin.'
print ''


# A mineral without phase transitions: stishovite
stv=minerals.HP_2011.stishovite()
stv.set_method("mtait")

Gstv=[['P \ T',25,500,1500],[0.001, -883.54, -909.13, -1020.78],[10.0, -869.55, -895.00, -1006.26],[100.0, -745.58, -769.83, -877.86]]

print 'G(stv)'

Pmaxerror=0.
Tmaxerror=0.
maxerror=0.
serror=0.
for pi in range(len(Gstv)-1):
    P=Gstv[pi+1][0]*1e8
    for ti in range(len(Gstv[0])-1):
        T=Gstv[0][ti+1]+273.15
        Gburnman=stv.calcgibbs(P,T)
        GHP=Gstv[pi+1][ti+1]*1000. # convert to J/mol
        Gdiff=Gburnman-GHP
#        print P, T, Gburnman, GHP, Gdiff
        serror=serror + pow(Gdiff,2)
        if abs(Gdiff) > abs(maxerror):
            maxerror = Gdiff
            Tmaxerror = T
            Pmaxerror = P

rmserror = np.sqrt(serror/((len(Gstv)-1)*(len(Gstv[0])-1)))
print 'Root mean squared error:', rmserror, 'J/mol'
print 'Maximum error:', maxerror, 'J/mol @', P/1.e9, 'GPa and', T, 'K'
print ''

# A mineral with a phase transition described by Landau theory: quartz
q=minerals.HP_2011.quartz()
q.set_method("mtait")

Gq=[['P \ T',25,500,1500],[0.001,-923.06,-957.10,-1088.42],[10.0,-900.62,-934.16,-1064.93],[100.0,-715.40,-746.62,-870.48]]

print 'G(q)'

Pmaxerror=0.
Tmaxerror=0.
maxerror=0.
serror=0.
for pi in range(len(Gq)-1):
    P=Gq[pi+1][0]*1e8
    for ti in range(len(Gq[0])-1):
        T=Gq[0][ti+1]+273.15
        Gburnman=q.calcgibbs(P,T)
        GHP=Gq[pi+1][ti+1]*1000. # convert to J/mol
        Gdiff=Gburnman-GHP
#        print P, T, Gburnman, GHP, Gdiff
        serror=serror + pow(Gdiff,2)
        if abs(Gdiff) > abs(maxerror):
            maxerror = Gdiff
            Tmaxerror = T
            Pmaxerror = P

rmserror = np.sqrt(serror/((len(Gq)-1)*(len(Gq[0])-1)))
print 'Root mean squared error:', rmserror, 'J/mol'
print 'Maximum error:', maxerror, 'J/mol @', P/1.e9, 'GPa and', T, 'K'
print ''


# A mineral with a phase transition described by Bragg-Williams theory: spinel
sp=minerals.HP_2011.spinel()
sp.set_method("mtait")

Gsp=[['P \ T',25,500,1500],[0.001,-2325.64,-2401.93,-2713.86],[10.0,-2285.96,-2361.80,-2672.59],[100.0,-1937.39,-2009.65,-2311.34]]

print 'G(sp)'

Pmaxerror=0.
Tmaxerror=0.
maxerror=0.
serror=0.
for pi in range(len(Gsp)-1):
    P=Gsp[pi+1][0]*1e8
    for ti in range(len(Gsp[0])-1):
        T=Gsp[0][ti+1]+273.15
        Gburnman=sp.calcgibbs(P,T)
        GHP=Gsp[pi+1][ti+1]*1000. # convert to J/mol
        Gdiff=Gburnman-GHP
        print P, T, Gburnman, GHP, Gdiff
        serror=serror + pow(Gdiff,2)
        if abs(Gdiff) > abs(maxerror):
            maxerror = Gdiff
            Tmaxerror = T
            Pmaxerror = P

rmserror = np.sqrt(serror/((len(Gsp)-1)*(len(Gsp[0])-1)))
print 'Root mean squared error:', rmserror, 'J/mol'
print 'Maximum error:', maxerror, 'J/mol @', P/1.e9, 'GPa and', T, 'K'
print ''
