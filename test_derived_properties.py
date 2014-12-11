import burnman
from burnman.minerals import SLB_2011
from burnman.minerals import HP_2011_ds62
import numpy as np
import matplotlib.pyplot as plt
###

f = open('perplex_output/fo_SLB2011_params.dat', 'r')
datalines = [ line.strip() for idx, line in enumerate(f.read().split('\n')) if line.strip() and idx>0 ]
f.close()
data = [ map(float,"%".join(line.split("%")[:1]).split()) for line in datalines ]
P, T, H, S, V, C_p, alpha, beta, rho = zip(*data)


fo = SLB_2011.fo()
fo.set_method('slb3')

print "P (GPa)  T(K)  Gibbs(burnman, J/mol)   Gibbs(SLB2011/PerpleX, J/mol)"
for i, pressure in enumerate(P):
    fo.set_state(P[i]*1.e5,T[i])
    print P[i]/1.e4, T[i], fo.gibbs, H[i]-T[i]*S[i], fo.gibbs - (H[i]-T[i]*S[i])

    print P[i], T[i], fo.H, H[i], (fo.H-H[i])/fo.H*100. # enthalpy
    print P[i], T[i], fo.S, S[i], (fo.S-S[i])/fo.S*100. # entropy
    print P[i], T[i], fo.V, V[i]*1.e-5 , (fo.V-V[i]*1.e-5)/fo.V*100. # volume
    print P[i], T[i], fo.C_p, C_p[i], (fo.C_p-C_p[i])/fo.C_p*100. # heat capacity
    print P[i], T[i], fo.alpha, alpha[i], (fo.alpha-alpha[i])/fo.alpha*100. # expansivity
    print P[i], T[i], fo.compressibility(), beta[i]*1.e-5, (fo.compressibility()-beta[i]*1.e-5)/fo.compressibility()*100. # compressibility
    print P[i], T[i], fo.density(), rho[i], (fo.density()-rho[i])/fo.density()*100. # density
    print ''
print ''


###

f = open('perplex_output/fo_HP2011_params.dat', 'r')
datalines = [ line.strip() for idx, line in enumerate(f.read().split('\n')) if line.strip() and idx>0 ]
f.close()
data = [ map(float,"%".join(line.split("%")[:1]).split()) for line in datalines ]

P, T, H, S, V, C_p, alpha, beta, rho = zip(*data)


fo = HP_2011_ds62.fo()
fo.set_method('mtait')

print "P (GPa)  T(K)  Gibbs(burnman, J/mol)   Gibbs(HP2011/PerpleX, J/mol)"
for i, pressure in enumerate(P):
    fo.set_state(P[i]*1.e5,T[i])
    print P[i]/1.e4, T[i], fo.gibbs, H[i]-T[i]*S[i], fo.gibbs - (H[i]-T[i]*S[i])

    print P[i], T[i], fo.H, H[i], (fo.H-H[i])/fo.H*100. # enthalpy
    print P[i], T[i], fo.S, S[i], (fo.S-S[i])/fo.S*100. # entropy
    print P[i], T[i], fo.V, V[i]*1.e-5 , (fo.V-V[i]*1.e-5)/fo.V*100. # volume
    print P[i], T[i], fo.C_p, C_p[i], (fo.C_p-C_p[i])/fo.C_p*100. # heat capacity
    print P[i], T[i], fo.alpha, alpha[i], (fo.alpha-alpha[i])/fo.alpha*100. # expansivity
    print P[i], T[i], fo.compressibility(), beta[i]*1.e-5, (fo.compressibility()-beta[i]*1.e-5)/fo.compressibility()*100. # compressibility
    print P[i], T[i], fo.density(), rho[i], (fo.density()-rho[i])/fo.density()*100. # density
    print ''


quartz = HP_2011_ds62.q()
quartz.set_method('mtait')
P=1.e6
T=1500.
quartz.set_state(P,T)
print quartz.S


temperatures=np.linspace(300., 1500., 1001.)
heat_capacities=np.empty_like(temperatures)
for i, T in enumerate(temperatures):
    quartz.set_state(P,T)
    heat_capacities[i]=quartz.C_p


plt.plot(temperatures, heat_capacities, 'g-')
plt.xlabel("Temperature (K)")
plt.ylabel("Cp quartz (J/K/mol)")
plt.title("Heat capacities of a mineral with a Landau transition")
plt.show()
