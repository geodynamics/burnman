from __future__ import absolute_import
# Benchmarks for the chemical potential functions
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy import optimize

# Equilibrium functions


def eqm_P_xMgB(A, B):
    def eqm(arg, T, xMgA):
        P = arg[0]
        xMgB = arg[1]

        A.set_composition([xMgA, 1.0 - xMgA])
        A.set_state(P, T)

        B.set_composition([xMgB, 1.0 - xMgB])
        B.set_state(P, T)

        diff_mu_Mg2SiO4 = A.partial_gibbs[0] - B.partial_gibbs[0]
        diff_mu_Fe2SiO4 = A.partial_gibbs[1] - B.partial_gibbs[1]
        return [diff_mu_Mg2SiO4, diff_mu_Fe2SiO4]
    return eqm


def eqm_P_xMgABC(A, B, C):
    def eqm(arg, T):
        P = arg[0]
        xMgA = arg[1]
        xMgB = arg[2]
        xMgC = arg[3]

        A.set_composition([xMgA, 1.0 - xMgA])
        A.set_state(P, T)

        B.set_composition([xMgB, 1.0 - xMgB])
        B.set_state(P, T)

        C.set_composition([xMgC, 1.0 - xMgC])
        C.set_state(P, T)

        diff_mu_Mg2SiO4_0 = A.partial_gibbs[0] - B.partial_gibbs[0]
        diff_mu_Fe2SiO4_0 = A.partial_gibbs[1] - B.partial_gibbs[1]
        diff_mu_Mg2SiO4_1 = A.partial_gibbs[0] - C.partial_gibbs[0]
        diff_mu_Fe2SiO4_1 = A.partial_gibbs[1] - C.partial_gibbs[1]

        return [diff_mu_Mg2SiO4_0, diff_mu_Fe2SiO4_0, diff_mu_Mg2SiO4_1, diff_mu_Fe2SiO4_1]
    return eqm


'''
Initialise solid solutions
'''
ol = burnman.minerals.SLB_2011.mg_fe_olivine()
wd = burnman.minerals.SLB_2011.mg_fe_wadsleyite()
rw = burnman.minerals.SLB_2011.mg_fe_ringwoodite()

'''
Temperature of phase diagram
'''
T = 1673.  # K

'''
Find invariant point
'''
invariant = optimize.fsolve(
    eqm_P_xMgABC(ol, wd, rw), [15.e9, 0.2, 0.3, 0.4], args=(T))
print(str(invariant[0] / 1.e9) + ' GPa')
print(invariant[1:4])

'''
Initialise arrays
'''
XMgA_ol_wad = np.linspace(invariant[1], 0.9999, 21)
XMgA_ol_rw = np.linspace(0.0001, invariant[1], 21)
XMgA_wad_rw = np.linspace(invariant[2], 0.9999, 21)

P_ol_wad = np.empty_like(XMgA_ol_wad)
XMgB_ol_wad = np.empty_like(XMgA_ol_wad)

P_ol_rw = np.empty_like(XMgA_ol_wad)
XMgB_ol_rw = np.empty_like(XMgA_ol_wad)

P_wad_rw = np.empty_like(XMgA_ol_wad)
XMgB_wad_rw = np.empty_like(XMgA_ol_wad)

'''
Find transition pressures
'''

for idx, XMgA in enumerate(XMgA_ol_wad):
    XMgB_guess = 1.0 - ((1.0 - XMgA_ol_wad[idx]) * 0.8)
    P_ol_wad[idx], XMgB_ol_wad[idx] = optimize.fsolve(
        eqm_P_xMgB(ol, wd), [5.e9, XMgB_guess], args=(T, XMgA_ol_wad[idx]))
    XMgB_guess = 1.0 - ((1.0 - XMgA_ol_rw[idx]) * 0.8)
    P_ol_rw[idx], XMgB_ol_rw[idx] = optimize.fsolve(
        eqm_P_xMgB(ol, rw), [5.e9, XMgB_guess], args=(T, XMgA_ol_rw[idx]))
    XMgB_guess = 1.0 - ((1.0 - XMgA_wad_rw[idx]) * 0.8)
    P_wad_rw[idx], XMgB_wad_rw[idx] = optimize.fsolve(
        eqm_P_xMgB(wd, rw), [5.e9, XMgB_guess], args=(T, XMgA_wad_rw[idx]))

'''
Plot model
'''
fig1 = mpimg.imread('../../burnman/data/input_figures/slb_fig10a.png')
plt.imshow(fig1, extent=[0, 1, 0., 30.], aspect='auto')

plt.plot(1.0 - np.array([invariant[1], invariant[2], invariant[3]]), np.array(
    [invariant[0], invariant[0], invariant[0]]) / 1.e9, color='black', linewidth=2, label='invariant')

plt.plot(1.0 - XMgA_ol_wad, P_ol_wad / 1.e9,
         'r-', linewidth=2, label='wad-out (ol, wad)')
plt.plot(1.0 - XMgB_ol_wad, P_ol_wad / 1.e9,
         'g-', linewidth=2, label='ol-out (ol, wad)')

plt.plot(1.0 - XMgA_ol_rw, P_ol_rw / 1.e9,
         'r-',  linewidth=2, label='rw-out (ol, rw)')
plt.plot(1.0 - XMgB_ol_rw, P_ol_rw / 1.e9, 'b-',
         linewidth=2, label='ol-out (ol, rw)')

plt.plot(1.0 - XMgA_wad_rw, P_wad_rw / 1.e9,
         'g-',  linewidth=2, label='rw-out (wad, rw)')
plt.plot(1.0 - XMgB_wad_rw, P_wad_rw / 1.e9,
         'b-',  linewidth=2, label='wad-out (wad, rw)')


plt.title('Mg2SiO4-Fe2SiO4 phase diagram')
plt.xlabel("X_Fe")
plt.ylabel("Pressure (GPa)")
plt.legend(loc='upper right')
plt.show()
