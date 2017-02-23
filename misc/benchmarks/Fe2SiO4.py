from __future__ import absolute_import
from __future__ import print_function
import os.path
import sys
sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman.minerals import DKS_2013_liquids
from burnman.minerals import RS_2014_liquids
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

fa = burnman.minerals.SLB_2011.fayalite()
fa_liq = RS_2014_liquids.Fe2SiO4_liquid()

test=True
if test==True:
    TE = np.loadtxt(fname='Ramo.dat')

    '''
    for temperature, rho, Ee, Ei, Ex in TE:
        volume = fa_liq.params['molar_mass']/rho

        E_el = (fa_liq.method._F_el(temperature, volume, fa_liq.params)
                + temperature
                * fa_liq.method._S_el(temperature, volume, fa_liq.params))
        E_ig = (fa_liq.method._F_ig(temperature, volume, fa_liq.params)
                + temperature
                * fa_liq.method._S_ig(temperature, volume, fa_liq.params))

        E_xs = (fa_liq.method._F_xs(temperature, volume, fa_liq.params)
                + temperature
                * fa_liq.method._S_xs(temperature, volume, fa_liq.params))

        print(E_el/Ee, E_ig/Ei, E_xs/Ex)
    '''


    Pm = 1.e5
    Tm = 1573.


    dT = 1.
    dV = 1.e-8


    fa.set_state(Pm, Tm)
    fa_liq.set_state(Pm, Tm)

    print(fa.V, fa_liq.V, fa_liq.rho)

    Vm = fa_liq.V

    temperatures = np.linspace(Tm - dT, Tm + dT, 3)
    volumes = np.linspace(Vm - dV, Vm + dV, 3)

    F = np.empty_like(temperatures)
    S = np.empty_like(temperatures)
    C_v = np.empty_like(temperatures)
    P = np.empty_like(temperatures)
    K_T = np.empty_like(temperatures)
    alphaK_T = np.empty_like(temperatures)

    for i, T in enumerate(temperatures):
        F[i] = fa_liq.method._F_mag(T, Vm, fa_liq.params)
        S[i] = fa_liq.method._S_mag(T, Vm, fa_liq.params)
        C_v[i] = fa_liq.method._C_v_mag(T, Vm, fa_liq.params)

    F_mag = F[1]
    S_mag_analytic = S[1]
    S_mag_numeric = -np.gradient(F)[1]/dT
    C_v_mag_analytic = C_v[1]
    C_v_mag_numeric = Tm*np.gradient(S)[1]/dT

    '''
        def _alphaK_T_mag(self, temperature, volume, params):
            S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
            S = S_a*temperature + S_b
            d2FdVdT = (-2.*params['spin_b'][1]*temperature/(params['V_0']*np.power(VoverVx, 2.))
                       + numerator)/(2.*S + 1.) - 2.*temperature*S_a*numerator/(np.power((2.*S+1.), 2.))
            return -n_atoms*constants.gas_constant*d2FdVdT
    '''

    for i, V in enumerate(volumes):
        F[i] = fa_liq.method._F_mag(Tm, V, fa_liq.params)
        P[i] = fa_liq.method._P_mag(Tm, V, fa_liq.params)
        K_T[i] = fa_liq.method._K_T_mag(Tm, V, fa_liq.params)

        alphaK_T[i] = fa_liq.method._alphaK_T_mag(Tm, V, fa_liq.params)
        S[i] = fa_liq.method._S_mag(Tm, V, fa_liq.params)

    P_mag_analytic = P[1]
    P_mag_numeric = -np.gradient(F)[1]/dV
    K_T_mag_analytic = K_T[1]
    K_T_mag_numeric = -Vm*np.gradient(P)[1]/dV
    alphaK_T_mag_analytic = alphaK_T[1]
    alphaK_T_mag_numeric = np.gradient(S)[1]/dV

    print('P:', P_mag_analytic/P_mag_numeric)
    print('K_T:', K_T_mag_analytic/K_T_mag_numeric)
    print('alphaK_T:', alphaK_T_mag_analytic/alphaK_T_mag_numeric)
    print('S:', S_mag_analytic/S_mag_numeric)
    print('C_v:', C_v_mag_analytic/C_v_mag_numeric)


    pressures = np.linspace(1.e5, 200.e9, 51)
    rhos = np.empty_like(pressures)


    fig1 = mpimg.imread('figures/Fe2SiO4_liquid_PVT.png')
    plt.imshow(fig1, extent=[3.5, 7.5, 0., 200], aspect='auto')

    for T in [3000., 4000., 6000.]: #, 5000., 2000., 1000.]:
        for i, P in enumerate(pressures):
            fa_liq.set_state(P, T)

            #print(P, T, fa_liq.rho)
            rhos[i] = fa_liq.rho
        plt.plot(rhos/1.e3, pressures/1.e9, label=str(T)+' K')

    plt.legend(loc='upper left')
    plt.show()

    fig1 = mpimg.imread('figures/Fe2SiO4_liquid_hugoniot.png')
    plt.imshow(fig1, extent=[3.5, 7.5, 0, 200], aspect='auto')


    temperatures, volumes = burnman.tools.hugoniot(fa_liq, 1.e5, 1573., pressures)
    plt.plot(fa_liq.molar_mass/volumes/1.e3, pressures/1.e9)


    plt.show()
