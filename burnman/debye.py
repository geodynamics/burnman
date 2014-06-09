# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.integrate as integrate
from numpy.polynomial.chebyshev import Chebyshev

"""
Functions for the Debye model.  Note that this is not Mie-Grueneisen-Debye,
just Debye, so is pretty limited.  Combine this with Mie-Grueneisen and
Birch-Murnaghan to get a full EOS
"""

R = 8.314462175



def debye_fn(x):
    """
    Evaluate the Debye function.  Takes the parameter
    xi = Debye_T/T
    """
    sol = integrate.quad( lambda xi: pow(xi,3.)/(np.exp(xi)-1.) , 0.0, x) # EQ B3
    return 3.*sol[0]/pow(x,3.)



chebyshev_representation = Chebyshev( [ 2.707737068327440945/2.0, 0.340068135211091751, -0.12945150184440869e-01, \
                                     0.7963755380173816e-03, -0.546360009590824e-04, 0.39243019598805e-05, \
                                    -0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08, \
                                     0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, -0.607797e-13, \
                                     0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17] )
eps = np.finfo(np.float).eps
sqrt_eps = np.sqrt(np.finfo(np.float).eps)
log_eps = np.log(np.finfo(np.float).eps)

def debye_fn_cheb(x):
    """
    Evaluate the Debye function using a Chebyshev series expansion coupled with
    asymptotic solutions of the function.  Shamelessly adapted from the GSL implementation
    of the same function (Itself adapted from Collected Algorithms from ACM).
    Should give the same result as debye_fn(x) to near machine-precision.
    """
    val_infinity = 19.4818182068004875;
    xcut = -log_eps

    assert(x > 0.0) #check for invalid x

    if x < 2.0*np.sqrt(2.0)*sqrt_eps:
        return 1.0 - 3.0*x/8.0 + x*x/20.0;
    elif x <= 4.0 :
        t = x*x/8.0 - 1.0;
        c = chebyshev_representation(t)
        return c - 0.375*x;
    elif x < -(np.log(2.0) + log_eps ):
        nexp = int(np.floor(xcut/x));
        ex  = np.exp(-x);
        xk  = nexp * x;
        rk  = nexp;
        sum = 0.0;
        for i in range(nexp,0,-1):
            xk_inv = 1.0/xk;
            sum *= ex;
            sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        return val_infinity/(x*x*x) - 3.0 * sum * ex;
    elif x < xcut:
        x3 = x*x*x;
        sum = 6.0 + 6.0*x + 3.0*x*x + x3;
        return  (val_infinity - 3.0 * sum * np.exp(-x)) / x3;
    else:
        return ((val_infinity/x)/x)/x;


def thermal_energy(T, debye_T, n):
    """
    calculate the thermal energy of a substance.  Takes the temperature,
    the Debye temperature, and n, the number of atoms per molecule.
    Returns thermal energy in J/mol
    """
    if T == 0:
        return 0
    E_th = 3.*n*R*T * debye_fn_cheb(debye_T/T)
    return E_th

def heat_capacity_v(T,debye_T,n):
    """
    Heat capacity at constant volume.  In J/K/mol
    """
    if T ==0:
        return 0
    x = debye_T/T
    C_v = 3.0*n*R* ( 4.0*debye_fn_cheb(x) - 3.0*x/(np.exp(x)-1.0) )
    return C_v




if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import time

    def old_thermal(T, debye_T, n):
        if T == 0:
            return 0
        return 3.*n*R*T * debye_fn(debye_T/T)

    def old_heat(T, debye_T, n):
        if T ==0:
            return 0
        deb = integrate.quad( lambda x : pow(x,4.)*np.exp(x)/pow((np.exp(x)-1.),2.), 0.0, debye_T/T)
        return 9.*n*R*deb[0]/pow(debye_T/T,3.)

    temperatures = np.linspace(100,5000, 10000)
    Debye_T = 1000.
    old = np.empty_like(temperatures)
    start = time.clock()
    for i in range(len(temperatures)):
        old[i] = old_heat(temperatures[i], Debye_T, 1.0)
    time_old = time.clock()-start

    new = np.empty_like(temperatures)
    start = time.clock()
    for i in range(len(temperatures)):
        new[i] = heat_capacity_v(temperatures[i], Debye_T, 1.0)
    time_new = time.clock()-start

    print "error %e"%np.linalg.norm((old-new)/new)
    print "time old %g, time new %g"%(time_old,time_new)



    temperatures = np.linspace(0,5000, 200)
    vibrational_energy = np.empty_like(temperatures)
    heat_capacity = np.empty_like(temperatures)
    Debye_T = 1000.
    for i in range(len(temperatures)):
      vibrational_energy[i] = thermal_energy(temperatures[i], Debye_T, 1.0)
      heat_capacity[i] = heat_capacity_v(temperatures[i], Debye_T, 1.0)

    plt.subplot(121)
    plt.plot(temperatures, vibrational_energy)
    plt.subplot(122)
    plt.plot(temperatures, heat_capacity)
    plt.show()



