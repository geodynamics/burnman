# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, 2014, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import constants

T_0=298.15 # Standard temperature = 25 C
P_0=1.e5 # Standard pressure = 1.e5 Pa

"""
Functions for order-disorder thermodynamic properties of endmembers
"""

def landau_ordering(P, T, params):
    """
    Returns critical temperature [K] at pressure [Pa], 
    and the states of order [unitless] at the 
    reference temperature and temperature of interest  

    from Holland and Powell, 1996 
    """
    Tcstar=params['landau_Tc'] + (params['landau_Vmax']/params['landau_Smax'])*P
            # Q_0 is Q at T0, P0? 
    Q_0=np.power((params['landau_Tc']-T_0)/params['landau_Tc'],1./4.)
    
    # Find state of ordering
    # Note that Q > 1 where Vmax*P > Smax*T. 
    if Tcstar-T > 0.:
        Q=np.power((Tcstar-T)/params['landau_Tc'],1./4.)
    else:
        Q=0.0

    return Tcstar, Q_0, Q

def gibbs_disorder_Landau(P, T, params):
    """
    Returns the Gibbs energy of disordering [J/mol] relative 
    to the equilibrium state of disorder (Q [unitless]) 
    based on the landau model

    from Holland and Powell, 1996, corrected using
    landaunote.pdf on Tim Holland's esc web page
    """
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
        
    # Vt is defined as a "function shaping the excess volume of disordering" in landaunote.pdf.
    # A sensible function to use in this case would be the expression in brackets in  EQ13 of Holland and Powell (2011) 
    # Here we follow tc337L, where Vt=1 ("the simplest way to proceed"; Roger Powell, pers. comm. 2014/08/22). 
    Vt=1.# EQ 13 would give: 1. - a + (a*(pow((1.-b*Pth), 1.-c) - pow((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure))
    Gdisord=params['landau_Tc']*params['landau_Smax']*(np.power(Q_0,2) - np.power(Q_0,6)/3.0) - params['landau_Smax']*(Tcstar*np.power(Q,2) - params['landau_Tc']*np.power(Q,6)/3.0) - T*(params['landau_Smax']*(np.power(Q_0,2) - np.power(Q,2))) + P*(params['landau_Vmax']*np.power(Q_0,2)*Vt)

    return Gdisord

def entropy_disorder_Landau(P, T, params):
    """
    Returns the entropy of disordering [J/K/mol] relative 
    to the equilibrium state of disorder (Q [unitless]) 
    based on the landau model

    from Holland and Powell, 1996
    """
    # N.B. Assumes Vt==1, see above
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
    Sdisord=params['landau_Smax']*(np.power(Q_0, 2.) - np.power(Q, 2.))
    return Sdisord

def enthalpy_disorder_Landau(P, T, params):
    """
    Returns the enthalpy of disordering [J/mol] relative 
    to the equilibrium state of disorder (Q [unitless]) 
    based on the landau model

    from Holland and Powell, 1996, corrected using
    landaunote.pdf on Tim Holland's esc web page
    """
    # N.B. Assumes Vt==1, see above
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
    Hdisord=gibbs_disorder_Landau(P, T, params) + T*entropy_disorder_Landau(P, T, params)
    return Hdisord

def heat_capacity_p_disorder_Landau(P, T, params):
    """
    Returns the heat capacity of disordering [J/mol] relative 
    to the equilibrium state of disorder (Q [unitless]) 
    based on the landau model

    from Holland and Powell, 1996, corrected using
    landaunote.pdf on Tim Holland's esc web page
    """
    # N.B. Assumes Vt==1, see above
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
    Hdisord=gibbs_disorder_Landau(P, T, params) + T*entropy_disorder_Landau(P, T, params)
    if Q==0.:
        Cp_disord=0.
    else:

        Cp_disord=params['landau_Smax']/(2.*np.sqrt((Tcstar-T)*params['landau_Tc']))
    return Cp_disord


# see derivation in thermodynamic_introduction.pdf
def lnxord(n,Q):
    return np.log(1.+n*Q) + n*np.log(n+Q) - (1.+n)*np.log(1.+n)

# see derivation in thermodynamic_introduction.pdf
def lnxdisord(n,Q):
    return (1./(1.+n))*np.log(1.+n*Q) + (n/(1.+n))*np.log(1.-Q) + (n/(1.+n))*np.log(n*(1.-Q)) + (n*n/(1.+n))*np.log(n+Q) - n*np.log(n)

# see derivation in thermodynamic_introduction.pdf
def entropydisorder(n):
    BW_deltaS=R*((1.+n)*np.log(1.+n) - n*np.log(n))
    return BW_deltaS

# see Holland and Powell, 1996
# params['BW_factor'] is a proportional reduction of the configurational energy of mixing, so feeds into the calculation of Q.
def equilibrium_Q(Q, deltaS, P, T, params):
    n=params['BW_n']
    W=params['BW_W'] + P*params['BW_Wv']
    if Q>1.0:
        Q=0.9 # A simple catch to make sure the optimisation doesn't fail
    return params['BW_deltaH'] - params['BW_factor']*T*deltaS + P*params['BW_deltaV'] + params['BW_factor']*constants.gas_constant*T*(lnxdisord(n,Q) - lnxord(n,Q)) + (2.*Q - 1.)*W

# Energy of disordering from Bragg-Williams symmetric model; see Holland and Powell, 1996
def gibbs_disorder_BW(P, T, params):
    n=params['BW_n']
    deltaS=entropydisorder(n)
    Q=opt.fsolve(equilibrium_Q, 0.999995, args=(deltaS, P, T, params))[0]
    W=params['BW_W'] + P*params['BW_Wv']
    ideal=(1.-Q)*(params['BW_deltaH'] - params['BW_factor']*T*entropydisorder(n) + P*params['BW_deltaV'] + params['BW_factor']*constants.gas_constant*T*lnxdisord(n,Q)) + params['BW_factor']*Q*(R*T*lnxord(n,Q))
    nonideal=(1.-Q)*Q*W
    Edisord=ideal+nonideal
    return Edisord

def entropy_disorder_BW(P, T, params):
    n=params['BW_n']
    deltaS=entropydisorder(n)
    Q=opt.fsolve(equilibrium_Q, 0.999995, args=(deltaS, P, T, params))[0]
    Sdisord=params['BW_factor']*((1.-Q)*(entropydisorder(n) - R*lnxdisord(n,Q)) - Q*(R*lnxord(n,Q)))
    return Sdisord
    
def enthalpy_disorder_BW(P, T, params):
    n=params['BW_n']
    deltaS=entropydisorder(n)
    Q=opt.fsolve(equilibrium_Q, 0.999995, args=(deltaS, P, T, params))[0]
    W=params['BW_W'] + P*params['BW_Wv']
    ideal=(1.-Q)*(params['BW_deltaH'] + P*params['BW_deltaV'])
    nonideal=(1.-Q)*Q*W
    Hdisord=ideal+nonideal
    return Hdisord
