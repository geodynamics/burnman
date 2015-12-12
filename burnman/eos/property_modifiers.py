# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import
import numpy as np
from burnman.constants import gas_constant

"""
Functions for modifying the thermodynamic properties of minerals
Currently includes modifications for second order transitions
(landau, landau_hp), order-disorder (bragg_williams), 
magnetism (magnetic_chs), and a linear gibbs energy modification 
with respect to pressure and temperature (dqf). 
"""


def _landau_excesses(pressure, temperature, params):
    """
    Applies a tricritical Landau correction to the properties 
    of an endmember which undergoes a displacive phase transition. 
    This correction follows Putnis (1992), and is done relative to 
    the completely *ordered* state (at 0 K). 
    It therefore differs in implementation from both 
    Stixrude and Lithgow-Bertelloni (2011) and 
    Holland and Powell (2011), who compute properties relative to 
    the completely disordered state and standard states respectively.

    The current implementation is preferred, as the excess
    entropy (and heat capacity) terms are equal to zero at 0 K.

    N.B. The excesses are for a *completely relaxed* mineral;
    i.e. the seismic wave propagation is *slow* compared to the 
    rate of reaction. 
    """

    Tc = params['Tc_0'] + params['V_D']*pressure/params['S_D']

    G_disordered = params['S_D']*((temperature - Tc) + params['Tc_0']/3.)
    dGdT_disordered = params['S_D']
    dGdP_disordered = -params['V_D']

    if temperature < Tc:
        # Wolfram input to check partial differentials
        # x = T, y = P, a = S, c = Tc0, d = V
        # D[D[a ((x - c - d*y/a)*(1 - x/(c + d*y/a))^0.5 + c/3*(1 - x/(c + d*y/a))^1.5), x], x]
        Q2 = np.sqrt(1. - temperature/Tc)
        G = params['S_D']*((temperature - Tc)*Q2 + params['Tc_0']*Q2*Q2*Q2/3.) - G_disordered
        dGdP = -params['V_D']*Q2*(1. + 0.5*temperature/Tc*(1. - params['Tc_0']/Tc)) - dGdP_disordered
        dGdT = params['S_D']*Q2*(1.5 - 0.5*params['Tc_0']/Tc) - dGdT_disordered
        d2GdP2 = params['V_D']*params['V_D']*temperature/(params['S_D']*Tc*Tc*Q2) \
            * (temperature*(1. + params['Tc_0']/Tc)/(4.*Tc) 
               + Q2*Q2*(1. - params['Tc_0']/Tc) - 1.)
        d2GdT2 = -params['S_D']/(Tc*Q2)*(0.75 - 0.25*params['Tc_0']/Tc)
        d2GdPdT = params['V_D']/(2.*Tc*Q2)*(1. + (temperature / (2.*Tc) - Q2*Q2)
                                            *(1. - params['Tc_0']/Tc))
                                         
    else:
        Q = 0.
        G = -G_disordered
        dGdT = -dGdT_disordered
        dGdP = -dGdP_disordered
        d2GdT2 = 0.
        d2GdP2 = 0.
        d2GdPdT = 0.

        
    excesses = {'G': G, 'dGdT': dGdT, 'dGdP': dGdP,
                'd2GdT2': d2GdT2, 'd2GdP2': d2GdP2, 'd2GdPdT': d2GdPdT}
    
    return excesses

def _landau_hp_excesses(pressure, temperature, params):
    """
    Applies a tricritical Landau correction to the properties 
    of an endmember which undergoes a displacive phase transition. 
    This correction is done relative to the standard state, as per
    Holland and Powell (1998).

    Includes the correction published within landaunote.pdf
    (Holland, pers. comm), which 'corrects' the terms involving 
    the critical temperature Tc / Tc*

    Note that this formalism is still inconsistent, as it predicts that 
    the order parameter can be greater than one. For this reason
    _landau_excesses is preferred.

    N.B. The excesses are for a *completely relaxed* mineral;
    i.e. the seismic wave propagation is *slow* compared to the 
    rate of reaction. 
    """
    params = mineral.landau_HP
    P = mineral.pressure
    T = mineral.temperature
    P_0 = mineral.params['P_0']
    T_0 = mineral.params['T_0']

    if T_0 < params['Tc_0']:
        Q_0 = np.power((params['Tc_0'] - T_0)/params['Tc_0'], 0.25)
    else:
        Q_0 = 0.

    Tc = params['Tc_0'] + params['V_D']*(P-P_0)/params['S_D']
    if T < Tc:
        Q = np.power((Tc - T)/params['Tc_0'], 0.25)
    else:
        Q = 0.

    # Gibbs
    G = params['Tc_0']*params['S_D']*(Q_0*Q_0 - np.power(Q_0, 6.)/3.) \
        - params['S_D']*(Tc*Q*Q - params['Tc_0']*np.power(Q, 6.)/3.) \
        - T*params['S_D']*(Q_0*Q_0 - Q*Q) + (P-P_0)*params['V_D']*Q_0*Q_0
    
    
    dGdT = params['S_D']*(Q*Q - Q_0*Q_0)
    dGdP = -params['V_D']*(Q*Q - Q_0*Q_0)

    if Q > 1.e-12:
        d2GdT2 = -params['S_D']/(2.*params['Tc_0']*Q*Q)
        d2GdP2 = -params['V_D']*params['V_D']/(2.*params['S_D']*params['Tc_0']*Q*Q)
        d2GdPdT = -params['V_D']/(2.*params['Tc_0']*Q*Q)
    else:
        d2GdT2 = 0.
        d2GdP2 = 0.
        d2GdPdT = 0.

    excesses = {'G': G, 'dGdT': dGdT, 'dGdP': dGdP,
                'd2GdT2': d2GdT2, 'd2GdP2': d2GdP2, 'd2GdPdT': d2GdPdT}
    
    return excesses

def _dqf_excesses(pressure, temperature, params):
    """
    Applies a 'Darken's quadratic formalism' correction (Powell, 1987)
    to the thermodynamic properties of a mineral endmember.
    This correction is relative to P = 0 and T = 0 and linear in P and T
    and therefore corresponds to a constant volume and entropy correction.

    Applying either a volume or entropy term will generally break
    equations of state (i.e. the properties of the mineral will 
    no longer obey the equation of state defined in the 
    params dictionary. However, this form of excess is extremely 
    useful as a first order tweak to free energies 
    (especially in solid solution calculations)
    """

    G = mineral.dqf['H'] \
        - (mineral.temperature)*mineral.dqf['S'] \
        + (mineral.pressure)*mineral.dqf['V']
    dGdT = -mineral.dqf['S']
    dGdP = mineral.dqf['V']
    d2GdT2 = 0.
    d2GdP2 = 0.
    d2GdPdT = 0.

    excesses = {'G': G, 'dGdT': dGdT, 'dGdP': dGdP,
                'd2GdT2': d2GdT2, 'd2GdP2': d2GdP2, 'd2GdPdT': d2GdPdT}
    
    return excesses

def _bragg_williams_excesses(pressure, temperature, params):
    """
    Applies a Bragg-Williams type correction to the thermodynamic
    properties of a mineral endmember. Used for modelling 
    order-disorder processes.
    Expressions are from Holland and Powell (1996).

    N.B. The excesses are for a *completely relaxed* mineral;
    i.e. the seismic wave propagation is *slow* compared to the 
    rate of reaction. 

    This may not be reasonable for order-disorder, especially 
    for slow or coupled diffusers (Si-Al, for example).
    The completely *unrelaxed* mineral (in terms of order-disorder)
    can be calculated with a solid solution model.
    """

    R = constants.gas_constant
    n=params['n']
    f=params['factor']
    deltaS = entropydisorder(n)

    lnxord = lambda n, Q: np.log(1.+n*Q) + n*np.log(n+Q) - (1.+n)*np.log(1.+n)
    lnxdisord = lambda n, Q: (1./(1.+n))*np.log(1.+n*Q) + (n/(1.+n))*np.log(1.-Q) \
                + (n/(1.+n))*np.log(n*(1.-Q)) + (n*n/(1.+n))*np.log(n+Q) - n*np.log(n)
    
    def reaction_bragg_williams(Q, gibbs_disorder, temperature, n, f, W):
        if Q>1.0:
            Q=0.9 # A simple catch to make sure the optimisation doesn't fail
        return gibbs_disorder + (2.*Q - 1.)*W \
            f*R*temperature*(lnxdisord(n,Q) - lnxord(n,Q))

    def order_gibbs(pressure, temperature, params):
        W=params['Wh'] + pressure*params['Wv']
        gibbs_disorder = params['deltaH'] - f*temperature*deltaS + pressure*params['deltaV']
        Q = opt.fsolve(reaction_bragg_williams, 0.999995,
                       args=(gibbs_disorder, pressure, temperature, params))[0]
        G = (1.-Q) * (gibbs_disorder + f*R*temperature*lnxdisord(n,Q)) \
            + f*Q*(R*temperature*lnxord(n,Q)) + (1.-Q)*Q*W

        return Q, G

    # Calculating partial differentials with respect to P and T
    # are complicated by the fact that Q changes with P and T
    # Since there's no analytical solution for Q(P, T), we are
    # unfortunately driven to numerical differentiation. Schade.
    dT = 1.
    dP = 1000.
    
    Q, G = order_gibbs(pressure, temperature, params)
    Q, GsubPsubT = order_gibbs(pressure - dP, temperature - dT, params)
    Q, GsubPaddT = order_gibbs(pressure - dP, temperature + dT, params)
    Q, GaddPsubT = order_gibbs(pressure + dP, temperature - dT, params)
    Q, GaddPaddT = order_gibbs(pressure + dP, temperature + dT, params)
    Q, GsubP = order_gibbs(pressure - dP, temperature, params)
    Q, GaddP = order_gibbs(pressure + dP, temperature, params)
    Q, GsubT = order_gibbs(pressure, temperature - dT, params)
    Q, GaddT = order_gibbs(pressure, temperature + dT, params)

    dGdT = (GaddT - GsubT)/(2.*dT)
    dGdP = (GaddP - GsubP)/(2.*dP)
    d2GdT2 = (GaddT + GsubT - 2.*G) / (dT*dT)
    d2GdP2 = (GaddP + GsubP - 2.*G) / (dP*dP)
    d2GdPdT = (GaddPaddT - GsubPaddT - GaddPsubT + GsubPsubT)/(4.*dT*dP)
    
    excesses = {'G': G, 'dGdT': dGdT, 'dGdP': dGdP,
                'd2GdT2': d2GdT2, 'd2GdP2': d2GdP2, 'd2GdPdT': d2GdPdT}
    
    return excesses

def _magnetic_excesses_chs(pressure, temperature, params):
    """
    Applies a magnetic contribution to the thermodynamic 
    properties of a mineral endmember.
    The expression for the gibbs energy contribution is that 
    used by Chin, Hertzman and Sundman (1987) as reported 
    in the Journal of Phase Equilibria (Sundman, 1991).
    """
    
    structural_parameter=mineral.magnetic['magnetic_structural_parameter']
    tau=temperature/(mineral.magnetic['curie_temperature'][0] + pressure*mineral.magnetic['curie_temperature'][1])
    magnetic_moment=mineral.magnetic['magnetic_moment'][0] + pressure*mineral.magnetic['magnetic_moment'][1]

    A = (518./1125.) + (11692./15975.)*((1./structural_parameter) - 1.)
    if tau < 1: 
        f=1.-(1./A)*(79./(140.*structural_parameter*tau) + (474./497.)*(1./structural_parameter - 1.)*(np.power(tau, 3.)/6. + np.power(tau, 9.)/135. + np.power(tau, 15.)/600.))
    else:
        f=-(1./A)*(np.power(tau,-5)/10. + np.power(tau,-15)/315. + np.power(tau, -25)/1500.)


    G = gas_constant*temperature*np.log(magnetic_moment + 1.)*f
    dGdT = 0.
    dGdP = 0.
    d2GdT2 = 0.
    d2GdP2 = 0.
    d2GdPdT = 0.


    excesses = {'G': G, 'dGdT': dGdT, 'dGdP': dGdP,
                'd2GdT2': d2GdT2, 'd2GdP2': d2GdP2, 'd2GdPdT': d2GdPdT}
    
    return excesses


def _modify_properties(mineral, excesses):
    """
    Modifies the properties of a mineral based on
    excess contributions to thermodynamic properties
    (gibbs and the first and second derivatives with 
    respect to pressure and temperature) contained
    within the dictionary 'excesses'.
    """
    
    # Gibbs
    mineral.gibbs = mineral.gibbs + excesses['G']

    # Second derivatives first
    mineral.C_p = mineral.C_p - mineral.temperature*excesses['d2GdT2'] # -T*d2G/dT2
    mineral.K_T = - (mineral.V + excesses['dGdP']) / (excesses['d2GdP2'] - (mineral.V / mineral.K_T)) # - excesses['dGdP / (excesses['d2GdP2'])
    mineral.alpha = ((mineral.alpha*mineral.V) + excesses['d2GdPdT']) / (mineral.V + excesses['dGdP']) # d2GdPdT / dGdP

    # Now first derivatives 
    mineral.S = mineral.S - excesses['dGdT'] # dGdT
    mineral.V = mineral.V + excesses['dGdP'] # dGdP
    mineral.H = mineral.gibbs + mineral.temperature*mineral.S # H = G + TS

    # Now volume derivatives
    mineral.helmholtz = mineral.gibbs - mineral.pressure*mineral.V
    mineral.C_v = mineral.C_p - mineral.V*mineral.temperature*mineral.alpha*mineral.alpha*mineral.K_T
    mineral.gr = mineral.alpha*mineral.K_T*mineral.V/mineral.C_v
    mineral.K_S = mineral.K_T*mineral.C_p/mineral.C_v

    return None

def apply_property_modifiers(mineral):
    """
    Modifies the properties of a mineral according to 
    a list of modifiers.
    """
    for modifier in mineral.property_modifiers:
        if modifier[0] == 'landau':
            _modify_properties(mineral, _landau_excesses(mineral.pressure, mineral.temperature, modifier[1])
        if modifier[0] == 'landau_hp':
            _modify_properties(mineral, _landau_hp_excesses(mineral.pressure, mineral.temperature, modifier[1])
        if modifier[0] == 'dqf':
            _modify_properties(mineral, _dqf_excesses(mineral.pressure, mineral.temperature, modifier[1])
        if modifier[0] == 'bragg_williams':
            _modify_properties(mineral, _bragg_williams_excesses(mineral.pressure, mineral.temperature, modifier[1])
        if modifier[0] == 'magnetic_chs':
            _modify_properties(mineral, _magnetic_excesses_chs(mineral.pressure, mineral.temperature, modifier[1])
    
    return None



