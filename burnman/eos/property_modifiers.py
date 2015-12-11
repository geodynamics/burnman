import numpy as np
from burnman.constants import gas_constant

def _landau_excesses(pressure, temperature, params):
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
