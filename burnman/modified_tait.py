# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


#TO DO: Correct heat capacity, volume where internal order-disorder is implemented (Landau and Bragg-Williams models)

import numpy as np
import scipy.optimize as opt

import burnman.equation_of_state as eos

T_0=298.15 # Standard temperature = 25 C
P_0=1.e5 # Standard pressure = 1.e5 Pa
R=8.3145 # J/K/mol

# see Holland and Powell, 2011
def einst(S, n):
    """
    Einstein temperature
    base of p.346, para.1
    """
    return 10636./(S/n + 6.44)

# see Holland and Powell, 2011
def ksi(u):
    """
    Einstein function to describe behaviour of ak
    EQ 11+1
    """
    return pow(u,2.)*np.exp(u)/pow((np.exp(u)-1.), 2.)

# see Holland and Powell, 2011
def tait_constants(params):
    a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
    b=params['Kprime_0']/params['K_0'] - params['Kdprime_0']/(1. + params['Kprime_0'])
    c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])
    return a, b, c

def gibbs_disord_Landau(P, T, params):
    # From Holland and Powell, 1996, corrected using
    # landaunote.pdf on Tim Holland's esc web page
    Tcstar=params['landau_Tc'] + (params['landau_Vmax']/params['landau_Smax'])*P
            # Q_0 is Q at T0, P0? 
    Q_0=pow((params['landau_Tc']-T_0)/params['landau_Tc'],1./4.)
    
    # Find state of ordering
    # Note that Q > 1 where Vmax*P > Smax*T. 
    if Tcstar-T > 0.:
        Q=pow((Tcstar-T)/params['landau_Tc'],1./4.)
    else:
        Q=0.0
        
    # Vt is defined as a "function shaping the excess volume of disordering" in landaunote.pdf.
    # A sensible function to use in this case would be the expression in brackets in  EQ13 of Holland and Powell (2011) 
    # Here we follow tc337L, where Vt=1 ("the simplest way to proceed"; Roger Powell, pers. comm. 2014/08/22). 
    Vt=1.# EQ 13 would give: 1. - a + (a*(pow((1.-b*Pth), 1.-c) - pow((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure))
    Gdisord=params['landau_Tc']*params['landau_Smax']*(pow(Q_0,2) - pow(Q_0,6)/3.0) - params['landau_Smax']*(Tcstar*pow(Q,2) - params['landau_Tc']*pow(Q,6)/3.0) - T*(params['landau_Smax']*(pow(Q_0,2) - pow(Q,2))) + P*(params['landau_Vmax']*pow(Q_0,2)*Vt)
    return Gdisord

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
    return params['BW_deltaH'] - params['BW_factor']*T*deltaS + P*params['BW_deltaV'] + params['BW_factor']*R*T*(lnxdisord(n,Q) - lnxord(n,Q)) + (2.*Q - 1.)*W

# Energy of disordering from Bragg-Williams symmetric model; see Holland and Powell, 1996
def gibbs_disord_BW(P, T, params):
    n=params['BW_n']
    deltaS=entropydisorder(n)
    Q=opt.fsolve(equilibrium_Q, 0.999995, args=(deltaS, P, T, params))[0]
    W=params['BW_W'] + P*params['BW_Wv']
    ideal=(1.-Q)*(params['BW_deltaH'] - params['BW_factor']*T*entropydisorder(n) + P*params['BW_deltaV'] + params['BW_factor']*R*T*lnxdisord(n,Q)) + params['BW_factor']*Q*(R*T*lnxord(n,Q))
    nonideal=(1.-Q)*Q*W
    Edisord=ideal+nonideal
    return Edisord

# see Holland and Powell, 2011
def modified_tait(x, params):
    """
    equation for the modified Tait equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    EQ 2
    """
    a, b, c = tait_constants(params)

    return (pow((x + a - 1.) / a, -1./c) - 1.)/b


class MTaitBase(eos.EquationOfState):
    """
    Base class for a generic modified Tait equation of state.  
    References for this can be found in Huang and Chow (1974) 
    and Holland and Powell (2011; followed here).
    """

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume.
        Not a part of the Modified Tait EoS, currently returns 0.
        Note that this means we can't calculate Cv or Ks yet.
        gamma = V*(dP/dE)|_V = (alpha*K_S)/(Cp*rho) = (alpha*K_T)/(Cv*rho)
        """
        return 0.

    def volume(self, pressure,temperature,params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ 12
        """

        a, b, c = tait_constants(params)
        Pth=self.__rel_thermal_pressure(temperature,params)
        x = 1 - a*( 1. - pow(( 1. + b*(pressure-Pth)), -1.0*c))
        return x*params['V_0']

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        a, b, c = tait_constants(params)
        Pth=self.__rel_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        return params['K_0']*(1. + b*(psubpth))*(a + (1.-a)*pow((1. + b*(psubpth)), c))

    #calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented. 
        Returns 0. 
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.

    # Cv, heat capacity at constant volume
    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature, and volume [J/K/mol].
        Not yet implemented, returns 0.
        """
        return 0.

    def thermal_expansivity(self, pressure, temperature, volume , params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        Replace -Pth in EQ 13+1 with P-Pth for non-ambient temperature 
        """
        a, b, c = tait_constants(params)
        Pth=self.__rel_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        ein=einst(params['S_0'], params['n'])
        alpha = params['a_0']*ksi(ein/temperature)/ksi(ein/T_0)*1./((1.+b*psubpth)*(a + (1.-a)*pow((1+b*psubpth), c)))
 
        return alpha

    # Heat capacity at ambient pressure
    # N.B. Cp=-T*(d2G/dT2)|p
    def heat_capacity_p0(self,temperature,params):
        """
        Returns heat capacity at ambient pressure as a function of temperature [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in Holland and Powell, 2011
        """
        Cp = params['Cp'][0] + params['Cp'][1]*temperature + params['Cp'][2]*pow(temperature,-2.) + params['Cp'][3]*pow(temperature,-0.5)
        return Cp


    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        Not yet implemented. Returns 0.
        """
        #alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        #gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        #C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        #C_p = C_v*(1. + gr * alpha * temperature)
        return 0.


    def adiabatic_bulk_modulus(self,pressure,temperature,volume,params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  
        Not yet implemented. Returns 0.
        """
        #K_T= self.isothermal_bulk_modulus(pressure,temperature,volume,params)
        #alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        #gr = self.__grueneisen_parameter(params['V_0']/volume, params)
        #K_S = K_T*(1. + gr * alpha * temperature)
        return 0.

    def gibbs(self,pressure,temperature,params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
    # Calculate temperature and pressure integrals
        a, b, c = tait_constants(params)
        Pth=self.__rel_thermal_pressure(temperature,params)
        psubpth=pressure-Pth

        intCpdT = (params['Cp'][0]*temperature + 0.5*params['Cp'][1]*pow(temperature,2.) - params['Cp'][2]/temperature + 2.*params['Cp'][3]*np.sqrt(temperature)) - (params['Cp'][0]*T_0 + 0.5*params['Cp'][1]*T_0*T_0 - params['Cp'][2]/T_0 + 2.0*params['Cp'][3]*np.sqrt(T_0))

        intCpoverTdT = (params['Cp'][0]*np.log(temperature) + params['Cp'][1]*temperature - 0.5*params['Cp'][2]/pow(temperature,2.) - 2.0*params['Cp'][3]/np.sqrt(temperature)) - (params['Cp'][0]*np.log(T_0) + params['Cp'][1]*T_0 - 0.5*params['Cp'][2]/(T_0*T_0) - 2.0*params['Cp'][3]/np.sqrt(T_0))

        # EQ 13
        intVdP = pressure*params['V_0']*(1. - a + (a*(pow((1.-b*Pth), 1.-c) - pow((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure)))

        # Add order-disorder terms if required
        if params.has_key('landau_Tc'): # For a phase transition described by Landau term
            Gdisord=gibbs_disord_Landau(pressure, temperature, params)
        else:
            if params.has_key('BW_deltaH'): # Add Bragg-Williams disordering
                Gdisord=gibbs_disord_BW(pressure, temperature, params) - gibbs_disord_BW(P_0, T_0, params)
            else:
                Gdisord=0.0


        return params['H_0'] + intCpdT - temperature*(params['S_0'] + intCpoverTdT) + intVdP + Gdisord

    # calculate P = P(T0) + Pth
    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        return modified_tait(params['V_0']/volume, params) + \
                self.__rel_thermal_pressure(temperature, params)


    #calculate relative thermal pressure (relative to T_0), see EQ 12 - 1
    def __rel_thermal_pressure(self,T, params):
        ein=einst(params['S_0'],params['n'])
        u=ein/T
        u_0=ein/T_0
        P_th = params['a_0']*params['K_0']*ein/ksi(u_0)*((1./(np.exp(u)-1.))-(1./(np.exp(u_0)-1.)))
        return P_th


class MT(MTaitBase):
    """
    Standard MT equation of state. 
    This class currently exists for consistency with the MGD, 
    SLB and BM class set structures.
    """

