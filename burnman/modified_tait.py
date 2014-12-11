# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


#TO DO: Correct heat capacity, volume where internal order-disorder is implemented (Landau and Bragg-Williams models)

import numpy as np
import scipy.optimize as opt

import burnman.equation_of_state as eos
import burnman.einstein as einstein
from burnman.constants import R

T_0=298.15 # Standard temperature = 25 C
P_0=1.e5 # Standard pressure = 1.e5 Pa

# see Holland and Powell, 2011
def einstein_temperature(S, n):
    """
    Empirical Einstein temperature
    base of p.346, para.1
    """
    return 10636./(S/n + 6.44)

# see Holland and Powell, 2011
def tait_constants(params):
    a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
    b=params['Kprime_0']/params['K_0'] - params['Kdprime_0']/(1. + params['Kprime_0'])
    c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])
    return a, b, c

def landau_ordering(P, T, params):
    # From Holland and Powell, 1996, corrected using
    # landaunote.pdf on Tim Holland's esc web page
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
    # From Holland and Powell, 1996, corrected using
    # landaunote.pdf on Tim Holland's esc web page
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
        
    # Vt is defined as a "function shaping the excess volume of disordering" in landaunote.pdf.
    # A sensible function to use in this case would be the expression in brackets in  EQ13 of Holland and Powell (2011) 
    # Here we follow tc337L, where Vt=1 ("the simplest way to proceed"; Roger Powell, pers. comm. 2014/08/22). 
    Vt=1.# EQ 13 would give: 1. - a + (a*(pow((1.-b*Pth), 1.-c) - pow((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure))
    Gdisord=params['landau_Tc']*params['landau_Smax']*(np.power(Q_0,2) - np.power(Q_0,6)/3.0) - params['landau_Smax']*(Tcstar*np.power(Q,2) - params['landau_Tc']*np.power(Q,6)/3.0) - T*(params['landau_Smax']*(np.power(Q_0,2) - np.power(Q,2))) + P*(params['landau_Vmax']*np.power(Q_0,2)*Vt)

    return Gdisord

def entropy_disorder_Landau(P, T, params):
    # N.B. Assumes Vt==1, see above
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
    Sdisord=params['landau_Smax']*(np.power(Q_0, 2.) - np.power(Q, 2.))
    return Sdisord

def enthalpy_disorder_Landau(P, T, params):
    # N.B. Assumes Vt==1, see above
    Tcstar, Q_0, Q = landau_ordering(P, T, params)
    Hdisord=gibbs_disorder_Landau(P, T, params) + T*entropy_disorder_Landau(P, T, params)
    return Hdisord

def heat_capacity_p_disorder_Landau(P, T, params):
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
    return params['BW_deltaH'] - params['BW_factor']*T*deltaS + P*params['BW_deltaV'] + params['BW_factor']*R*T*(lnxdisord(n,Q) - lnxord(n,Q)) + (2.*Q - 1.)*W

# Energy of disordering from Bragg-Williams symmetric model; see Holland and Powell, 1996
def gibbs_disorder_BW(P, T, params):
    n=params['BW_n']
    deltaS=entropydisorder(n)
    Q=opt.fsolve(equilibrium_Q, 0.999995, args=(deltaS, P, T, params))[0]
    W=params['BW_W'] + P*params['BW_Wv']
    ideal=(1.-Q)*(params['BW_deltaH'] - params['BW_factor']*T*entropydisorder(n) + P*params['BW_deltaV'] + params['BW_factor']*R*T*lnxdisord(n,Q)) + params['BW_factor']*Q*(R*T*lnxord(n,Q))
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


# see Holland and Powell, 2011
def modified_tait(x, params):
    """
    equation for the modified Tait equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    EQ 2
    """
    a, b, c = tait_constants(params)

    return (np.power((x + a - 1.) / a, -1./c) - 1.)/b


class MT(eos.EquationOfState):
    """
    Base class for a generic modified Tait equation of state.  
    References for this can be found in Huang and Chow (1974) 
    and Holland and Powell (2011; followed here).
    """

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume.
        """
        alpha = self.thermal_expansivity (pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus (pressure, temperature, volume, params)
        C_V = self.heat_capacity_v( pressure, temperature, volume, params)
        return alpha * K_T * volume / C_V

    def volume(self, pressure,temperature,params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ 12
        """
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)
        x = 1 - a*( 1. - np.power(( 1. + b*(pressure-Pth)), -1.0*c))
        return x*params['V_0']

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        return params['K_0']*(1. + b*(psubpth))*(a + (1.-a)*np.power((1. + b*(psubpth)), c))

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
        """
        C_p=self.heat_capacity_p(pressure, temperature, volume, params)
        V=self.volume(pressure,temperature,params)
        alpha=self.thermal_expansivity(pressure, temperature, volume , params)
        K_T=self.isothermal_bulk_modulus(pressure,temperature,volume, params)
        return C_p - V*temperature*alpha*alpha*K_T

    def thermal_expansivity(self, pressure, temperature, volume , params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        Replace -Pth in EQ 13+1 with P-Pth for non-ambient temperature 
        """
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        einstein_T=einstein_temperature(params['S_0'], params['n'])
        C_V0 = einstein.heat_capacity_v( T_0, einstein_T, params['n'] )
        C_V =  einstein.heat_capacity_v(temperature, einstein_T,params['n'])
        alpha = params['a_0'] * (C_V/C_V0) *1./((1.+b*psubpth)*(a + (1.-a)*np.power((1+b*psubpth), c)))
 
        return alpha

    def heat_capacity_p0(self,temperature,params):
        """
        Returns heat capacity at ambient pressure as a function of temperature [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in Holland and Powell, 2011
        """
        Cp = params['Cp'][0] + params['Cp'][1]*temperature + params['Cp'][2]*np.power(temperature,-2.) + params['Cp'][3]*np.power(temperature,-0.5)
        return Cp


    def heat_capacity_p_einstein(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v*(1. + gr * alpha * temperature)
        return C_p


    def adiabatic_bulk_modulus(self,pressure,temperature,volume,params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  
        """
        K_T= self.isothermal_bulk_modulus(pressure,temperature,volume,params)
        alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        K_S = K_T*(1. + gr * alpha * temperature)
        return K_S

    def gibbs_free_energy(self,pressure,temperature,params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
       # Calculate temperature and pressure integrals
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)

        psubpth=pressure-Pth

        # EQ 13
        intVdP = pressure*params['V_0']*(1. - a + (a*(np.power((1.-b*Pth), 1.-c) - np.power((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure)))

        # Add order-disorder terms if required
        if params.has_key('landau_Tc'): # For a phase transition described by Landau term
            Gdisord=gibbs_disorder_Landau(pressure, temperature, params)
        else:
            if params.has_key('BW_deltaH'): # Add Bragg-Williams disordering
                Gdisord=gibbs_disorder_BW(pressure, temperature, params) - gibbs_disorder_BW(P_0, T_0, params)
            else:
                Gdisord=0.0

        return params['H_0'] + self.__intCpdT(temperature, params) - temperature*(params['S_0'] + self.__intCpoverTdT(temperature, params)) + intVdP + Gdisord


    def entropy(self,pressure,temperature,params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)

        einstein_T=einstein_temperature(params['S_0'], params['n'])
        ksi_over_ksi_0=einstein.heat_capacity_v( temperature, einstein_T, params['n'] )/einstein.heat_capacity_v( T_0, einstein_T, params['n'] )

        dintVdpdx=(params['V_0']*params['a_0']*params['K_0']*a*ksi_over_ksi_0)*(np.power((1.+b*(pressure-Pth)), 0.-c) - np.power((1.-b*Pth), 0.-c))

        # Add order-disorder terms if required
        if params.has_key('landau_Tc'): # For a phase transition described by Landau term
            Sdisord=entropy_disorder_Landau(pressure, temperature, params)
        else:
            if params.has_key('BW_deltaH'): # Add Bragg-Williams disordering
                Sdisord=entropy_disorder_BW(pressure, temperature, params) - entropy_disorder_BW(P_0, T_0, params)
            else:
                Sdisord=0.0

        return params['S_0'] + self.__intCpoverTdT(temperature, params) + dintVdpdx + Sdisord

    def enthalpy(self, pressure, temperature, params):
        """
        Returns the enthalpy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        gibbs=self.gibbs_free_energy(pressure,temperature,params)
        entropy=self.entropy(pressure,temperature,params)
        
        # Add order-disorder terms if required
        if params.has_key('landau_Tc'): # For a phase transition described by Landau term
            Hdisord=enthalpy_disorder_Landau(pressure, temperature, params)
        else:
            if params.has_key('BW_deltaH'): # Add Bragg-Williams disordering
                Hdisord=enthalpy_disorder_BW(pressure, temperature, params) - enthalpy_disorder_BW(P_0, T_0, params)
            else:
                Hdisord=0.0

        return gibbs + temperature*entropy + Hdisord

    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)

        einstein_T=einstein_temperature(params['S_0'], params['n'])
        ksi_over_ksi_0=einstein.heat_capacity_v( temperature, einstein_T, params['n'] )/einstein.heat_capacity_v( T_0, einstein_T, params['n'] )

        dSdT=params['V_0']*params['K_0']*np.power((ksi_over_ksi_0*params['a_0']),2.0)*(np.power((1.+b*(pressure-Pth)), -1.-c) - np.power((1.-b*Pth), -1.-c))

        # Add order-disorder terms if required
        if params.has_key('landau_Tc'): # For a phase transition described by Landau term
            Cpdisord=heat_capacity_p_disorder_Landau(pressure, temperature, params)
        else:
            Cpdisord=0.0

        return self.heat_capacity_p0(temperature,params) + temperature*dSdT + Cpdisord

    # calculate P = P(T0) + Pth
    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        return modified_tait(params['V_0']/volume, params) + \
               self.__relative_thermal_pressure(temperature, params)
                  


    #thermal pressure (relative to T_0), see EQ 12 - 1
    def __thermal_pressure(self,T,params):

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is 
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR.  I don't know why they don't use that, but anyhow...

        einstein_T=einstein_temperature(params['S_0'],params['n'])
        E_th = einstein.thermal_energy( T, einstein_T, params['n'] )
        C_V0 = einstein.heat_capacity_v( T_0, einstein_T, params['n'] )
        P_th = params['a_0']*params['K_0'] / C_V0 * E_th
        return P_th

    #calculate relative thermal pressure (relative to T_0), see EQ 12 - 1
    def __relative_thermal_pressure( self, T, params):
        return self.__thermal_pressure(T, params) - \
               self.__thermal_pressure(T_0, params)

    def __intCpdT (self, temperature, params):
        return (params['Cp'][0]*temperature + 0.5*params['Cp'][1]*np.power(temperature,2.) - params['Cp'][2]/temperature + 2.*params['Cp'][3]*np.sqrt(temperature)) - (params['Cp'][0]*T_0 + 0.5*params['Cp'][1]*T_0*T_0 - params['Cp'][2]/T_0 + 2.0*params['Cp'][3]*np.sqrt(T_0))

    def __intCpoverTdT (self, temperature, params):
        return (params['Cp'][0]*np.log(temperature) + params['Cp'][1]*temperature - 0.5*params['Cp'][2]/np.power(temperature,2.) - 2.0*params['Cp'][3]/np.sqrt(temperature)) - (params['Cp'][0]*np.log(T_0) + params['Cp'][1]*T_0 - 0.5*params['Cp'][2]/(T_0*T_0) - 2.0*params['Cp'][3]/np.sqrt(T_0))
