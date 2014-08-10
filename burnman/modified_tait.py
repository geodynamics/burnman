# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt

import burnman.equation_of_state as eos

T_0=298.15 # Standard temperature = 25 C

def einst(S, n):
    "
    Einstein temperature
    base of p.346, para.1
    "
    return 10636./(S/n + 6.44)

def ksi(u):
    "
    Einstein function to describe behaviour of ak
    EQ 11+1
    "
    return pow(u,2.)*exp(u)/pow((exp(u)-1.), 2.)

def modified_tait(x, params):
    """
    equation for the modified Tait equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    EQ 2
    """

    a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
    b=(params['Kprime_0'])/(params['K_0'] - (params['Kdprime_0'])/(1. + params['Kprime_0']))
    c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])

    return (pow((x + a - 1.) / a, -1./c) - 1.)/b


class MTaitBase(eos.EquationOfState):
    """
    Base class for a generic modified Tait equation of state.  
    References for this can be found in Huang and Chow (1974) 
    and Holland and Powell (2011; followed here).
    """

    #def __init__(self):
    #    pass

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
        a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
        b=(params['Kprime_0'])/(params['K_0'] - (params['Kdprime_0'])/(1. + params['Kprime_0']))
        c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])

        Pth=self.__rel_thermal_pressure(temperature,params)
        x = 1 - a*( 1. - pow(( 1. + b*(pressure-Pth)), -1.0*c))
        return x*params['V_0']

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
        b=(params['Kprime_0'])/(params['K_0'] - (params['Kdprime_0'])/(1. + params['Kprime_0']))
        c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])

        Pth=self.__rel_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        return params['K_0']*(1. + b(psubpth))*(a + (1.-a)*pow((1. + b(psubpth)), c))

    #calculate the mgd shear modulus as a function of P, V, and T
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
        a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
        b=(params['Kprime_0'])/(params['K_0'] - (params['Kdprime_0'])/(1. + params['Kprime_0']))
        c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])
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
        a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
        b=(params['Kprime_0'])/(params['K_0'] - (params['Kdprime_0'])/(1. + params['Kprime_0']))
        c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])
        Pth=self.__rel_thermal_pressure(temperature,params)
        psubpth=pressure-Pth

        intCpdT = (params['Cp'][0]*temperature + 0.5*params['Cp'][1]*pow(temperature,2.) - params['Cp'][2]/temperature + 2.*params['Cp'][3]*sqrt(temperature)) - (params['Cp'][0]*T_0 + 0.5*params['Cp'][1]*T_0*T_0 - params['Cp'][2]/T_0 + 2.0*params['Cp'][3]*sqrt(T_0))

        intCpoverTdT = (params['Cp'][0]*log(temperature) + params['Cp'][1]*temperature - 0.5*params['Cp'][2]/pow(temperature,2.) - 2.0*params['Cp'][3]/sqrt(temperature)) - (params['Cp'][0]*log(T_0) + params['Cp'][1]*T_0 - 0.5*params['Cp'][2]/(T_0*T_0) - 2.0*params['Cp'][3]/sqrt(T_0))

        # EQ 13
        intVdP = pressure*params['V_0']*(1. - a + (a*(pow((1.-b*Pth), 1.-c) - pow((1. + b*(pressure-Pth)), 1.-c))/(b*(c-1.)*pressure)))
        
        return params['H_0'] + intCpdT - temperature*(params['S_0'] + intCpoverTdT) + intVdP

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
        P_th = params['a_0']*params['k_0']*ein/ksi(T_0)*((1./(exp(u))-1.))-(1./(exp(u_0))-1.))
        return P_th


class MT(MTaitBase):
    """
    Standard MT equation of state. 
    This class currently exists for consistency with the MGD, 
    SLB and BM class set structures.
    """
