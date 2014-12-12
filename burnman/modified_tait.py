# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt

import burnman.equation_of_state as eos
import burnman.einstein as einstein
from burnman.endmemberdisorder import *
import constants

T_0=298.15 # Standard temperature = 25 C
P_0=1.e5 # Standard pressure = 1.e5 Pa


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
        a, b, c = self.__tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)
        x = 1 - a*( 1. - np.power(( 1. + b*(pressure-Pth)), -1.0*c))
        return x*params['V_0']

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        a, b, c = self.__tait_constants(params)
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
        a, b, c = self.__tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)
        psubpth=pressure-Pth
        einstein_T=self.__einstein_temperature(params['S_0'], params['n'])
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
        Returns heat capacity at constant pressure at the pressure, temperature, and volume, using the C_v and Einstein model [J/K/mol]
        WARNING: Only for comparison with internally self-consistent C_p
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
        a, b, c = self.__tait_constants(params)
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
        a, b, c = self.__tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)

        einstein_T=self.__einstein_temperature(params['S_0'], params['n'])
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
        a, b, c = self.__tait_constants(params)
        Pth=self.__relative_thermal_pressure(temperature,params)

        einstein_T=self.__einstein_temperature(params['S_0'], params['n'])
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
        return self.__modified_tait(params['V_0']/volume, params) + \
               self.__relative_thermal_pressure(temperature, params)
                  

    def __einstein_temperature(self, S, n):
        """
        Empirical Einstein temperature
        Holland and Powell, 2011; base of p.346, para.1
        """
        return 10636./(S/n + 6.44)
    
    def __tait_constants(self, params):
        """
        returns parameters for the modified Tait equation of state
        derived from K_T and its two first pressure derivatives
        EQ 4 from Holland and Powell, 2011
        """
        a=(1.+params['Kprime_0'])/(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])
        b=params['Kprime_0']/params['K_0'] - params['Kdprime_0']/(1. + params['Kprime_0'])
        c=(1. + params['Kprime_0'] + params['K_0']*params['Kdprime_0'])/(params['Kprime_0']*params['Kprime_0'] + params['Kprime_0'] - params['K_0']*params['Kdprime_0'])
        return a, b, c
    
    def __modified_tait(self, x, params):
        """
        equation for the modified Tait equation of state, returns
        pressure in the same units that are supplied for the reference bulk
        modulus (params['K_0'])
        EQ 2 from Holland and Powell, 2011
        """
        a, b, c = self.__tait_constants(params)
        
        return (np.power((x + a - 1.) / a, -1./c) - 1.)/b
    
    def __thermal_pressure(self,T,params):
        """
        Returns thermal pressure [Pa] as a function of T [K] 
        EQ 12 - 1 of Holland and Powell, 2011 
        """

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is 
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR.  I don't know why they don't use that, but anyhow...

        einstein_T=self.__einstein_temperature(params['S_0'],params['n'])
        E_th = einstein.thermal_energy( T, einstein_T, params['n'] )
        C_V0 = einstein.heat_capacity_v( T_0, einstein_T, params['n'] )
        P_th = params['a_0']*params['K_0'] / C_V0 * E_th
        return P_th

    def __relative_thermal_pressure( self, T, params):
        """
        Returns relative thermal pressure [Pa] as a function of T-T_0 [K] 
        EQ 12 - 1 of Holland and Powell, 2011 
        """
        return self.__thermal_pressure(T, params) - \
               self.__thermal_pressure(T_0, params)

    def __intCpdT (self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return (params['Cp'][0]*temperature + 0.5*params['Cp'][1]*np.power(temperature,2.) - params['Cp'][2]/temperature + 2.*params['Cp'][3]*np.sqrt(temperature)) - (params['Cp'][0]*T_0 + 0.5*params['Cp'][1]*T_0*T_0 - params['Cp'][2]/T_0 + 2.0*params['Cp'][3]*np.sqrt(T_0))

    def __intCpoverTdT (self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return (params['Cp'][0]*np.log(temperature) + params['Cp'][1]*temperature - 0.5*params['Cp'][2]/np.power(temperature,2.) - 2.0*params['Cp'][3]/np.sqrt(temperature)) - (params['Cp'][0]*np.log(T_0) + params['Cp'][1]*T_0 - 0.5*params['Cp'][2]/(T_0*T_0) - 2.0*params['Cp'][3]/np.sqrt(T_0))
