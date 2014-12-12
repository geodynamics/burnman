# BurnMan - a lower mantle toolkit
# Copyright (C) 2012-2014, Myhill, R., Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.


#TO DO: Correct heat capacity, volume where internal order-disorder is implemented (Landau and Bragg-Williams models)

import numpy as np
import scipy.optimize as opt

import burnman.equation_of_state as eos
import constants

T_0=298.15 # Standard temperature = 25 C
P_0=1.e5 # Standard pressure = 1.e5 Pa

def cork_variables(cork, cork_P, cork_T, temperature):
    a=cork[0][0]*cork_T**(2.5)/cork_P + cork[0][1]*cork_T**(1.5)/cork_P*temperature
    b=cork[1][0]*cork_T/cork_P
    c=cork[2][0]*cork_T/cork_P**(1.5) + cork[2][1]/cork_P**(1.5)*temperature
    d=cork[3][0]*cork_T/cork_P**(2.0) + cork[3][1]/cork_P**(2.0)*temperature
    return [a, b, c, d]

class CORK(eos.EquationOfState):
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
        return 0.

    def volume(self, pressure,temperature,params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        Eq. 7 in Holland and Powell, 1991
        """
        cork=cork_variables(params['cork_params'], params['cork_P'], params['cork_T'], temperature)
        V = R*temperature/pressure + (cork[1] - cork[0]*R*np.sqrt(temperature)/((R*temperature + cork[1]*pressure)*(R*temperature + 2.*cork[1]*pressure)) + cork[2]*np.sqrt(pressure) + cork[3]*pressure)
        return V

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        return 0.

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
        return 0.

    def thermal_expansivity(self, pressure, temperature, volume , params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        Replace -Pth in EQ 13+1 with P-Pth for non-ambient temperature 
        """
        return 0.

    # Heat capacity at ambient pressure
    def heat_capacity_p0(self,temperature,params):
        """
        Returns heat capacity at ambient pressure as a function of temperature [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in Holland and Powell, 2011
        """
        Cp = params['Cp'][0] + params['Cp'][1]*temperature + params['Cp'][2]*np.power(temperature,-2.) + params['Cp'][3]*np.power(temperature,-0.5)
        return Cp


    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        return 0


    def adiabatic_bulk_modulus(self,pressure,temperature,volume,params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  
        """
        return 0.

    def gibbs_free_energy(self,pressure,temperature,params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
       # Calculate temperature and pressure integrals
        intCpdT = (params['Cp'][0]*temperature + 0.5*params['Cp'][1]*np.power(temperature,2.) - params['Cp'][2]/temperature + 2.*params['Cp'][3]*np.sqrt(temperature)) - (params['Cp'][0]*T_0 + 0.5*params['Cp'][1]*T_0*T_0 - params['Cp'][2]/T_0 + 2.0*params['Cp'][3]*np.sqrt(T_0))

        intCpoverTdT = (params['Cp'][0]*np.log(temperature) + params['Cp'][1]*temperature - 0.5*params['Cp'][2]/np.power(temperature,2.) - 2.0*params['Cp'][3]/np.sqrt(temperature)) - (params['Cp'][0]*np.log(T_0) + params['Cp'][1]*T_0 - 0.5*params['Cp'][2]/(T_0*T_0) - 2.0*params['Cp'][3]/np.sqrt(T_0))

        if params['cork_T'] == 0:
            RTlnf = 0.
        else:
            cork=cork_variables(params['cork_params'], params['cork_P'], params['cork_T'], temperature)
        
            RTlnf = R*temperature*np.log(1e-5*pressure) + cork[1]*pressure + cork[0]/(cork[1]*np.sqrt(temperature))*(np.log(R*temperature + cork[1]*pressure) - np.log(R*temperature + 2.*cork[1]*pressure)) + 2./3.*cork[2]*pressure*np.sqrt(pressure) + cork[3]/2.*pressure*pressure  # Eq. 8 in Holland and Powell, 1991

        return params['H_0'] + intCpdT - temperature*(params['S_0'] + intCpoverTdT) + RTlnf

    # calculate P = P(T0) + Pth
    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        return 0.
                  
