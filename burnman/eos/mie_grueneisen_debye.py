# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import warnings

import equation_of_state as eos
import birch_murnaghan as bm
import burnman.debye as debye
import burnman.constants as constants

class MGDBase(eos.EquationOfState):
    """
    Base class for a generic finite-strain-mie-grueneisen-debye
    equation of state.  References for this can be found in many
    places, such as Shim, Duffym and Kenichi (2002) and Jackson and Rigedn
    (1996).  Here we mostly follow the appendices of Matas et al (2007)
    """

    #def __init__(self):
    #    pass

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume (EQ B6)
        """
        return self.__grueneisen_parameter(params['V_0']/volume, params)

    def volume(self, pressure,temperature,params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ B7
        """
        T_0 = self.reference_temperature( params )
        func = lambda x: bm.birch_murnaghan(params['V_0']/x, params) + \
            self.__thermal_pressure(temperature, x, params) - \
            self.__thermal_pressure(T_0, x, params) - pressure
        V = opt.brentq(func, 0.5*params['V_0'], 1.5*params['V_0'])
        return V

    def pressure(self, temperature, volume, params):
        P = bm.birch_murnaghan(params['V_0']/volume, params) + \
            self.__thermal_pressure(temperature, volume, params) - \
            self.__thermal_pressure(300., volume, params)
        return P

    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B8
        """
        T_0 = self.reference_temperature( params )
        K_T = bm.bulk_modulus(volume, params) + \
            self.__thermal_bulk_modulus(temperature,volume, params) - \
            self.__thermal_bulk_modulus(T_0,volume, params)  #EQB13
        return K_T

    #calculate the mgd shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ B11
        """
        T_0 = self.reference_temperature( params )
        if self.order==2:
            return bm.shear_modulus_second_order(volume,params) + \
                self.__thermal_shear_modulus(temperature,volume, params) - \
                self.__thermal_shear_modulus(T_0,volume, params) # EQ B11
        elif self.order==3:
            return bm.shear_modulus_third_order(volume,params) + \
                self.__thermal_shear_modulus(temperature,volume, params) - \
                self.__thermal_shear_modulus(T_0,volume, params) # EQ B11
        else:
            raise NotImplementedError("")

    #heat capacity at constant volume
    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature, and volume [J/K/mol]
        """
        Debye_T = self.__debye_temperature(params['V_0']/volume, params)
        C_v = debye.heat_capacity_v(temperature, Debye_T, params['n'])
        return C_v

    def thermal_expansivity(self, pressure, temperature, volume , params):
        """
        Returns thermal expansivity at the pressure, temperature, and volume [1/K]
        """
        C_v = self.heat_capacity_v(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['V_0']/volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume ,params)
        alpha = gr * C_v / K / volume
        return alpha

    #heat capacity at constant pressure
    def heat_capacity_p(self,pressure, temperature,volume,params):
        """
        Returns heat capacity at constant pressure at the pressure, temperature, and volume [J/K/mol]
        """
        alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['V_0']/volume, params)
        C_v = self.heat_capacity_v(pressure,temperature,volume,params)
        C_p = C_v*(1. + gr * alpha * temperature)
        return C_p

    def adiabatic_bulk_modulus(self,pressure,temperature,volume,params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ D6
        """
        K_T= self.isothermal_bulk_modulus(pressure,temperature,volume,params)
        alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['V_0']/volume, params)
        K_S = K_T*(1. + gr * alpha * temperature)
        return K_S

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        T_0 = self.reference_temperature( params )
        return bm.birch_murnaghan(params['V_0']/volume, params) + \
                self.__thermal_pressure(temperature,volume, params) - \
                self.__thermal_pressure(T_0,volume, params)

    #calculate the thermal correction to the shear modulus as a function of V, T
    def __thermal_shear_modulus(self, T, V, params):
        gr = self.__grueneisen_parameter(params['V_0']/V, params)
        Debye_T = self.__debye_temperature(params['V_0']/V, params)
        G_th= 3./5. * ( self.__thermal_bulk_modulus(T,V,params) - \
                 6*constants.gas_constant*T*params['n']/V * gr * debye.debye_fn(Debye_T/T) ) # EQ B10
        return G_th

    #compute the Debye temperature in K.  Takes the
    #parameter x, which is V_0/V (molar volumes).
    #Depends on the reference grueneisen parameter,
    #the reference Debye temperature, and the factor
    #q_0, see Matas eq B6
    def __debye_temperature(self, x, params):
        return params['Debye_0']*np.exp((params['grueneisen_0']- \
            self.__grueneisen_parameter(x, params))/params['q_0'])

    #compute the grueneisen parameter with depth, according
    #to q_0.  Takes x=V_0/V. See Matas eq B6
    def __grueneisen_parameter(self, x, params):
        return params['grueneisen_0']*pow(1./x, params['q_0'])

    #calculate isotropic thermal pressure, see
    # Matas et. al. (2007) eq B4
    def __thermal_pressure(self,T,V, params):
        Debye_T = self.__debye_temperature(params['V_0']/V, params)
        gr = self.__grueneisen_parameter(params['V_0']/V, params)
        P_th = gr * debye.thermal_energy(T,Debye_T, params['n'])/V
        return P_th


    #calculate the thermal correction for the mgd
    #bulk modulus (see matas et al, 2007)
    def __thermal_bulk_modulus(self, T,V,params):
        gr = self.__grueneisen_parameter(params['V_0']/V, params)
        Debye_T = self.__debye_temperature(params['V_0']/V, params)
        K_th = 3.*params['n']*constants.gas_constant*T/V * gr * \
            ((1. - params['q_0'] - 3.*gr)*debye.debye_fn(Debye_T/T)+3.*gr*(Debye_T/T)/(np.exp(Debye_T/T) - 1.)) # EQ B5
        return K_th


    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        #if G and Gprime are not included this is presumably deliberate,
        #as we can model density and bulk modulus just fine without them,
        #so just add them to the dictionary as nans
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')
  
        #check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0', 'molar_mass', 'n', 'Debye_0', 'grueneisen_0', 'q_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)
        
        #now check that the values are reasonable.  I mostly just
        #made up these values from experience, and we are only 
        #raising a warning.  Better way to do this? [IR]
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-3:
            warnings.warn( 'Unusual value for V_0' , stacklevel=2 )
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn( 'Unusual value for K_0', stacklevel=2 )
        if params['Kprime_0'] < -5. or params['Kprime_0'] > 10.:
            warnings.warn( 'Unusual value for Kprime_0', stacklevel=2 )
        if params['G_0'] < 0. or params['G_0'] > 1.e13:
            warnings.warn( 'Unusual value for G_0', stacklevel=2 )
        if params['Gprime_0'] < -5. or params['Gprime_0'] > 10.:
            warnings.warn( 'Unusual value for Gprime_0', stacklevel=2 )
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 1.:
            warnings.warn( 'Unusual value for molar_mass' , stacklevel=2)
        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn( 'Unusual value for n' , stacklevel=2)
        if params['Debye_0'] < 1. or params['Debye_0'] > 10000.:
            warnings.warn( 'Unusual value for Debye_0' , stacklevel=2)
        if params['grueneisen_0'] < 0. or params['grueneisen_0'] > 10.:
            warnings.warn( 'Unusual value for grueneisen_0' , stacklevel=2)
        if params['q_0'] < -10. or params['q_0'] > 10.:
            warnings.warn( 'Unusual value for q_0' , stacklevel=2)


class MGD3(MGDBase):
    """
    MGD equation of state with third order finite strain expansion for the
    shear modulus (this should be preferred, as it is more thermodynamically
    consistent.
    """
    def __init__(self):
        self.order=3


class MGD2(MGDBase):
    """
    MGD equation of state with second order finite strain expansion for the
    shear modulus.  In general, this should not be used, but sometimes
    shear modulus data is fit to a second order equation of state.  In that
    case, you should use this.  The moral is, be careful!
    """
    def __init__(self):
        self.order=2

