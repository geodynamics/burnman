# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import equation_of_state as eos
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt
import birch_murnaghan as bm
import debye


class mie_grueneisen_debye(eos.equation_of_state):
    """
    """

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        return self.__grueneisen_parameter(params['ref_V']/volume, params)

    #invert the mie-grueneisen-debye eqn of state with thermal corrections
    #for the volume
    def volume(self, pressure,temperature,params):
        func = lambda x: bm.birch_murnaghan(params['ref_V']/x, params) + \
            self.__thermal_pressure(temperature, x, params) - \
            self.__thermal_pressure(300., x, params) - pressure
        V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V'])
        return V

    #calculate the mgd bulk modulus (K_T) as a function of P, T, and V
    def isothermal_bulk_modulus(self, pressure,temperature,volume, params):
        K_T = bm.bulk_modulus(volume, params) + \
            self.__thermal_bulk_modulus(temperature,volume, params) - \
            self.__thermal_bulk_modulus(300.,volume, params)  #EQB13
        return K_T

    #calculate the mgd shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        mu = bm.shear_modulus_second_order(volume,params) + \
            self.__thermal_shear_modulus(temperature,volume, params) - \
            self.__thermal_shear_modulus(300.,volume, params) # EQ B11
        return mu

    #heat capacity at constant volume
    def heat_capacity_v(self, pressure, temperature, volume, params):
         Debye_T = self.__debye_temperature(params['ref_V']/volume, params)
         C_v = debye.heat_capacity_v(temperature, Debye_T, params['n'])
         return C_v

    def thermal_expansivity(self, pressure, temperature, volume , params):
        C_v = self.heat_capacity_v(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['ref_V']/volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume ,params)
        alpha = gr * C_v / K / volume
        return alpha

    #heat capacity at constant pressure
    def heat_capacity_p(self,pressure, temperature,volume,params):
        alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['ref_V']/volume, params)
        C_v = self.heat_capacity_v(pressure,temperature,volume,params)
        C_p = C_v*(1. + gr * alpha * temperature)
        return C_p

    #calculate the adiabatic bulk modulus (K_S) as a function of P, T, and V
    # alpha is basically 1e-5
    def adiabatic_bulk_modulus(self,pressure,temperature,volume,params):
        K_T= self.isothermal_bulk_modulus(pressure,temperature,volume,params)
        alpha = self.thermal_expansivity(pressure,temperature,volume,params)
        gr = self.__grueneisen_parameter(params['ref_V']/volume, params)
        K_S = K_T*(1. + gr * alpha * temperature)
        return K_S

    #do the forward problem of the mie-grueneisen-debye equation of state, i.e.
    #get pressure from volume
    def pressure(self, temperature, volume, params):
        return bm.birch_murnaghan(params['ref_V']/volume, params) + \
                self.__thermal_pressure(temperature,volume, params) - \
                self.__thermal_pressure(300.,volume, params)

    #calculate the thermal correction to the shear modulus as a function of V, T
    def __thermal_shear_modulus(self, T, V, params):
        gr = self.__grueneisen_parameter(params['ref_V']/V, params)
        Debye_T = self.__debye_temperature(params['ref_V']/V, params) 
        mu_th= 3./5. * ( self.__thermal_bulk_modulus(T,V,params) - \
                 6*debye.R*T*params['n']/V * gr * debye.debye_fn(Debye_T/T) ) # EQ B10
        return mu_th

    #compute the Debye temperature in K.  Takes the
    #parameter x, which is ref_V/V (molar volumes).
    #Depends on the reference grueneisen parameter,
    #the reference Debye temperature, and the factor
    #q0, see Matas eq B6
    def __debye_temperature(self, x, params):
        return params['ref_Debye']*np.exp((params['ref_grueneisen']- \
            self.__grueneisen_parameter(x, params))/params['q0'])

    #compute the grueneisen parameter with depth, according
    #to q0.  Takes x=ref_V/V. See Matas eq B6
    def __grueneisen_parameter(self, x, params):
        return params['ref_grueneisen']*pow(1./x, params['q0'])

    #calculate isotropic thermal pressure, see
    # Matas et. al. (2007) eq B4
    def __thermal_pressure(self,T,V, params):
        Debye_T = self.__debye_temperature(params['ref_V']/V, params) 
        gr = self.__grueneisen_parameter(params['ref_V']/V, params)
        P_th = gr * debye.thermal_energy(T,Debye_T, params['n'])/V
        return P_th


    #calculate the thermal correction for the mgd
    #bulk modulus (see matas et al, 2007)
    def __thermal_bulk_modulus(self, T,V,params):
        gr = self.__grueneisen_parameter(params['ref_V']/V, params)
        Debye_T = self.__debye_temperature(params['ref_V']/V, params) 
        K_th = 3.*params['n']*debye.R*T/V * gr * \
            ((1. - params['q0'] - 3.*gr)*debye.debye_fn(Debye_T/T)+3.*gr*(Debye_T/T)/(np.exp(Debye_T/T) - 1.)) # EQ B5
        return K_th


