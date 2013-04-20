# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import scipy.optimize as opt
import birch_murnaghan as bm
import debye
import numpy as np
from equation_of_state import equation_of_state

## Based on Stixrude & Lithgow-Bertelloni (2005), all equation numbers refer to this paper. 
class slb_base(equation_of_state):
    def __debye_temperature(self,x,params):
        """ x = ref_vol/vol"""
        f = 1./2. * (pow(x, 2./3.) - 1.)
        a1_ii = 6. * params['ref_grueneisen'] # EQ 47
        a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
        return params['ref_Debye'] * np.sqrt(1. + a1_ii * f + 1./2. * a2_iikk*f*f) 

    def __volume_dependent_q(self, x, params):
        f = 1./2. * (pow(x, 2./3.) - 1.)
        a1_ii = 6. * params['ref_grueneisen'] # EQ 47
        a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
        nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # EQ 41
        gr = 1./6./nu_o_nu0_sq * (2.*f+1.) * ( a1_ii + a2_iikk*f )
        q = 1./9.*(18.*gr - 6. - 1./2. / nu_o_nu0_sq * (2.*f+1.)*(2.*f+1.)*a2_iikk/gr)
        return q
    
    def __isotropic_eta_s(self, x, params):
        f = 1./2. * (pow(x, 2./3.) - 1.)
        a2_s = -2.*params['ref_grueneisen'] - 2.*params['eta_0s'] # EQ 47 
        a1_ii = 6. * params['ref_grueneisen'] # EQ 47
        a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # EQ 47
        nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # EQ 41
        gr = 1./6./nu_o_nu0_sq * (2.*f+1.) * ( a1_ii + a2_iikk*f )
        eta_s = - gr - (1./2. * pow(nu_o_nu0_sq,-1.) * pow((2.*f)+1.,2.)*a2_s) # EQ 46 NOTE the typo from Stixrude 2005
        return eta_s


    def volume(self, pressure, temperature, params):
        debye_T = lambda x : self.__debye_temperature(params['ref_V']/x, params)
        gr = lambda x : self.grueneisen_parameter(pressure, temperature, x, params)
        E_th =  lambda x : debye.thermal_energy(temperature, debye_T(x), params['n']) #thermal energy at temperature T
        E_th_ref = lambda x : debye.thermal_energy(300., debye_T(x), params['n']) #thermal energy at reference temperature
      
        b_iikk= 9.*params['ref_K'] # EQ 28
        b_iikkmm= 27.*params['ref_K']*(params['K_prime']-4.) # EQ 29
        f = lambda x: 0.5*(pow(params['ref_V']/x,2./3.)-1.) # EQ 24
        func = lambda x: (1./3.)*(pow(1.+2.*f(x),5./2.))*((b_iikk*f(x)) \
            +(0.5*b_iikkmm*pow(f(x),2.))) + gr(x)*(E_th(x) - E_th_ref(x))/x - pressure #EQ 21 
    
        V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V']) 
        return V

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        x = params['ref_V'] / volume
        f = 1./2. * (pow(x, 2./3.) - 1.)
        ref_gruen = params['ref_grueneisen']
        a1_ii = 6. * ref_gruen # EQ 47
        a2_iikk = -12.*ref_gruen + 36.*ref_gruen*ref_gruen - 18.*params['q0']*ref_gruen # EQ 47
        nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * f*f # EQ 41
        return 1./6./nu_o_nu0_sq * (2.*f+1.) * ( a1_ii + a2_iikk*f )

    def isothermal_bulk_modulus(self, pressure,temperature, volume, params):

        debye_T = self.__debye_temperature(params['ref_V']/volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)

        E_th = debye.thermal_energy(temperature, debye_T, params['n']) #thermal energy at temperature T
        E_th_ref = debye.thermal_energy(300.,debye_T, params['n']) #thermal energy at reference temperature
    
        C_v = debye.heat_capacity_v(temperature, debye_T, params['n']) #heat capacity at temperature T
        C_v_ref = debye.heat_capacity_v(300.,debye_T, params['n']) #heat capacity at reference temperature

        #f =.5*(pow(params['ref_V']/V,2./3.)-1) # EQ 24

        q = self.__volume_dependent_q(params['ref_V']/volume, params)
    
        K = bm.bulk_modulus(volume, params) \
            + (gr + 1.-q)* ( gr / volume ) * (E_th - E_th_ref) \
            - ( pow(gr , 2.) / volume )*(C_v*temperature - C_v_ref*300.)   
    
        return K
    
    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        #calculate the adiabatic bulk modulus (K_S) as a function of P, T, and V
        # alpha is basically 1e-5
        K_T=self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        K_S = K_T*(1. + gr * alpha * temperature)
        return K_S

    def shear_modulus(self, pressure, temperature, volume, params):
        debye_T = self.__debye_temperature(params['ref_V']/volume, params)
        #gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        eta_s = self.__isotropic_eta_s(params['ref_V']/volume, params)

        E_th = debye.thermal_energy(temperature ,debye_T, params['n'])
        E_th_ref = debye.thermal_energy(300.,debye_T, params['n'])

        #f =.5*(pow(params['ref_V']/V,2./3.)-1.) # EQ 24
  
        if self.order==2:
            return bm.shear_modulus_second_order(volume, params) - eta_s * (E_th-E_th_ref) / volume
        elif self.order==3:
            return bm.shear_modulus_third_order(volume, params) - eta_s * (E_th-E_th_ref) / volume
        else:
            raise NotImplementedError("")

    def heat_capacity_v(self, pressure, temperature, volume, params):
        debye_T = self.__debye_temperature(params['ref_V']/volume, params)
        return debye.heat_capacity_v(temperature, debye_T,params['n'])

    def heat_capacity_p(self, pressure, temperature, volume, params):
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v*(1. + gr * alpha * temperature)
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = gr * C_v / K / volume
        return alpha
    
    
    
class slb3(slb_base):
    def __init__(self):
        self.order=3

class slb2(slb_base):
    def __init__(self):
        self.order=2
    
