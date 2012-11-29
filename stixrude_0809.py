#system libs:
import numpy
import scipy.optimize as opt
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt

#own libs:
from tools import *
import birch_murnaghan as bm
import mie_grueneisen_debye as mgd

# TODO: add up weight percent and check <100 and tell them how much

gas_constant=8.314462175

def stx_volume(p,T,params):
 	
    P_th = lambda x: mgd.mgd_thermal_pressure(x, T, params) #From Stixrude 2005: delta U (T)*gamma*ro (1/V)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(x, 300., params)#From Stixrude 2005: delta U (300 K)*gamma*ro (1/V)
      
    b_iikk= 9.*params['ref_K']
    b_iikkmm= 27.*params['ref_K']*(params['K_prime']-4.)

    f = lambda x: 0.5*(pow(params['ref_V']/x,2./3.)-1.)
	
    func = lambda x: (1./3.)*(pow(1.+2.*f(x),5./2.))*((b_iikk*f(x)) \
    		+(0.5*b_iikkmm*pow(f(x),2.)))*1.e9 + P_th(x) - P_th_ref(x) - p #EQ 21 Stixrude 2005
    
    V = opt.brentq(func, 0.09*params['ref_V'], 3.*params['ref_V']) 
    return V

def stx_shear_modulus(p,T,V,params):

        # low T:
        # C_v = 234. * n * Av * boltzmann_constant (T/theta_D) ^ 3.  (3.43 from Poirier/EarthInterior)
        # high T:
        # C_v = 3. * n * gas_constant
        # n = number of moles
    #C_v = 3. * n * gas_constant # in J * mol / K
    P_th = lambda x: mgd.mgd_thermal_pressure(x, T, params)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(x, 300., params)
        
    f =.5*(pow(params['ref_V']/V,2./3.)-1) 
    gamma = params['ref_grueneisen'] * (pow(params['ref_V']/V,-params['q0']))
    a2_s = -2.*params['ref_grueneisen'] - 2.*params['eta_0s'] # eq 47 (Stixrude)
    a1_ii = 6. * params['ref_grueneisen'] # eq 47
    a2_iikk = -12.*params['ref_grueneisen']+36.*pow(params['ref_grueneisen'],2.) - 18.*params['q0']*params['ref_grueneisen'] # eq 47
    nu_o_nu0_sq = 1.+ a1_ii*f + (1./2.)*a2_iikk * pow(f,2.) # eq 41
    eta_s = - gamma - (1./2. * nu_o_nu0_sq * pow(2.*f+1.,2.)*a2_s) # eq 46
        
    #delta_Uqn = C_v* (T-300.)  #NEED TO KNOW POTENTIAL TEMPERATURE
    #delta_Uq = delta_Uqn
    #print "f	density	total_params['molar_mass']	eta_s	delta_u	p	T"
    #print f,density, params['molar_mass'],eta_s, delta_U, p, T
    
    delta_U = (P_th(V) - P_th_ref(V))*(V/gamma)
    G = pow(1.+2.*f, 5./2.) * (params['ref_mu'] + (3.*(params['ref_K']*params['mu_prime']) - 5.*params['ref_mu'])*f \
                                   + (6.*params['ref_K']*params['mu_prime']) - 24.*params['ref_K'] -14.*params['ref_mu'] + 9./2.*params['ref_K']*params['mu_prime'])*pow(f,2.) \
				   - (eta_s*(1./V)*delta_U/1.e9)
    return G


def stx_bulk_modulus(p,T,V,params):
    # simple model:
    #K = params['ref_K'] + params['K_prime'] * p + dKdT*(T-300.)
    #K = K * (1. + alpha * params['ref_grueneisen'] * T)        # formula Matas, D6
    func_int = lambda t: pow(t,3.)/(math.exp(t)-1.)
    
    P_th = lambda x: mgd.mgd_thermal_pressure(x, T, params)
    P_th_ref = lambda x: mgd.mgd_thermal_pressure(x, 300., params)
    
    
    func_int_cv = lambda t: math.exp(t)*pow(t,4.)/(math.exp(t)-1.)
    Dcv = integrate.quad(func_int_cv,0.,params['ref_Debye']/T) 
    Dcv_300 = integrate.quad(func_int_cv,0.,params['ref_Debye']/300.)

    cv_300 = 9.*params['n']*gas_constant*(pow(300./params['ref_Debye'],3.))*Dcv_300[0]
    Cv = (9.*params['n']*gas_constant*(pow(T/params['ref_Debye'],3.))*Dcv[0])-cv_300 #units of R

    #density_0 = (params['molar_mass'] / (params['ref_V']*1e-6) )
    #Kth_0 = params['ref_K'] - (pow(params['ref_grueneisen'],2.) * density_0 * delta_U_300/params['molar_mass'])/1.e9 
    f =.5*(pow(params['ref_V']/V,2./3.)-1)
    gamma = params['ref_grueneisen'] * (pow(params['ref_V']/V,-params['q0']))
    
    delta_U = (P_th(V) - P_th_ref(V))*(V/gamma) #convert into Stixrudian units

    K = pow(1.+2.*f, 5./2.) * ( params['ref_K'] + (3*params['ref_K']*params['K_prime'] -5.*params['ref_K'])*f+27./2.*(params['ref_K']*params['K_prime']-(4.*params['ref_K']))*pow(f,2.)) \
        + ((gamma+1.-params['q0'])*(gamma/V)  *delta_U )/1.e9 \
        - ((pow(gamma,2.) / V )*(Cv*T - cv_300*300.) )/1.e9   # eq 32

    return K


