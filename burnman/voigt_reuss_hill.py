# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import seismic

def voigt_reuss_hill(pressure, temperature, phases, molar_abundances):
    
    """ Compute Voigt-Reuss-Hill average.
    
    Inputs: 
    pressure: [Pa]
    temperature:
    phases: list of n objects from type material
    molar_abundances: list of n molar fractions (sums up to 1.0)
     
    Returns: rho, v_p, v_s, v_phi, K_vrh, mu_vrh (all scalars)    
    """
    
    # This is based on Matas 2007, Appendix C

    n_phase = len(phases)
    assert n_phase == len(molar_abundances)
    assert sum(molar_abundances)  > 0.999
    assert sum(molar_abundances)  < 1.001
    it = range(n_phase)
    
    #determine the partial molar volumes of the phases
    for i in it:
        phases[i].set_state(pressure, temperature)

    V_i = [(molar_abundances[i]*phases[i].molar_volume()) for i in it]
    V_tot= sum(V_i)
    #calculate the density of the phase assemblage
    rho = (1./V_tot) * sum( molar_abundances[i]*phases[i].molar_mass() for i in it)

    #calculate the voigt and reuss averages of the phase assemblage for K and mu
    K_i = [phases[i].adiabatic_bulk_modulus() for i in it]
    K_vrh = vhr_average(V_i,K_i)
    mu_i = [phases[i].shear_modulus() for i in it]
    mu_vrh = vhr_average(V_i,mu_i)
    
    #compute seismic velocities
    v_s = np.sqrt( mu_vrh / rho)
    v_p = np.sqrt( (K_vrh + 4./3.*mu_vrh) / rho)
    v_phi = np.sqrt( (K_vrh) / rho)

    return rho, v_p, v_s, v_phi, K_vrh, mu_vrh

def vhr_average(phase_volume, X):

    """ calculate the voigt and reuss averages of the phase assemblage for 
    a certain property X.
    
    Input: phase_volume: array of n volumes, X: array of n properties
    Returns: mixture of property X
    
    Source: Matas 2007, Appendix D """

    it = range(len(phase_volume))
    
    V_i = phase_volume
    V_tot = sum(V_i)

    X_voigt = sum( V_i[i]/V_tot * X[i] for i in it)
    X_reuss = 1./sum( V_i[i]/V_tot / X[i] for i in it)
    X_vrh = (X_voigt+X_reuss)/2.
    return X_vrh

def attenuation_correction(v_p,v_s,v_phi,Qs,Qphi):
  
    beta = 0.3 # Matas et al. (2007) page 4        
    Qp  = 3./4.*pow((v_p/v_s),2.)*Qs    # Matas et al. (2007) page 4


    cot=1./np.tan(beta*np.pi/2.)
    v_p  = v_p*(1.-1./2.*cot*1./Qp)    # Matas et al. (2007) page 1
    v_s  = v_s*(1.-1./2.*cot*1./Qs)
    v_phi= v_phi*(1.-1./2.*cot*1./Qphi)
    return v_p, v_s, v_phi
