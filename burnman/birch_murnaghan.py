# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import scipy.optimize as opt

def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in Pa.
    """

    x = params['ref_V']/volume
    f = 0.5*(pow(x, 2./3.) - 1.0)

    K = pow(1. + 2*f, 5./2.)* (params['ref_K'] + (3. * params['ref_K'] * params['K_prime'] -  \
           5*params['ref_K'] ) * f + 27./2. * (params['ref_K']*params['K_prime'] - 4* params['ref_K'])*f*f)
    return K

def birch_murnaghan(x, params):
    """
    equation for the third order birch-murnaghan equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['ref_K'])
    """

    return 3.*params['ref_K']/2. * (pow(x, 7./3.) - pow(x, 5./3.)) \
    * (1 - .75*(4-params['K_prime'] )*(pow(x, 2./3.) - 1))

def density(pressure, params):
    """ 
    Get the birch-murnaghan density at a reference temperature for a given
    pressure (in Pa). Returns density in kg/m^3
    """
    return params['molar_mass']/volume(pressure,params)

def volume(pressure, params):
    """
    Get the birch-murnaghan volume at a reference temperature for a given
    pressure (Pa). Returns molar volume in m^3
    """

    func = lambda x: birch_murnaghan(params['ref_V']/x, params) - pressure
    V = opt.brentq(func, 0.5*params['ref_V'], 1.5*params['ref_V'])
    return V

def shear_modulus_second_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in Pa (the same units as in
    params['ref_mu']).  This uses a second order finite strain expansion
    """

    x = params['ref_V']/volume
    G=params['ref_mu'] * pow(x,5./3.)*(1.-0.5*(pow(x,2./3.)-1.)*(5.-3.*params['mu_prime']*params['ref_K']/params['ref_mu']))
    return G 

def shear_modulus_third_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in Pa (the same units as in
    params['ref_mu']).  This uses a third order finite strain expansion
    """

    x = params['ref_V']/volume
    f = 0.5*(pow(x, 2./3.) - 1.0)
    G = pow((1. + 2*f), 5./2.)*(params['ref_mu']+(3.*params['ref_K']*params['mu_prime'] - 5.*params['ref_mu'])*f + (6.*params['ref_K']*params['mu_prime']-24.*params['ref_K']-14.*params['ref_mu']+9./2. * params['ref_K']*params['K_prime'])*f*f)
    return G 

def shear_modulus(volume, params):
    """
    More mineral physics studies seem to fit their shear modulus data to
    the second order finite strain expansion (can we confirm this more?)
    so the default shear_modulus function wraps the second order function.
    """
    return shear_modulus_second_order(volume, params)
