# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import numpy as np
import scipy.optimize as opt
import scipy.integrate as integrate
import matplotlib.pylab as plt
import mie_grueneisen_debye as mgd

def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in Pa.
    """

    x = params['ref_V']/volume
    pressure = birch_murnaghan(x, params)
    ## Matas et al. equation A2
    A= -2.*(params['K_prime']-4.)+8./3.
    B=-4.*(params['K_prime']-4.)+8./3.
    C= B-A
    top=5./3.*A*pow(x,5./3.)-7./3.*B*pow(x,7./3.)+3.*C*pow(x,3.)
    bottom=A*pow(x,5./3.)-B*pow(x,7./3.)+C*pow(x,3.)
    if bottom == 0.0 :
        return params['ref_K']
    K=pressure*top/bottom
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

    ratio = opt.brentq(lambda x: birch_murnaghan(x, params)-pressure, 0.8, 5.0)
    return ratio*params

def volume(pressure, params):
    """
    Get the birch-murnaghan volume at a reference temperature for a given
    pressure (Pa). Returns molar volume in m^3
    """

    ratio = opt.brentq(lambda x: birch_murnaghan(x, params)-pressure, 0.8, 5.0)
    return 1./ratio * params['ref_V']

def shear_modulus(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in Pa (the same units as in
    params['ref_mu'])
    """

    x = params['ref_V']/volume
    pressure = birch_murnaghan(x, params)
    G=params['ref_mu'] * pow(x,5./3.)*(1.-0.5*(pow(x,2./3.)-1.)*(5.-3.*params['mu_prime']*params['ref_K']/params['ref_mu']))
    f = 0.5*(pow(x, 2./3.) - 1.0)
    #G = pow((1. + 2*f), 5./2.)*(params['ref_mu']+(3.*params['ref_K']*params['mu_prime'] - 5.*params['ref_mu'])*f + (6.*params['ref_K']*params['mu_prime']-24.*params['ref_K']-14.*params['ref_mu']+9./2. * params['ref_K']*params['K_prime'])*f*f)
    return G 
