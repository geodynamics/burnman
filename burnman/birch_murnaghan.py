# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

import scipy.optimize as opt
import burnman.equation_of_state as eos

def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in Pa.
    """

    x = params['V_0']/volume
    f = 0.5*(pow(x, 2./3.) - 1.0)

    K = pow(1. + 2.*f, 5./2.)* (params['K_0'] + (3. * params['K_0'] * params['Kprime_0'] -  \
           5*params['K_0'] ) * f + 27./2. * (params['K_0']*params['Kprime_0'] - 4.* params['K_0'])*f*f)
    return K

def birch_murnaghan(x, params):
    """
    equation for the third order birch-murnaghan equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    """

    return 3.*params['K_0']/2. * (pow(x, 7./3.) - pow(x, 5./3.)) \
    * (1 - .75*(4-params['Kprime_0'] )*(pow(x, 2./3.) - 1))

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

    func = lambda x: birch_murnaghan(params['V_0']/x, params) - pressure
    V = opt.brentq(func, 0.5*params['V_0'], 1.5*params['V_0'])
    return V

def shear_modulus_second_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in Pa (the same units as in
    params['G_0']).  This uses a second order finite strain expansion
    """

    x = params['V_0']/volume
    G=params['G_0'] * pow(x,5./3.)*(1.-0.5*(pow(x,2./3.)-1.)*(5.-3.*params['Gprime_0']*params['K_0']/params['G_0']))
    return G

def shear_modulus_third_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in Pa (the same units as in
    params['G_0']).  This uses a third order finite strain expansion
    """

    x = params['V_0']/volume
    f = 0.5*(pow(x, 2./3.) - 1.0)
    G = pow((1. + 2*f), 5./2.)*(params['G_0']+(3.*params['K_0']*params['Gprime_0'] - 5.*params['G_0'])*f + (6.*params['K_0']*params['Gprime_0']-24.*params['K_0']-14.*params['G_0']+9./2. * params['K_0']*params['Kprime_0'])*f*f)
    return G


class BirchMurnaghanBase(eos.EquationOfState):
    """
    Base class for the isothermal Birch Murnaghan equation of state.  This is third order in strain, and
    has no temperature dependence.  However, the shear modulus is sometimes fit to a second order 
    function, so if this is the case, you should use that.  For more see :class:`burnman.birch_murnaghan.BM2` and :class:`burnman.birch_murnaghan.BM3`.
    """
    def volume(self,pressure, temperature, params):
        return volume(pressure,params)

    def isothermal_bulk_modulus(self,pressure,temperature, volume, params):
        return bulk_modulus(volume, params)

    def adiabatic_bulk_modulus(self,pressure, temperature, volume, params):
        return bulk_modulus(volume,params)

    def shear_modulus(self,pressure, temperature, volume, params):
        if(self.order == 2):
          return shear_modulus_second_order(volume,params)
        elif(self.order == 3):
          return shear_modulus_third_order(volume,params)

    def heat_capacity_v(self,pressure, temperature, volume, params):
        """
        Since this equation of state does not have temperature effects, simply return a very large number.
        """
        return 1.e99

    def heat_capacity_p(self,pressure, temperature, volume, params):
        """
        Since this equation of state does not have temperature effects, simply return a very large number.
        """
        return 1.e99

    def thermal_expansivity(self,pressure, temperature, volume, params):
        """
        Since this equation of state does not have temperature effects, simply return zero.
        """
        return 0.

    def grueneisen_parameter(self,pressure,temperature,volume,params):
        """
        Since this equation of state does not have temperature effects, simply return zero.
        """
        return 0.


class BM3(BirchMurnaghanBase):
    """
    Third order Birch Murnaghan isothermal equation of state.  
    This uses the third order expansion for shear modulus.
    """
    def __init__(self):
        self.order=3


class BM2(BirchMurnaghanBase):
    """
    Third order Birch Murnaghan isothermal equation of state.  
    This uses the third order expansion for shear modulus.
    """
    def __init__(self):
        self.order=2
