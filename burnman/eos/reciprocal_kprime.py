from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from scipy.special import gamma, gammainc
from . import equation_of_state as eos
from ..utils.math import bracket
import warnings
import numpy as np

# Try to import the jit from numba.  If it is
# not available, just go with the standard
# python interpreter
try:
    from numba import jit
except ImportError:
    def jit(fn):
        return fn


@jit
def _delta_PoverK_from_P(PoverK, pressure, K_0, Kprime_0, Kprime_inf):
    return PoverK - (pressure/K_0)*np.power((1. - Kprime_inf*PoverK), Kprime_0/Kprime_inf) # eq. 58

@jit
def _delta_PoverK_from_V(PoverK, V, V_0, K_0, Kprime_0, Kprime_inf):
    Kprime_ratio = Kprime_0 / Kprime_inf
    return ( np.log( V_0 / V ) +
             Kprime_ratio / Kprime_inf * np.log(1. - Kprime_inf * PoverK) +
             (Kprime_ratio - 1.) * PoverK ) # eq. 61

def _upper_incomplete_gamma(z, a):
    """
    An implementation of the non-regularised upper incomplete gamma
    function. Computed using the relationship with the regularised
    lower incomplete gamma function (scipy.special.gammainc).
    Uses the recurrence relation wherever z<0.
    """
    n = int(-np.floor(z))
    if n > 0:
        z = z + n
        u_gamma = (1. - gammainc(z, a))*gamma(z)

        for i in range(n):
            z = z - 1.
            u_gamma = (u_gamma - np.power(a, z)*np.exp(-a))/z
        return u_gamma
    else:
        return (1. - gammainc(z, a))*gamma(z)

def _PoverK_from_P(pressure, params):
    """
    Calculates the pressure:bulk modulus ratio
    from a given pressure using brentq optimization
    """
    args = ((pressure - params['P_0']), params['K_0'],
            params['Kprime_0'], params['Kprime_inf'])
    return opt.brentq(_delta_PoverK_from_P,
                      1./(params['Kprime_inf'] - params['Kprime_0']) + np.finfo(float).eps,
                      1./params['Kprime_inf'] - np.finfo(float).eps,
                      args=args)

def _PoverK_from_V(volume, params):
    """
    Calculates the pressure:bulk modulus ratio
    from a given volume using brentq optimization
    """
    args = (volume, params['V_0'], params['K_0'],
            params['Kprime_0'], params['Kprime_inf'])
    return opt.brentq(_delta_PoverK_from_V,
                      1./(params['Kprime_inf'] - params['Kprime_0']) + np.finfo(float).eps,
                      1./params['Kprime_inf'] - np.finfo(float).eps,
                      args=args)

def bulk_modulus(pressure, params):
    """
    Returns the bulk modulus at a given pressure
    """
    PoverK = _PoverK_from_P(pressure, params)
    K = params['K_0']*np.power((1. - params['Kprime_inf']*PoverK), -
                               params['Kprime_0']/params['Kprime_inf'])
    return K

def shear_modulus(pressure, params):
    """
    Returns the shear modulus at a given pressure
    """
    G = ( params['G_0']/params['K_0'] * bulk_modulus(pressure, params) -
          (params['G_0']/params['K_0']*params['Kprime_inf'] - params['Gprime_inf']) * pressure )
    return G # eq. 78


class RKprime(eos.EquationOfState):

    """
    Class for the isothermal reciprocal K-prime equation of state
    detailed in :cite:`StaceyDavis2004`.  This equation of state is
    a development of work by :cite:`Keane1954` and :cite:`Stacey2000`,
    making use of the fact that :math:`K'` typically varies smoothly
    as a function of :math:`P/K`, and is thermodynamically required to
    exceed 5/3 at infinite pressure.

    It is worth noting that this equation of state rapidly becomes
    unstable at negative pressures, so should not be trusted to provide
    a good *HT-LP* equation of state using a thermal pressure
    formulation. The negative root of :math:`dP/dK`
    can be found at :math:`K/P = K'_{\infty} - K'_0`,
    which corresponds to a bulk modulus of
    :math:`K = K_0 ( 1 - K'_{\infty}/K'_0 )^{K'_0/K'_{\infty}}`
    and a volume of
    :math:`V = V_0 ( K'_0 / (K'_0 - K'_{\infty}) )^{K'_0/{K'}^2_{\infty}} \exp{(-1/K'_{\infty})}`.

    This equation of state has no temperature dependence.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        Kprime_ratio = params['Kprime_0']/params['Kprime_inf']
        PoverK = _PoverK_from_P(pressure, params)

        V = params['V_0'] * np.exp( Kprime_ratio/params['Kprime_inf'] *
                                    np.log(1. - params['Kprime_inf'] * PoverK) +
                                    (Kprime_ratio - 1.) * PoverK ) # Eq. 61

        return V

    def pressure(self, temperature, volume, params):
        """
        Returns pressure :math:`[Pa]` as a function of volume :math:`[m^3]`.
        """
        PoverK = _PoverK_from_V(volume, params)
        return params['P_0'] + ( params['K_0'] * PoverK *
                                 np.power(1. - params['Kprime_inf'] * PoverK,
                                          -params['Kprime_0']/params['Kprime_inf']) )

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(pressure, params)

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus :math:`K_s` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(pressure, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return shear_modulus(pressure, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        """
        return 0.

    def _intVdP(self, xi, params):

        a = params['Kprime_inf']
        b = (params['Kprime_0']/params['Kprime_inf']/params['Kprime_inf'] -
             params['Kprime_0']/params['Kprime_inf'] - 1.)
        c = params['Kprime_0'] - params['Kprime_inf']
        f = (params['Kprime_0']/params['Kprime_inf'] - 1.)

        i1 = float( params['V_0'] * params['K_0'] *
                    np.exp(f / a) * np.power(a, b - 1.) /
                    np.power(f, b + 2.) *
                    ( f * params['Kprime_0'] * _upper_incomplete_gamma( b + 1. ,
                                                                        f * (1./a - xi) ) -
                      a * c * _upper_incomplete_gamma( b + 2., f * (1./a - xi) ) ) )

        return i1

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        # G = E0 + int VdP (when S = 0)
        K = self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        return params['E_0'] + params['P_0']*params['V_0'] + self._intVdP((pressure - params['P_0'])/K, params) - self._intVdP(0., params)

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        # E = G - PV (+ TS)
        return ( self.gibbs_free_energy(pressure, temperature, volume, params) - pressure*volume)

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return a very large number. :math:`[J/K/mol]`
        """
        return 1.e99

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[1/K]`
        """
        return 0.

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Since this equation of state does not contain temperature effects, simply return zero. :math:`[unitless]`
        """
        return 0.

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters.
        The value for :math:`K'_{\infty}` is thermodynamically bounded
        between 5/3 and :math:`K'_0` :cite:`StaceyDavis2004`.
        """

        if 'E_0' not in params:
            params['E_0'] = 0.
        if 'P_0' not in params:
            params['P_0'] = 0.

        # If G and Gprime_inf are not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_inf' not in params:
            params['Gprime_inf'] = float('nan')

        # Check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'K_0', 'Kprime_0', 'Kprime_inf', 'G_0', 'Gprime_inf']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
        if params['P_0'] < 0.:
            warnings.warn('Unusual value for P_0', stacklevel=2)
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-3:
            warnings.warn('Unusual value for V_0', stacklevel=2)
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn('Unusual value for K_0', stacklevel=2)
        if params['Kprime_0'] < 0. or params['Kprime_0'] > 10.:
            warnings.warn('Unusual value for Kprime_0', stacklevel=2)
        if params['Kprime_inf'] < 5./3. or params['Kprime_inf'] > params['Kprime_0']:
            warnings.warn('Unusual value for Kprime_inf', stacklevel=2) # eq. 17
        if params['G_0'] < 0.0 or params['G_0'] > 1.e13:
            warnings.warn('Unusual value for G_0', stacklevel=2)
        if params['Gprime_inf'] < -5. or params['Gprime_inf'] > 10.:
            warnings.warn('Unusual value for Gprime_inf', stacklevel=2)
