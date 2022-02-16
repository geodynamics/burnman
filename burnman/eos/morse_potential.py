from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
from ..utils.math import bracket
import warnings
import numpy as np

def bulk_modulus(volume, params):
    """
    Compute the bulk modulus as per the Morse potential
    equation of state.
    Returns bulk modulus in the same units as
    the reference bulk modulus.
    Pressure must be in :math:`[Pa]`.
    """

    VoverV0 = volume / params['V_0']
    x = (params['Kprime_0']  - 1.)*(1. - np.power(VoverV0, 1./3.))
    K = params['K_0']*( ( 2./(params['Kprime_0']  - 1.) *
                          np.power(VoverV0, -2./3.) *
                          (np.exp(2.*x) - np.exp(x)) ) +
                        ( np.power(VoverV0, -1./3.) *
                          (2.*np.exp(2.*x) - np.exp(x)) ) )
    return K

def shear_modulus(volume, params):
    """
    Shear modulus not currently implemented for this equation of state
    """
    return 0.

def morse_potential(VoverV0, params):
    """
    Equation for the Morse Potential equation of state,
    returns pressure in the same units that are supplied
    for the reference bulk modulus (params['K_0'])
    """
    x = (params['Kprime_0']  - 1.)*(1. - np.power(VoverV0, 1./3.))
    return ( 3. * params['K_0'] / (params['Kprime_0']  - 1.) *
             np.power(VoverV0, -2./3.) *
             (np.exp(2.*x) - np.exp(x)) ) + params['P_0']

def volume(pressure, params):
    """
    Get the Morse Potential volume at a
    reference temperature for a given pressure :math:`[Pa]`.
    Returns molar volume in :math:`[m^3]`
    """
    func = lambda V: morse_potential(V / params['V_0'], params) - pressure
    try:
        sol = bracket(func, params['V_0'], 1.e-2 * params['V_0'])
    except:
        raise ValueError(
            'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
    return opt.brentq(func, sol[0], sol[1])



class Morse(eos.EquationOfState):

    """
    Class for the isothermal Morse Potential equation of state
    detailed in :cite:`Stacey1981`.
    This equation of state has no temperature dependence.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume(pressure, params)

    def pressure(self, temperature, volume, params):
        return morse_potential(volume / params['V_0'], params)

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`K_T` :math:`[Pa]` as a function of pressure :math:`[Pa]`,
        temperature :math:`[K]` and volume :math:`[m^3]`.
        """
        return bulk_modulus(volume, params)

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus :math:`K_s` of the mineral. :math:`[Pa]`.
        """
        return bulk_modulus(volume, params)

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus :math:`G` of the mineral. :math:`[Pa]`
        """
        return shear_modulus(volume, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        """
        return 0.

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\mathcal{E}` of the mineral. :math:`[J/mol]`
        """

        x = (params['Kprime_0'] - 1)*(1 - np.power(volume/params['V_0'], 1./3.))
        intPdV = ( 9./2. * params['V_0'] * params['K_0'] /
                   np.power(params['Kprime_0'] - 1., 2.) *
                   (2.*np.exp(x) - np.exp(2.*x) - 1.) )

        return -intPdV + params['E_0']

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        return self.molar_internal_energy(pressure, temperature, volume, params) + volume*pressure

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
        Check for existence and validity of the parameters
        """

        if 'E_0' not in params:
            params['E_0'] = 0.
        if 'P_0' not in params:
            params['P_0'] = 0.

        # If G and Gprime are not included this is presumably deliberate,
        # as we can model density and bulk modulus just fine without them,
        # so just add them to the dictionary as nans
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')

        # Check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'K_0', 'Kprime_0', 'G_0', 'Gprime_0']
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
        if params['G_0'] < 0.0 or params['G_0'] > 1.e13:
            warnings.warn('Unusual value for G_0', stacklevel=2)
        if params['Gprime_0'] < -5. or params['Gprime_0'] > 10.:
            warnings.warn('Unusual value for Gprime_0', stacklevel=2)
