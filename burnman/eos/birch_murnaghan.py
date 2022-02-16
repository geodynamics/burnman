from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import scipy.optimize as opt
from . import equation_of_state as eos
from ..utils.math import bracket
import warnings


def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in :math:`[Pa]`.
    """

    x = params['V_0'] / volume
    f = 0.5 * (pow(x, 2. / 3.) - 1.0)

    K = pow(1. + 2. * f, 5. / 2.) * (params['K_0'] + (3. * params['K_0'] * params['Kprime_0'] -
                                                      5 * params['K_0']) * f + 27. / 2. * (params['K_0'] * params['Kprime_0'] - 4. * params['K_0']) * f * f)
    return K


def birch_murnaghan(x, params):
    """
    equation for the third order birch-murnaghan equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0'])
    """

    return 3. * params['K_0'] / 2. * (pow(x, 7. / 3.) - pow(x, 5. / 3.)) \
        * (1. - .75 * (4. - params['Kprime_0']) * (pow(x, 2. / 3.) - 1.)) + params['P_0']


def volume(pressure, params):
    """
    Get the birch-murnaghan volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    """

    func = lambda x: birch_murnaghan(params['V_0'] / x, params) - pressure
    try:
        sol = bracket(func, params['V_0'], 1.e-2 * params['V_0'])
    except:
        raise ValueError(
            'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
    return opt.brentq(func, sol[0], sol[1])


def shear_modulus_second_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in :math:`[Pa]` (the same units as in
    params['G_0']).  This uses a second order finite strain expansion
    """

    x = params['V_0'] / volume
    G = params['G_0'] * pow(x, 5. / 3.) * (
        1. - 0.5 * (pow(x, 2. / 3.) - 1.) * (5. - 3. * params['Gprime_0'] * params['K_0'] / params['G_0']))
    return G


def shear_modulus_third_order(volume, params):
    """
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume.  Returns shear modulus in :math:`[Pa]` (the same units as in
    params['G_0']).  This uses a third order finite strain expansion
    """

    x = params['V_0'] / volume
    f = 0.5 * (pow(x, 2. / 3.) - 1.0)
    G = pow((1. + 2. * f), 5. / 2.) * (params['G_0'] + (3. * params['K_0'] * params['Gprime_0'] - 5. * params['G_0']) * f + (
        6. * params['K_0'] * params['Gprime_0'] - 24. * params['K_0'] - 14. * params['G_0'] + 9. / 2. * params['K_0'] * params['Kprime_0']) * f * f)
    return G


class BirchMurnaghanBase(eos.EquationOfState):

    """
    Base class for the isothermal Birch Murnaghan equation of state.  This is third order in strain, and
    has no temperature dependence.  However, the shear modulus is sometimes fit to a second order
    function, so if this is the case, you should use that.  For more see :class:`burnman.birch_murnaghan.BM2` and :class:`burnman.birch_murnaghan.BM3`.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume(pressure, params)

    def pressure(self, temperature, volume, params):
        return birch_murnaghan(params['V_0'] / volume, params)

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
        if(self.order == 2):
            return shear_modulus_second_order(volume, params)
        elif(self.order == 3):
            return shear_modulus_third_order(volume, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        """
        return 0.

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        x = np.power(volume/params['V_0'], -1./3.)
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        x8 = x4*x4

        xi1 = 3.*(4. - params['Kprime_0'])/4.

        intPdV = (-9./2. * params['V_0'] * params['K_0'] *
                  ((xi1 + 1.)*(x4/4. - x2/2. + 1./4.) -
                   xi1*(x6/6. - x4/4. + 1./12.)))

        return - intPdV + params['E_0']

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy :math:`\mathcal{G}` of the mineral. :math:`[J/mol]`
        """
        # G = int VdP = [PV] - int PdV = E + PV

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


class BM3(BirchMurnaghanBase):

    """
    Third order Birch Murnaghan isothermal equation of state.
    This uses the third order expansion for shear modulus.
    """

    def __init__(self):
        self.order = 3


class BM2(BirchMurnaghanBase):

    """
    Third order Birch Murnaghan isothermal equation of state.
    This uses the second order expansion for shear modulus.
    """

    def __init__(self):
        self.order = 2
