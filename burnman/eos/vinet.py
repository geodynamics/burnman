# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2017 by the BurnMan team, released under the GNU
# GPL v2 or later.


import scipy.optimize as opt
from . import equation_of_state as eos
import warnings
from math import exp


def bulk_modulus(volume, params):
    """
    compute the bulk modulus as per the
    Vinet equation of state.  Reference bulk
    modulus should be in :math:`[Pa]`.
    """

    x = volume / params['V_0']
    eta = (3. / 2.) * (params['Kprime_0'] - 1.)

    K = (params['K_0'] * pow(x, -2. / 3.)) * \
        (1 + ((eta * pow(x, 1. / 3.) + 1.) * (1. - pow(x, 1. / 3.)))) * \
        exp(eta * (1. - pow(x, 1. / 3.)))
    return K


def vinet(x, params):
    """
    equation for the  Vinet equation of state, returns
    pressure in the same units that are supplied for the reference bulk
    modulus (params['K_0']), which should be in math:`[Pa]`.
    """
    eta = (3. / 2.) * (params['Kprime_0'] - 1.)
    return 3. * params['K_0'] * (pow(x, -2. / 3.)) * (1. - (pow(x, 1. / 3.))) \
        * exp(eta * (1. - pow(x, 1. / 3.))) + params['P_0']


def volume(pressure, params):
    """
    Get the Vinet volume at a reference temperature for a given
    pressure :math:`[Pa]`. Returns molar volume in :math:`[m^3]`
    """

    func = lambda x: vinet(x / params['V_0'], params) - pressure
    V = opt.brentq(func, 0.1 * params['V_0'], 1.5 * params['V_0'])
    return V


class Vinet(eos.EquationOfState):

    """
    Base class for the isothermal Vinet equation of state.
    References for this equation of state are :cite:`vinet1986`
    and :cite:`vinet1987`. This equation of state actually
    predates Vinet by 55 years :cite:`Rydberg1932`,
    and was investigated further by :cite:`Stacey1981`.
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume :math:`[m^3]` as a function of pressure :math:`[Pa]`.
        """
        return volume(pressure, params)

    def pressure(self, temperature, volume, params):
        return vinet(volume / params['V_0'], params)

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
        Currently not included in the Vinet EOS, so omitted.
        """
        return 0.

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the molar entropy :math:`\mathcal{S}` of the mineral. :math:`[J/K/mol]`
        """
        return 0.

    def molar_internal_energy(self, pressure, temperature, volume, params):
        """
        Returns the internal energy :math:`\mathcal{E}` of the mineral. :math:`[J/mol]`
        """
        x = pow(volume/params['V_0'], 1./3.)
        eta = (3. / 2.) * (params['Kprime_0'] - 1.)


        intPdV = (9.* params['V_0'] * params['K_0'] / (eta*eta) *
                  ((1. - eta*(1. - x))*exp(eta*(1. - x)) - 1.))


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

        # G is not included in the Vinet EOS so we shall set them to NaN's
        if 'G_0' not in params:
            params['G_0'] = float('nan')
        if 'Gprime_0' not in params:
            params['Gprime_0'] = float('nan')

        # check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'K_0', 'Kprime_0']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # now check that the values are reasonable.  I mostly just
        # made up these values from experience, and we are only
        # raising a warning.  Better way to do this? [IR]
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-3:
            warnings.warn('Unusual value for V_0', stacklevel=2)
        if params['K_0'] < 1.e9 or params['K_0'] > 1.e13:
            warnings.warn('Unusual value for K_0', stacklevel=2)
        if params['Kprime_0'] < -5. or params['Kprime_0'] > 10.:
            warnings.warn('Unusual value for Kprime_0', stacklevel=2)
