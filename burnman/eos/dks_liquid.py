# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

from __future__ import absolute_import

from os import path
import numpy as np
import scipy.optimize as opt
try:
    # scipy's factorial was moved to special in scipy 1.3.0+
    from scipy.special import factorial
except ImportError:
    from scipy.misc import factorial

import warnings

from . import equation_of_state as eos
from .. import constants as constants
from ..utils.chemistry import read_masses
from ..utils.math import bracket

atomic_masses=read_masses()
# energy_states should provide the energies and degeneracies of each electronic level in a variety of elements


class DKS_L(eos.EquationOfState):
    """
    Base class for the finite strain liquid equation of state detailed
    in :cite:`deKoker2013` (supplementary materials).
    """

    """
    Ideal gas contributions (translational and electronic)
    to thermodynamic properties
    """

    def _ln_partition_function(self, mass, temperature):
        """
        Calculates the natural log of the partition function
        """
        return 3./2.*np.log(temperature) \
            + 3./2.*np.log(mass*constants.Boltzmann \
                           /(2*np.pi*constants.Dirac*constants.Dirac)) \

    def _F_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the helmholtz free energy
        Eq. S6, see also eq. 16.72 of Callen., 1985; p. 373
        """

        V = volume/constants.Avogadro
        figoverRT=0.
        for element, N in params['formula'].items(): # N is a.p.f.u
            if N > 1.e-5:
                mass = atomic_masses[element]/constants.Avogadro
                figoverRT += -N*(np.log(V) + self._ln_partition_function(mass, temperature) \
                                     + 1.) + N*np.log(N)
        return constants.gas_constant*temperature*figoverRT


    def _S_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the entropy
        """

        V = volume/constants.Avogadro
        entropy_sum=0.
        for element, N in params['formula'].items(): # N is a.p.f.u
            if N > 1.e-5:
                mass = atomic_masses[element]/constants.Avogadro
                entropy_sum -= -N*(np.log(V) + self._ln_partition_function(mass, temperature) \
                                     + 5./2.) + N*np.log(N)
        return constants.gas_constant*entropy_sum

    def _C_v_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the heat capacity
        """

        n_atoms=0
        for element, N in params['formula'].items():
            n_atoms += N
        return 1.5*constants.gas_constant*n_atoms

    def _P_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the pressure
        PV = nRT
        """

        n_atoms=0
        for element, N in params['formula'].items():
            n_atoms += N
        return n_atoms*constants.gas_constant*temperature / volume

    def _K_T_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the isothermal bulk modulus
        V * d/dV(-nRT/V) = V*nRT/V^2
        """
        n_atoms=0
        for element, N in params['formula'].items():
            n_atoms += N
        return n_atoms*constants.gas_constant*temperature / volume

    def _alphaK_T_ig(self, temperature, volume, params):
        """
        The ideal gas contribution to the product of the
        thermal expansivity and isothermal bulk modulus
        d/dT(nRT/V) = nR/V
        """

        n_atoms=0
        for element, N in params['formula'].items():
            n_atoms += N
        return n_atoms*constants.gas_constant / volume

    """
    Electronic contributions to thermodynamic properties
    """

    def _zeta(self, temperature, volume, params): # eq. S5a, beta in deKoker thesis (3.34)
        return params['zeta_0']*(np.power(volume/params['el_V_0'], params['xi']))

    def _dzetadV(self, temperature, volume, params):
        return params['zeta_0']*params['xi']*(np.power(volume/params['el_V_0'], params['xi']))/volume

    def _d2zetadV2(self, temperature, volume, params):
        return params['zeta_0'] \
            * params['xi'] * (params['xi'] - 1.) \
            * (np.power(volume/params['el_V_0'], params['xi'])) \
            / volume / volume

    def _Tel(self, temperature, volume, params): # eq. S5b
        return params['Tel_0']*(np.power(volume/params['el_V_0'], params['eta']))

    def _dTeldV(self, temperature, volume, params):
        return params['Tel_0'] * params['eta'] \
            * (np.power(volume/params['el_V_0'], params['eta'])) \
            / volume

    def _d2TeldV2(self, temperature, volume, params):
        return params['Tel_0'] \
            * params['eta'] * (params['eta'] - 1.) \
            * (np.power(volume/params['el_V_0'], params['eta'])) \
            / volume / volume

    def _gimel(self, temperature_el, temperature, volume, params): # -F_el/zeta, 3.30 in de Koker thesis
        return 0.5*(temperature*temperature - temperature_el*temperature_el) \
            - temperature*temperature_el*np.log(temperature/temperature_el)

    def _dgimeldTel(self, temperature_el, temperature, volume, params):
        return (temperature-temperature_el) - temperature*np.log(temperature/temperature_el)

    def _dgimeldT(self, temperature_el, temperature, volume, params):
        return (temperature-temperature_el) - temperature_el*np.log(temperature/temperature_el)

    def _d2gimeldTdTel(self, temperature_el, temperature, volume, params):
        return -np.log(temperature/temperature_el)

    def _d2gimeldTel2(self, temperature_el, temperature, volume, params):
        return (temperature/temperature_el)  - 1.

    def _F_el(self, temperature, volume, params): # F_el
        temperature_el = self._Tel(temperature, volume, params)
        if temperature < temperature_el:
            F_el = 0
        else:
            F_el = -self._zeta(temperature, volume, params) \
                * self._gimel(temperature_el, temperature, volume, params)
        return F_el

    def _S_el(self, temperature, volume, params): # S_el
        temperature_el = self._Tel(temperature, volume, params)
        if temperature < temperature_el:
            S_el = 0
        else:
            S_el = self._zeta(temperature, volume, params) \
                * self._dgimeldT(temperature_el, temperature, volume, params)
        return S_el


    def _P_el(self, temperature, volume, params): # P_el
        temperature_el = self._Tel(temperature, volume, params)
        if temperature < temperature_el:
            P_el = 0
        else:
            P_el =  self._dzetadV(temperature, volume, params) \
                * self._gimel(temperature_el, temperature, volume, params) \
                + self._zeta(temperature, volume, params) \
                * self._dTeldV(temperature, volume, params) \
                * self._dgimeldTel(temperature_el, temperature, volume, params)
        return P_el

    def _K_T_el(self, temperature, volume, params): # K_T_el
        temperature_el = self._Tel(temperature, volume, params)
        if temperature < temperature_el:
            K_T_el = 0
        else:
            K_T_el =  -volume \
                * ( self._d2zetadV2(temperature, volume, params) \
                        * self._gimel(temperature_el, temperature, volume, params) \
                        + 2. * self._dzetadV(temperature, volume, params) \
                        * self._dgimeldTel(temperature_el, temperature, volume, params) \
                        * self._dTeldV(temperature, volume, params) \
                        + self._zeta(temperature, volume, params) \
                        * ( self._d2TeldV2(temperature, volume, params) \
                                * self._dgimeldTel(temperature_el, temperature, volume, params) \
                                + self._dTeldV(temperature, volume, params) \
                                * self._dTeldV(temperature, volume, params) \
                                * self._d2gimeldTel2(temperature_el, temperature, volume, params)))
        return K_T_el

    def _alphaK_T_el(self, temperature, volume, params): # (alphaK_T)_el
        temperature_el = self._Tel(temperature, volume, params)
        if temperature < temperature_el:
            alphaK_T_el = 0
        else:
            alphaK_T_el = self._dzetadV(temperature, volume, params) \
                * self._dgimeldT(temperature_el, temperature, volume, params) \
                + self._zeta(temperature, volume, params) \
                * self._d2gimeldTdTel(temperature_el, temperature, volume, params) \
                * self._dTeldV(temperature, volume, params)
        return alphaK_T_el

    def _C_v_el(self, temperature, volume, params): # C_el, eq. 3.28 of de Koker thesis
        temperature_el = self._Tel(temperature, volume, params)
        zeta = self._zeta(temperature, volume, params)

        if temperature > temperature_el:
            Cv_el = zeta*(temperature - temperature_el)
        else:
            Cv_el = 0.
        return Cv_el



    """
    Excess (bonding) contributions to thermodynamic properties
    """

    # Finite strain
    def _finite_strain(self, temperature, volume, params): # f(V), eq. S3a
        return (1./2.)*(np.power(params['V_0']/volume, 2./3.) - 1.0)

    def _dfdV(self, temperature, volume, params): # f(V), eq. S3a
        return (-1./3.)*np.power(params['V_0']/volume, 2./3.)/volume

    def _d2fdV2(self,temperature, volume, params):
        return (5./9.)*np.power(params['V_0']/volume, 2./3.)/volume/volume

    # Temperature
    def _theta(self, temperature, volume, params): # theta, eq. S3b
        return np.power(temperature/params['T_0'], params['m']) - 1.

    def _dthetadT(self, temperature, volume, params):
        return params['m']*np.power(temperature/params['T_0'], params['m']) \
            / temperature

    def _d2thetadT2(self, temperature, volume, params):
        return params['m']*(params['m']-1.)*np.power(temperature/params['T_0'], params['m']) \
            / temperature / temperature

    def _F_xs(self, temperature, volume, params): # F_xs, eq. S2
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        energy = 0.
        for i in range(len(params['a'])):
            ifact=factorial(i, exact=False)
            for j in range(len(params['a'][0])):
                jfact=factorial(j, exact=False)
                energy += params['a'][i][j]*np.power(f, i)*np.power(theta, j)/ifact/jfact
        return energy

    def _S_xs(self, temperature, volume, params): # F_xs, eq. 3.18
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        entropy = 0.
        for i in range(len(params['a'])):
            ifact = factorial(i, exact=False)
            for j in range(len(params['a'][0])):
                if j > 0:
                    jfact = factorial(j, exact=False)
                    entropy += j*params['a'][i][j]*np.power(f, i)*np.power(theta, j-1.)/ifact/jfact
        return -self._dthetadT(temperature, volume, params)*entropy

    def _P_xs(self, temperature, volume, params): # P_xs, eq. 3.17 of de Koker thesis
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        pressure=0.
        for i in range(len(params['a'])):
            ifact=factorial(i, exact=False)
            if i > 0:
                for j in range(len(params['a'][0])):
                    jfact=factorial(j, exact=False)
                    pressure += float(i)*params['a'][i][j]*np.power(f, float(i)-1.)*np.power(theta, float(j))/ifact/jfact
        return -self._dfdV(temperature, volume, params)*pressure

    def _K_T_xs(self, temperature, volume, params): # K_T_xs, eq. 3.20 of de Koker thesis
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        K_ToverV=0.
        for i in range(len(params['a'])):
            ifact=factorial(i, exact=False)
            for j in range(len(params['a'][0])):
                if i > 0:
                    jfact=factorial(j, exact=False)
                    prefactor = float(i) * params['a'][i][j] \
                        * np.power(theta, float(j)) / ifact / jfact
                    K_ToverV += prefactor*self._d2fdV2(temperature, volume, params) \
                        * np.power(f, float(i-1))
                if i > 1:
                    dfdV = self._dfdV(temperature, volume, params)
                    K_ToverV += prefactor * dfdV * dfdV \
                        * float(i-1) * np.power(f, float(i-2))
        return volume*K_ToverV

    def _alphaK_T_xs(self, temperature, volume, params): # eq. 3.21 of de Koker thesis
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        sum_factors = 0.
        for i in range(len(params['a'])):
            ifact=factorial(i, exact=False)
            if i > 0:
                for j in range(len(params['a'][0])):
                    if j > 0:
                        jfact=factorial(j, exact=False)
                        sum_factors += float(i)*float(j)*params['a'][i][j] \
                            * np.power(f, float(i-1)) * np.power(theta, float(j-1)) \
                            / ifact / jfact

        return -self._dfdV(temperature, volume, params) \
            * self._dthetadT(temperature, volume, params) \
            * sum_factors


    def _C_v_xs(self, temperature, volume, params): # Cv_xs, eq. 3.22 of de Koker thesis
        f = self._finite_strain(temperature, volume, params)
        theta = self._theta(temperature, volume, params)
        C_voverT=0.
        for i in range(len(params['a'])):
            ifact=factorial(i, exact=False)
            for j in range(len(params['a'][0])):
                if j > 0:
                    jfact=factorial(j, exact=False)
                    prefactor = float(j)*params['a'][i][j]*np.power(f, float(i))/ifact/jfact
                    C_voverT += prefactor * self._d2thetadT2(temperature, volume, params) \
                        * np.power(theta, float(j-1))
                if j > 1:
                    dthetadT = self._dthetadT(temperature, volume, params)
                    C_voverT += prefactor * dthetadT * dthetadT \
                        * float(j-1) * np.power(theta, float(j-2))
        return -temperature*C_voverT


    """
    Magnetic contributions to thermodynamic properties
    (as found in Ramo and Stixrude, 2014)
    """

    def _spin(self, temperature, volume, params):
        S_a = 0.
        S_b = 0.
        numerator = 0.
        numerator_2 = 0.
        n_atoms = 0.
        if 'spin_a' in params:
            for element, N in params['formula'].items():
                if element == 'Fe':
                    n_atoms += N

            VoverVx = volume/params['V_0']
            S_a = params['spin_a'][0] + params['spin_a'][1]*VoverVx
            S_b = (params['spin_b'][0]
                   + params['spin_b'][1]/VoverVx
                   + params['spin_b'][2]/(np.power(VoverVx, 2.))
                   + params['spin_b'][3]/(np.power(VoverVx, 3.)))

            # S = S_a*T + S_b
            # d(2S + 1)/dV
            numerator=-2.*(-params['spin_a'][1]*temperature
                           + params['spin_b'][1]/(np.power(VoverVx, 2.))
                           + 2.*params['spin_b'][2]/(np.power(VoverVx, 3.))
                           + 3.*params['spin_b'][3]/(np.power(VoverVx, 4.)))/params['V_0']

            # d2S/dV2
            numerator_2 = 2.*((2.*params['spin_b'][1]/(np.power(VoverVx, 3.))
                               + 6.*params['spin_b'][2]/(np.power(VoverVx, 4.))
                               + 12.*params['spin_b'][3]/(np.power(VoverVx, 5.)))
                              /np.power(params['V_0'], 2.))
        return S_a, S_b, numerator, numerator_2, n_atoms


    def _F_mag(self, temperature, volume, params):
        S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
        S = S_a*temperature + S_b
        return -n_atoms*constants.gas_constant*temperature*np.log(2.*S + 1.)


    def _S_mag(self, temperature, volume, params):
        S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
        S = S_a*temperature + S_b
        return n_atoms*constants.gas_constant * ((2.*S_a*temperature/(2.*S + 1.)
                                                   + np.log(2.*S + 1.)))

    def _P_mag(self, temperature, volume, params):
        S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
        S = S_a*temperature + S_b
        dFdV = -n_atoms*constants.gas_constant*temperature*numerator/(2.*S + 1.)
        return -dFdV


    def _K_T_mag(self, temperature, volume, params):
        S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
        S = S_a*temperature + S_b
        dFdV = numerator/(2.*S + 1.)
        d2FdV2 = numerator_2/(2.*S + 1.) - np.power(dFdV, 2.)

        return -volume*n_atoms*constants.gas_constant*temperature*d2FdV2

    def _alphaK_T_mag(self, temperature, volume, params): # WARNING: numeric differentiation a.t.m.
        return (self._P_mag(temperature + 0.5, volume, params)
                - self._P_mag(temperature - 0.5, volume, params))


    def _C_v_mag(self, temperature, volume, params):
        S_a, S_b, numerator, numerator_2, n_atoms = self._spin(temperature, volume, params)
        S = S_a*temperature + S_b
        return n_atoms * constants.gas_constant * temperature * 4.*S_a*(S_a*temperature + 2.*S_b + 1.)/np.power(2.*S + 1., 2.)


    def _aK_T(self, temperature, volume, params):
        aK_T =  (self._alphaK_T_ig(temperature, volume, params)
                 + self._alphaK_T_el(temperature, volume, params)
                 + self._alphaK_T_xs(temperature, volume, params)
                 + self._alphaK_T_mag(temperature, volume, params))
        return aK_T

    # Pressure
    def pressure(self, temperature, volume, params):
        P = (self._P_ig(temperature, volume, params)
             + self._P_el(temperature, volume, params)
             + self._P_xs(temperature, volume, params)
             + self._P_mag(temperature, volume, params))
        return P

    def volume(self, pressure, temperature, params):
        _delta_pressure = lambda x, pressure, temperature, params: pressure - self.pressure(temperature, x, params)

        # we need to have a sign change in [a,b] to find a zero. Let us start with a
        # conservative guess:
        args = (pressure, temperature, params)
        try:
            sol = bracket(_delta_pressure, params['V_0'],
                          1.e-2 * params['V_0'], args)
        except ValueError:
            raise Exception(
                'Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?')
        return opt.brentq(_delta_pressure, sol[0], sol[1], args=args)

    def isothermal_bulk_modulus(self, pressure,temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]`
        """
        K_T = (self._K_T_ig(temperature, volume, params)
               + self._K_T_el(temperature, volume, params)
               + self._K_T_xs(temperature, volume, params)
               + self._K_T_mag(temperature, volume, params))
        return K_T

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus. :math:`[Pa]`
        """
        K_S = (self.isothermal_bulk_modulus(pressure,temperature, volume, params)
               * ( 1. + temperature
                   * self.thermal_expansivity(pressure, temperature, volume, params)
                   * self.grueneisen_parameter(pressure, temperature, volume, params)))
        return K_S

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter. :math:`[unitless]`
        """
        gamma = (self._aK_T(temperature, volume, params)
                 * volume
                 / self.molar_heat_capacity_v(pressure, temperature, volume, params))
        return gamma

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus. :math:`[Pa]`
        Zero for fluids
        """
        return 0.

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]`
        """
        C_v = (self._C_v_ig(temperature, volume, params)
               + self._C_v_el(temperature, volume, params)
               + self._C_v_xs(temperature, volume, params)
               + self._C_v_mag(temperature, volume, params))
        return C_v

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]`
        """
        C_p = (self.molar_heat_capacity_v(pressure,temperature, volume, params)
               * ( 1. + temperature
                   * self.thermal_expansivity(pressure, temperature, volume, params)
                   * self.grueneisen_parameter(pressure, temperature, volume, params) ))
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        """
        alpha = (self._aK_T(temperature, volume, params)
                 / self.isothermal_bulk_modulus(0., temperature, volume, params))
        return alpha

    def gibbs_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        """
        G = self.helmholtz_free_energy( pressure, temperature, volume, params) + pressure * volume
        return G

    def entropy( self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """
        S = (self._S_ig(temperature, volume, params)
             + self._S_el(temperature, volume, params)
             + self._S_xs(temperature, volume, params)
             + self._S_mag(temperature, volume, params))
        return S

    def enthalpy( self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        """
        H = self.helmholtz_free_energy( pressure, temperature, volume, params) + \
            temperature * self.entropy( pressure, temperature, volume, params) + \
            pressure * self.volume( pressure, temperature, params)
        return H

    def helmholtz_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        """
        F = (self._F_ig(temperature, volume, params)
             + self._F_el(temperature, volume, params)
             + self._F_xs(temperature, volume, params)
             + self._F_mag(temperature, volume, params))
        return F

    def molar_internal_energy(self, pressure, temperature, volume, params):
        E = self.helmholtz_free_energy(pressure, temperature, volume, params) + \
            temperature*self.entropy(pressure, temperature, volume, params)
        return E

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Check that all the required keys are in the dictionary
        expected_keys = ['V_0', 'T_0', 'O_theta', 'O_f', 'm', 'a', 'zeta_0', 'xi', 'Tel_0', 'eta']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Sometimes the standard electronic volume is different to V_0.
        # If not, make it the same.
        if 'el_V_0' not in params:
            params['el_V_0'] = params['V_0']
