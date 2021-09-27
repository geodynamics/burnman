from __future__ import absolute_import
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import numpy as np
import warnings

from . import modified_tait as mt
from . import murnaghan as murn
from . import equation_of_state as eos

from . import einstein


class HP_TMT(eos.EquationOfState):

    """
    Base class for the thermal equation of state based on
    the generic modified Tait equation of state (class MT),
    as described in :cite:`HP2011`.


    An instance "m" of a Mineral can be assigned this
    equation of state with the command m.set_method('hp_tmt')
    (or by initialising the class with the param
    equation_of_state = 'hp_tmt'
    """

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        EQ 12
        """
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.volume(pressure - Pth, params)

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        EQ B7
        """
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.modified_tait(params['V_0'] / volume, params) + Pth

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume.
        """
        alpha = self.thermal_expansivity(
            pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(
            pressure, temperature, volume, params)
        C_V = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        return alpha * K_T * volume / C_V

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].  EQ 13+2
        """
        Pth = self.__relative_thermal_pressure(temperature, params)
        return mt.bulk_modulus(pressure - Pth, params)

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.

    # Cv, heat capacity at constant volume
    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature,
        and volume [J/K/mol].
        """
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        V = self.volume(pressure, temperature, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(
            pressure, temperature, volume, params)
        return C_p - V * temperature * alpha * alpha * K_T

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]. This function replaces -Pth in EQ 13+1
        with P-Pth for non-ambient temperature
        """
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)
        psubpth = pressure - params['P_0'] - Pth

        C_V0 = einstein.molar_heat_capacity_v(
            params['T_0'], params['T_einstein'], params['n'])
        C_V = einstein.molar_heat_capacity_v(
            temperature, params['T_einstein'], params['n'])
        alpha = params['a_0'] * (C_V / C_V0) * 1. / (
            (1. + b * psubpth) * (a + (1. - a) * np.power((1 + b * psubpth), c)))
        return alpha

    def molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]. Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP2011`.
        """
        Cp = params['Cp'][0] + params['Cp'][1] * temperature + params['Cp'][2] * \
            np.power(temperature, -2.) + params[
                'Cp'][3] * np.power(temperature, -0.5)
        return Cp

    def molar_heat_capacity_p_einstein(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure at the pressure,
        temperature, and volume, using the C_v and Einstein model [J/K/mol]
        WARNING: Only for comparison with internally self-consistent C_p
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v * (1. + gr * alpha * temperature)
        return C_p

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        K_T = self.isothermal_bulk_modulus(
            pressure, temperature, volume, params)
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        K_S = K_T * C_p / C_v
        return K_S

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # Calculate temperature and pressure integrals
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)

        psubpth = pressure - params['P_0'] - Pth

        # EQ 13
        if pressure != params['P_0']:
            intVdP = (pressure - params['P_0']) * params['V_0'] * (
                1. - a + (a * (np.power((1. - b * Pth), 1. - c) - np.power((1. + b * (psubpth)), 1. - c)) / (b * (c - 1.) * (pressure - params['P_0']))))
        else:
            intVdP = 0.
        return params['H_0'] + self.__intCpdT(temperature, params) - temperature * (params['S_0'] + self.__intCpoverTdT(temperature, params)) + intVdP

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        return self.gibbs_free_energy(pressure, temperature, volume, params) - pressure * self.volume(pressure, temperature, params)

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = mt.tait_constants(params)
        Pth = self.__relative_thermal_pressure(temperature, params)

        ksi_over_ksi_0 = einstein.molar_heat_capacity_v(temperature, params['T_einstein'], params[
                                                  'n']) / einstein.molar_heat_capacity_v(params['T_0'], params['T_einstein'], params['n'])

        dintVdpdT = (params['V_0'] * params['a_0'] * params['K_0'] * a * ksi_over_ksi_0) * (
            np.power((1. + b * (pressure - params['P_0'] - Pth)), 0. - c) - np.power((1. - b * Pth), 0. - c))
        return params['S_0'] + self.__intCpoverTdT(temperature, params) + dintVdpdT

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        gibbs = self.gibbs_free_energy(pressure, temperature, volume, params)
        entropy = self.entropy(pressure, temperature, volume, params)
        return gibbs + temperature * entropy

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        a, b, c = mt.tait_constants(params)
        T = temperature
        T_e = params['T_einstein']
        n = params['n']
        Pth = self.__relative_thermal_pressure(T, params)

        ksi_over_ksi_0 = einstein.molar_heat_capacity_v(T, T_e, n) \
                         / einstein.molar_heat_capacity_v(params['T_0'], T_e, n)

        dintVdpdT = (params['V_0'] * params['a_0'] * params['K_0'] * a * ksi_over_ksi_0) * (
            np.power((1. + b * (pressure - params['P_0'] - Pth)), 0. - c) - np.power((1. - b * Pth), 0. - c))

        dSdT0 = params['V_0'] * params['K_0'] * np.power((ksi_over_ksi_0 * params['a_0']), 2.0) * \
                (np.power((1. + b * (pressure - params['P_0'] - Pth)), -1. - c)
                 - np.power((1. + b * (-Pth)), -1. - c))

        x = T_e/T
        dCv_einstdT = -(einstein.molar_heat_capacity_v(T, T_e, n)
                        * (1 - 2./x + 2./(np.exp(x) - 1.)) * x/T)

        dSdT1 = -dintVdpdT * dCv_einstdT \
                / einstein.molar_heat_capacity_v(T, T_e, n)

        dSdT = dSdT0 + dSdT1
        return self.molar_heat_capacity_p0(temperature, params) + temperature * dSdT

    def __thermal_pressure(self, T, params):
        """
        Returns thermal pressure [Pa] as a function of T [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """

        # This is basically the mie-gruneisen equation of state for thermal
        # pressure using an Einstein model for heat capacity.  The additional
        # assumption that they make is that alpha*K/Cv, (or gamma / V) is
        # constant over a wide range of compressions.

        # Note that the xi function in HP2011 is just the Einstein heat capacity
        # divided by 3nR. This function is *not* used to calculate the
        # heat capacity - Holland and Powell (2011) prefer the additional
        # freedom provided by their polynomial expression.

        E_th = einstein.thermal_energy(T, params['T_einstein'], params['n'])
        C_V0 = einstein.molar_heat_capacity_v(
            params['T_0'], params['T_einstein'], params['n'])
        P_th = params['a_0'] * params['K_0'] / C_V0 * E_th
        return P_th

    def __relative_thermal_pressure(self, T, params):
        """
        Returns relative thermal pressure [Pa] as a function of T-params['T_0'] [K]
        EQ 12 - 1 of :cite:`HP2011`.
        """
        return self.__thermal_pressure(T, params) - \
            self.__thermal_pressure(params['T_0'], params)

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return (params['Cp'][0] * temperature + 0.5 * params['Cp'][1] * np.power(temperature, 2.) - params['Cp'][2] / temperature + 2. * params['Cp'][3] * np.sqrt(temperature)) - (params['Cp'][0] * params['T_0'] + 0.5 * params['Cp'][1] * params['T_0'] * params['T_0'] - params['Cp'][2] / params['T_0'] + 2.0 * params['Cp'][3] * np.sqrt(params['T_0']))

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return (params['Cp'][0] * np.log(temperature) + params['Cp'][1] * temperature - 0.5 * params['Cp'][2] / np.power(temperature, 2.) - 2.0 * params['Cp'][3] / np.sqrt(temperature)) - (params['Cp'][0] * np.log(params['T_0']) + params['Cp'][1] * params['T_0'] - 0.5 * params['Cp'][2] / (params['T_0'] * params['T_0']) - 2.0 * params['Cp'][3] / np.sqrt(params['T_0']))

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if 'T_0' not in params:
            params['T_0'] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if 'H_0' not in params:
            params['H_0'] = float('nan')
        if 'S_0' not in params:
            params['S_0'] = float('nan')

        # First, let's check the EoS parameters for Tref
        mt.MT.validate_parameters(mt.MT(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ['H_0', 'S_0', 'V_0', 'Cp', 'a_0', 'n', 'molar_mass']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # The following line estimates the Einstein temperature
        # according to the empirical equation of
        # Holland and Powell, 2011; base of p.346, para.1
        if 'T_einstein' not in params:
            params['T_einstein'] = 10636. / \
                (params['S_0'] / params['n'] + 6.44)

        # Finally, check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn('Unusual value for T_0', stacklevel=2)
        if params['G_0'] is not float('nan') and (params['G_0'] < 0. or params['G_0'] > 1.e13):
            warnings.warn('Unusual value for G_0', stacklevel=2)
        if params['Gprime_0'] is not float('nan') and (params['Gprime_0'] < -5. or params['Gprime_0'] > 10.):
            warnings.warn('Unusual value for Gprime_0', stacklevel=2)

        # no test for H_0
        if params['S_0'] is not float('nan') and params['S_0'] < 0.:
            warnings.warn('Unusual value for S_0', stacklevel=2)
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-2:
            warnings.warn('Unusual value for V_0', stacklevel=2)

        if self.molar_heat_capacity_p0(params['T_0'], params) < 0.:
            warnings.warn('Negative heat capacity at T_0', stacklevel=2)
        if self.molar_heat_capacity_p0(2000., params) < 0.:
            warnings.warn('Negative heat capacity at 2000K', stacklevel=2)

        if params['a_0'] < 0. or params['a_0'] > 1.e-3:
            warnings.warn('Unusual value for a_0', stacklevel=2)

        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn('Unusual value for n', stacklevel=2)
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 10.:
            warnings.warn('Unusual value for molar_mass', stacklevel=2)


class HP_TMTL(eos.EquationOfState):

    """
    Base class for the thermal equation of state
    described in :cite:`HP1998`, but with the Modified Tait as the static part,
    as described in :cite:`HP2011`.

    An instance "m" of a Mineral can be assigned this
    equation of state with the command m.set_method('hp_tmtL')
    (or by initialising the class with the param
    equation_of_state = 'hp_tmtL'
    """

    def _V_T_1bar(self, temperature, params):
        # Constant thermal expansivity at standard state pressure
        # (p.348 of HP2011)
        return params['V_0']*np.exp(params['a_0']*(temperature - params['T_0']))

    def _K_T_1bar(self, temperature, params):
        # Linear bulk modulus dependence as in HP1998 (p.348 of HP2011)
        return params['K_0'] + params['dKdT_0']*(temperature - params['T_0'])

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        self.static_params['V_0'] = self._V_T_1bar(temperature, params)
        self.static_params['K_0'] = self._K_T_1bar(temperature, params)
        return mt.volume(pressure, self.static_params)

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        self.static_params['V_0'] = self._V_T_1bar(temperature, params)
        self.static_params['K_0'] = self._K_T_1bar(temperature, params)
        return mt.pressure(volume, self.static_params)

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume.
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        C_V = self.molar_heat_capacity_v(pressure, temperature,
                                         volume, params)
        return alpha * K_T * volume / C_V

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        self.static_params['V_0'] = self._V_T_1bar(temperature, params)
        self.static_params['K_0'] = self._K_T_1bar(temperature, params)
        return mt.bulk_modulus(pressure, self.static_params)

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.

    # Cv, heat capacity at constant volume
    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature,
        and volume [J/K/mol].
        """
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        V = self.volume(pressure, temperature, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        return C_p - V * temperature * alpha * alpha * K_T

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]
        """
        # The derivation of the high pressure thermal expansivity is tedious,
        # so here we take a numerical derivative.
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        self.static_params['V_0'] = self._V_T_1bar(temperature+dT/2., params)
        self.static_params['K_0'] = self._K_T_1bar(temperature+dT/2., params)
        volume1 = mt.volume(pressure, self.static_params)
        self.static_params['V_0'] = self._V_T_1bar(temperature-dT/2., params)
        self.static_params['K_0'] = self._K_T_1bar(temperature-dT/2., params)
        volume0 = mt.volume(pressure, self.static_params)

        return 2.*(volume1 - volume0)/(volume1 + volume0)/dT

    def molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP1998`.
        """
        Cp = (params['Cp'][0]
              + params['Cp'][1] * temperature
              + params['Cp'][2] * np.power(temperature, -2.)
              + params['Cp'][3] * np.power(temperature, -0.5))
        return Cp

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        K_S = K_T * C_p / C_v
        return K_S

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        self.static_params['V_0'] = self._V_T_1bar(temperature, params)
        self.static_params['K_0'] = self._K_T_1bar(temperature, params)
        return (params['H_0'] + self.__intCpdT(temperature, params)
                - temperature * (params['S_0']
                                 + self.__intCpoverTdT(temperature, params))
                + mt.intVdP(pressure, self.static_params))

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        return (self.gibbs_free_energy(pressure, temperature, volume, params)
                - pressure * self.volume(pressure, temperature, params))

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The derivation of the entropy is tedious,
        # so here we take a numerical derivative.
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        G1 = self.gibbs_free_energy(pressure, temperature+dT/2., volume, params)
        G0 = self.gibbs_free_energy(pressure, temperature-dT/2., volume, params)

        return (G0 - G1)/dT

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        gibbs = self.gibbs_free_energy(pressure, temperature, volume, params)
        entropy = self.entropy(pressure, temperature, volume, params)
        return gibbs + temperature * entropy

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The differentiation is tedious, so for now we just take the
        # numerical derivative of S
        # TODO Derive and use the analytical derivative.
        dT = 0.1
        S1 = self.entropy(pressure, temperature+dT/2., volume, params)
        S0 = self.entropy(pressure, temperature-dT/2., volume, params)
        return temperature * (S1 - S0)/dT

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return ((params['Cp'][0] * temperature
                 + 0.5 * params['Cp'][1] * np.power(temperature, 2.)
                 - params['Cp'][2] / temperature
                 + 2. * params['Cp'][3] * np.sqrt(temperature))
                - (params['Cp'][0] * params['T_0']
                   + 0.5 * params['Cp'][1] * params['T_0'] * params['T_0']
                   - params['Cp'][2] / params['T_0']
                   + 2.0 * params['Cp'][3] * np.sqrt(params['T_0'])))

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return ((params['Cp'][0] * np.log(temperature)
                 + params['Cp'][1] * temperature
                 - 0.5 * params['Cp'][2] / np.power(temperature, 2.)
                 - 2.0 * params['Cp'][3] / np.sqrt(temperature))
                - (params['Cp'][0] * np.log(params['T_0'])
                   + params['Cp'][1] * params['T_0']
                   - 0.5 * params['Cp'][2] / (params['T_0'] * params['T_0'])
                   - 2.0 * params['Cp'][3] / np.sqrt(params['T_0'])))

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if 'T_0' not in params:
            params['T_0'] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if 'H_0' not in params:
            params['H_0'] = float('nan')
        if 'S_0' not in params:
            params['S_0'] = float('nan')

        # First, let's check the EoS parameters for Tref
        mt.MT.validate_parameters(mt.MT(), params)
        self.static_params = {'V_0': params['V_0'],
                              'K_0': params['K_0'],
                              'Kprime_0': params['Kprime_0'],
                              'Kdprime_0': params['Kdprime_0'],
                              'P_0': params['P_0']}

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ['H_0', 'S_0', 'V_0', 'Cp', 'a_0', 'dKdT_0',
                         'n', 'molar_mass']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn('Unusual value for T_0', stacklevel=2)
        if ((params['G_0'] is not float('nan')
             and (params['G_0'] < 0. or params['G_0'] > 1.e13))):
            warnings.warn('Unusual value for G_0', stacklevel=2)
        if ((params['Gprime_0'] is not float('nan')
             and (params['Gprime_0'] < -5. or params['Gprime_0'] > 10.))):
            warnings.warn('Unusual value for Gprime_0', stacklevel=2)

        # no test for H_0 or S_0 (several HP endmembers have S_0 < 0)
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-2:
            warnings.warn('Unusual value for V_0', stacklevel=2)

        if self.molar_heat_capacity_p0(params['T_0'], params) < 0.:
            warnings.warn('Negative heat capacity at T_0', stacklevel=2)
        if self.molar_heat_capacity_p0(2000., params) < 0.:
            warnings.warn('Negative heat capacity at 2000K', stacklevel=2)

        if params['a_0'] < 0. or params['a_0'] > 1.e-3:
            warnings.warn('Unusual value for a_0', stacklevel=2)

        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn('Unusual value for n', stacklevel=2)
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 10.:
            warnings.warn('Unusual value for molar_mass', stacklevel=2)


class HP98(eos.EquationOfState):

    """
    Base class for the thermal equation of state
    described in :cite:`HP1998`.

    An instance "m" of a Mineral can be assigned this
    equation of state with the command m.set_method('hp98')
    (or by initialising the class with the param
    equation_of_state = 'hp98'
    """

    def _V_T_1bar(self, temperature, params):
        return params['V_0']*(1. + params['a_0']*(temperature - params['T_0'])
                              - 20.*params['a_0']*(np.sqrt(temperature)
                                                   - np.sqrt(params['T_0'])))

    def _K_T_1bar(self, temperature, params):
        return params['K_0'] + params['dKdT_0']*(temperature - params['T_0'])

    def volume(self, pressure, temperature, params):
        """
        Returns volume [m^3] as a function of pressure [Pa] and temperature [K]
        """
        return murn.volume(pressure,
                           self._V_T_1bar(temperature, params),
                           self._K_T_1bar(temperature, params),
                           params['Kprime_0'])

    def pressure(self, temperature, volume, params):
        """
        Returns pressure [Pa] as a function of temperature [K] and volume[m^3]
        """
        return murn.pressure(volume,
                             self._V_T_1bar(temperature, params),
                             self._K_T_1bar(temperature, params),
                             params['Kprime_0'])

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter [unitless] as a function of pressure,
        temperature, and volume.
        """
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        C_V = self.molar_heat_capacity_v(pressure, temperature,
                                         volume, params)
        return alpha * K_T * volume / C_V

    def isothermal_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns isothermal bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        return murn.bulk_modulus(pressure,
                                 self._K_T_1bar(temperature, params),
                                 params['Kprime_0'])

    # calculate the shear modulus as a function of P, V, and T
    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Not implemented.
        Returns 0.
        Could potentially apply a fixed Poissons ratio as a rough estimate.
        """
        return 0.

    # Cv, heat capacity at constant volume
    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume at the pressure, temperature,
        and volume [J/K/mol].
        """
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        V = self.volume(pressure, temperature, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        return C_p - V * temperature * alpha * alpha * K_T

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity at the pressure, temperature,
        and volume [1/K]
        """
        VT = self._V_T_1bar(temperature, params)
        KT = self._K_T_1bar(temperature, params)
        volume = murn.volume(pressure, VT, KT, params['Kprime_0'])
        dVTdT = params['V_0']*params['a_0']*(1. - 10./np.sqrt(temperature))
        g = volume / VT
        dgdKT = (pressure
                 * np.power(1. + pressure*params['Kprime_0']/KT,
                            -1-1./params['Kprime_0'])
                 / (KT*KT))
        dVdT = dVTdT*g + VT*dgdKT*params['dKdT_0']
        return dVdT / volume

    def molar_heat_capacity_p0(self, temperature, params):
        """
        Returns heat capacity at ambient pressure as a function of temperature
        [J/K/mol]
        Cp = a + bT + cT^-2 + dT^-0.5 in :cite:`HP1998`.
        """
        Cp = (params['Cp'][0]
              + params['Cp'][1] * temperature
              + params['Cp'][2] * np.power(temperature, -2.)
              + params['Cp'][3] * np.power(temperature, -0.5))
        return Cp

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus [Pa] as a function of pressure [Pa],
        temperature [K], and volume [m^3].
        """
        K_T = self.isothermal_bulk_modulus(pressure, temperature,
                                           volume, params)
        C_p = self.molar_heat_capacity_p(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        K_S = K_T * C_p / C_v
        return K_S

    def gibbs_free_energy(self, pressure, temperature, volume, params):
        """
        Returns the gibbs free energy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        return (params['H_0'] + self.__intCpdT(temperature, params)
                - temperature * (params['S_0']
                                 + self.__intCpoverTdT(temperature, params))
                + murn.intVdP(pressure,
                              self._V_T_1bar(temperature, params),
                              self._K_T_1bar(temperature, params),
                              params['Kprime_0']))

    def helmholtz_free_energy(self, pressure, temperature, volume, params):
        return (self.gibbs_free_energy(pressure, temperature, volume, params)
                - pressure * self.volume(pressure, temperature, params))

    def entropy(self, pressure, temperature, volume, params):
        """
        Returns the entropy [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The entropy involves differentiating intdVdp
        # with respect to temperature.
        # We do this using the product and chain rules

        VT = self._V_T_1bar(temperature, params)
        KT = self._K_T_1bar(temperature, params)
        dVTdT = params['V_0']*params['a_0']*(1. - 10./np.sqrt(temperature))
        g = murn.intVdP(pressure, VT, KT, params['Kprime_0']) / VT
        dgdKT = (((pressure/KT + 1.)
                 * np.power(1. + pressure*params['Kprime_0']/KT,
                            -1./params['Kprime_0']) - 1.)
                 / (params['Kprime_0'] - 1.))
        dintVdpdT = dVTdT*g + VT*dgdKT*params['dKdT_0']
        return (params['S_0'] + self.__intCpoverTdT(temperature, params)
                - dintVdpdT)

    def enthalpy(self, pressure, temperature, volume, params):
        """
        Returns the enthalpy [J/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        gibbs = self.gibbs_free_energy(pressure, temperature, volume, params)
        entropy = self.entropy(pressure, temperature, volume, params)
        return gibbs + temperature * entropy

    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns the heat capacity [J/K/mol] as a function of pressure [Pa]
        and temperature [K].
        """
        # The differentiation is tedious, so for now we just take the
        # numerical derivative of S
        # TODO calculate the analytical derivative
        dT = 0.1
        S1 = self.entropy(pressure, temperature+dT/2., volume, params)
        S0 = self.entropy(pressure, temperature-dT/2., volume, params)
        return temperature * (S1 - S0)/dT

    def __intCpdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state enthalpy [J/mol]
        at ambient pressure [Pa]
        """
        return ((params['Cp'][0] * temperature
                 + 0.5 * params['Cp'][1] * np.power(temperature, 2.)
                 - params['Cp'][2] / temperature
                 + 2. * params['Cp'][3] * np.sqrt(temperature))
                - (params['Cp'][0] * params['T_0']
                   + 0.5 * params['Cp'][1] * params['T_0'] * params['T_0']
                   - params['Cp'][2] / params['T_0']
                   + 2.0 * params['Cp'][3] * np.sqrt(params['T_0'])))

    def __intCpoverTdT(self, temperature, params):
        """
        Returns the thermal addition to the standard state entropy [J/K/mol]
        at ambient pressure [Pa]
        """
        return ((params['Cp'][0] * np.log(temperature)
                 + params['Cp'][1] * temperature
                 - 0.5 * params['Cp'][2] / np.power(temperature, 2.)
                 - 2.0 * params['Cp'][3] / np.sqrt(temperature))
                - (params['Cp'][0] * np.log(params['T_0'])
                   + params['Cp'][1] * params['T_0']
                   - 0.5 * params['Cp'][2] / (params['T_0'] * params['T_0'])
                   - 2.0 * params['Cp'][3] / np.sqrt(params['T_0'])))

    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """
        if 'T_0' not in params:
            params['T_0'] = 298.15

        # If standard state enthalpy and entropy are not included
        # this is presumably deliberate, as we can model density
        # and bulk modulus just fine without them.
        # Just add them to the dictionary as nans.
        if 'H_0' not in params:
            params['H_0'] = float('nan')
        if 'S_0' not in params:
            params['S_0'] = float('nan')

        # First, let's check the EoS parameters for Tref
        murn.Murnaghan.validate_parameters(murn.Murnaghan(), params)

        # Now check all the required keys for the
        # thermal part of the EoS are in the dictionary
        expected_keys = ['H_0', 'S_0', 'V_0', 'Cp', 'a_0', 'dKdT_0',
                         'n', 'molar_mass']
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)

        # Finally, check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn('Unusual value for T_0', stacklevel=2)
        if ((params['G_0'] is not float('nan')
             and (params['G_0'] < 0. or params['G_0'] > 1.e13))):
            warnings.warn('Unusual value for G_0', stacklevel=2)
        if ((params['Gprime_0'] is not float('nan')
             and (params['Gprime_0'] < -5. or params['Gprime_0'] > 10.))):
            warnings.warn('Unusual value for Gprime_0', stacklevel=2)

        # no test for H_0 or S_0 (several HP endmembers have S_0 < 0)
        if params['V_0'] < 1.e-7 or params['V_0'] > 1.e-2:
            warnings.warn('Unusual value for V_0', stacklevel=2)

        if self.molar_heat_capacity_p0(params['T_0'], params) < 0.:
            warnings.warn('Negative heat capacity at T_0', stacklevel=2)
        if self.molar_heat_capacity_p0(2000., params) < 0.:
            warnings.warn('Negative heat capacity at 2000K', stacklevel=2)

        if params['a_0'] < 0. or params['a_0'] > 1.e-3:
            warnings.warn('Unusual value for a_0', stacklevel=2)

        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn('Unusual value for n', stacklevel=2)
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 10.:
            warnings.warn('Unusual value for molar_mass', stacklevel=2)
