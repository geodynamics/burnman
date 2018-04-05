# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from scipy.optimize import brentq
import warnings

from . import equation_of_state as eos
from ..constants import gas_constant

class AA(eos.EquationOfState):
    """
    Class for the :math`E-V-S` liquid metal EOS detailed in :cite:`AA1994`.
    Internal energy (:math:`E`) is first calculated 
    along a reference isentrope using a fourth order BM EoS
    (:math:`V_0`, :math:`KS`, :math:`KS'`, :math:`KS''`), 
    which gives volume as a function of pressure, 
    coupled with the thermodynamic identity:

    :math:`-\partial E/ \partial V |_S = P`.

    The temperature along the isentrope is calculated via

    :math:`\partial (\ln T)/\partial (\ln \\rho) |_S = \gamma`

    which gives:

    :math:`T_S/T_0 = \exp(\int( \gamma/\\rho ) d \\rho)`
    
    The thermal effect on internal energy is calculated at constant volume 
    using expressions for the kinetic, electronic and potential contributions 
    to the volumetric heat capacity, which can then be integrated with respect 
    to temperature:

    :math:`\partial E/\partial T |_V = C_V`

    :math:`\partial E/\partial S |_V = T`

    We note that :cite:`AA1994` also include a detailed description
    of the Gruneisen parameter as a function of volume and energy (Equation 15), 
    and use this to determine the temperature along the principal isentrope 
    (Equations B1-B10) and the thermal pressure away from that isentrope 
    (Equation 23). However, this expression is inconsistent with 
    the equation of state away from the principal isentrope. Here we choose 
    to calculate the thermal pressure and Grueneisen parameter thus:

    1) As energy and entropy are defined by the equation of state at any 
    temperature and volume, pressure can be found by via the expression:

    :math:`\partial E/\partial V |_S = P`

    2) The Grueneisen parameter can now be determined as
    :math:`\gamma = V \partial P/\partial E |_V`

    To reiterate: away from the reference isentrope, the Grueneisen parameter 
    calculated using these expressions is *not* equal to the 
    (thermodynamically inconsistent) analytical expression given by :cite:`AA1994`.

    A final note: the expression for :math:`\Lambda` (Equation 17).
    does not reproduce Figure 5. We assume here that the figure matches the model
    actually used by :cite:`AA1994`, which has the form:
    :math:`F(-325.23 + 302.07 (\\rho/\\rho_0) + 30.45 (\\rho/\\rho_0)^{0.4})`.
    """

    def _ABTheta(self, V, params):
        """
        Electronic heat capacity functions
        """
        Vfrac = V/params['V_0']
        
        A = params['a'][0] + params['a'][1]*Vfrac # A2
        B = params['b'][0] + params['b'][1]*Vfrac*Vfrac # A3
        Theta = params['Theta'][0]*np.power(Vfrac, -params['Theta'][1]) # A4

        return A, B, Theta

    def _lambdaxi(self, V, params):
        """
        Potential heat capacity functions
        """
        rhofrac = params['V_0']/V
        xi = params['xi_0']*np.power(rhofrac, -0.6) # A16
        F = 1./(1. + np.exp((rhofrac - params['F'][0])/params['F'][1])) # A18
        #lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac) + params['lmda'][2])*np.power(rhofrac, 0.4) # A17
        lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac + params['lmda'][2]*np.power(rhofrac, 0.4))) # this incorrect expression for lmda seems to provide a very close fit to figure 5

        return lmda, xi
    
    def _rhofracxksis(self, V, params):
        """
        Functions for the fourth order Birch-Murnaghan equation of state
        """
        rhofrac = params['V_0']/V # rho/rho0 = V0/V
        x = np.power(rhofrac, 1./3.) # equation 18
        ksi1 = 0.75*(4. - params['Kprime_S']) # equation 19
        ksi2 = 0.375*(params['K_S']*params['Kprime_prime_S'] +
                      params['Kprime_S']*(params['Kprime_S'] - 7.)) + 143./24. # equation 20
        return rhofrac, x, ksi1, ksi2
    

    def _isentropic_temperature(self, V, params):
        """
        Temperature along the reference isentrope
        """
        
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(V, params)
        
        # equation B6 -- B10
        a1 = ksi2 / 8.
        a2 = ( ksi1 + 3. * ksi2 ) / 6.
        a3 = ( 1. + 2.*ksi1 + 3.*ksi2 ) / 4.
        a4 = (1. + ksi1 + ksi2)/2.
        a5 = (6. + 4.*ksi1 + 3.*ksi2)/24.
    
        # equation B5
        Ts = params['T_0']*np.exp(params['grueneisen_0']*np.log(rhofrac)
                                + 13.5*params['grueneisen_prime']*params['V_0']*params['K_S'] *
                                (   (a1/(3*params['grueneisen_n'] + 8.))*(np.power(x,(3*params['grueneisen_n'] + 8.)) - 1.)
                                    - (a2/(3*params['grueneisen_n'] + 6.))*(np.power(x,(3*params['grueneisen_n'] + 6.)) - 1.)
                                + (a3/(3*params['grueneisen_n'] + 4.))*(np.power(x,(3*params['grueneisen_n'] + 4.)) - 1.)
                                - (a4/(3*params['grueneisen_n'] + 2.))*(np.power(x,(3*params['grueneisen_n'] + 2.)) - 1.)
                                + (a5/(3*params['grueneisen_n'] + 0.))*(np.power(x,(3*params['grueneisen_n'] + 0.)) - 1.)))
                                
        return Ts



    def _isentropic_pressure(self, V, params):
        """
        Pressure along the reference isentrope
        """
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(V, params)
        x2 = x*x
        x3 = x*x*x
        x5 = x3*x2
        x7 = x5*x2
    
        Ps = ( 1.5*params['K_S'] * (x7 - x5) *
               (1. + ksi1 - ksi1*x2 +
                ksi2 * (x2 - 1.) * (x2 - 1.)) ) # Eq. 17
    
        return Ps

    def _isentropic_energy_change(self, V, params):
        """
        Birch Murnaghan equation of state expression for the energy change along an isentrope
        """
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(V, params)
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        x8 = x4*x4
        
        E_S = 4.5*params['V_0']*params['K_S'] * ((ksi1 + 1.) * (x4/4. - x2/2. + 1./4.) -
                                                 ksi1*(x6/6. - x4/4. + 1./12.) +
                                                 ksi2*(x8/8. - x6/2. + 3.*x4/4. - x2/2. + 1./8.)) # Eq. 21
        return E_S

    def _isochoric_energy_change(self, Ts, T, V, params):
        """
        int Cv dT 
        """
        A, B, Theta = self._ABTheta(V, params)
        lmda, xi = self._lambdaxi(V, params)
        
        E_kin = 1.5*params['n']*gas_constant*(T - Ts)
        E_el = A*(T - Ts - Theta*(np.arctan(T/Theta) - np.arctan(Ts/Theta))) + 5./8*B*(np.power(T, 1.6) - np.power(Ts, 1.6)) # A5
        E_pot = (lmda*(T - Ts) + params['theta']*(xi - lmda)*np.log((params['theta'] + T)/(params['theta'] + Ts))) # A19

        return E_kin + E_el + E_pot

    
    def volume_dependent_q(self, x, params):
        """
        Finite strain approximation for :math:`q`, the isotropic volume strain
        derivative of the grueneisen parameter.
        """
        raise NotImplementedError("")

    def _isotropic_eta_s(self, x, params):
        """
        Finite strain approximation for :math:`eta_{s0}`, the isotropic shear
        strain derivative of the grueneisen parameter.
        Zero for a liquid
        """
        return 0.

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """

        _volume = lambda V, P, T, params: ( P -
                                            self.pressure(T, V, params) )
        
        return brentq(_volume, params['V_0']*0.1, params['V_0']*2., args=(pressure, temperature, params))

    def pressure( self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume [Pa]
        """
        
        '''
        Ts = self._isentropic_temperature(volume, params)
        
        
        dE = self._isochoric_energy_change(Ts, temperature, volume, params)
        E1 = self._isentropic_energy_change(volume, params) - params['E_0']
        E2 = E1 + dE

        # Integrate at constant volume (V \int dP = \int gr dE)
        dP = (params['grueneisen_0']*(E2 - E1) +
              (0.5*params['grueneisen_prime'] *
               np.power(params['V_0']/volume, params['grueneisen_n']) *
               (E2*E2 - E1*E1))) / volume # eq. 23

        P = self._isentropic_pressure(volume, params) + dP
        '''

        dV = volume*1.e-4
        S = self.entropy(0., temperature, volume, params)

        delta_S = lambda T, S, V: S - self.entropy(0., T, V, params)
        
        T0 = brentq(delta_S, temperature*0.97, temperature*1.03, args=(S, volume - 0.5*dV))
        T1 = brentq(delta_S, temperature*0.97, temperature*1.03, args=(S, volume + 0.5*dV))

        E0 = self.molar_internal_energy(0., T0, volume - 0.5*dV, params)
        E1 = self.molar_internal_energy(0., T1, volume + 0.5*dV, params)

        P = -(E1 - E0)/dV # |S
        
        return P

        

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]` 
        """
        '''
        gr = (params['grueneisen_0'] +
              params['grueneisen_prime'] *
              (np.power(params['V_0']/volume, params['grueneisen_n']) *
               (self.molar_internal_energy(pressure, temperature, volume, params) -
                params['E_0'])))
        '''
        dT = 1.
        dE = (self.molar_internal_energy(0., temperature + 0.5*dT, volume, params) -
              self.molar_internal_energy(0., temperature - 0.5*dT, volume, params))
        dP = (self.pressure(temperature + 0.5*dT, volume, params) -
              self.pressure(temperature - 0.5*dT, volume, params))
        gr = volume*dP/dE
        
        return gr

    def isothermal_bulk_modulus(self, pressure,temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]` 
        """
        # K_T = -V * dP/dV
        dV = volume*1.e-3
        P0 = self.pressure(temperature, volume - 0.5*dV, params)
        P1 = self.pressure(temperature, volume + 0.5*dV, params)

        K_T = -volume*(P1 - P0)/dV
        return K_T

    def adiabatic_bulk_modulus(self, pressure, temperature, volume, params):
        """
        Returns adiabatic bulk modulus. :math:`[Pa]` 
        """
        K_T=self.isothermal_bulk_modulus(pressure, temperature, volume, params)
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        K_S = K_T*(1. + gr * alpha * temperature)
        return K_S

    def shear_modulus(self, pressure, temperature, volume, params):
        """
        Returns shear modulus. :math:`[Pa]` 
        Zero for a liquid
        """
        return 0.

    def molar_heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]` 
        """

        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)
        
        C_kin = 1.5*params['n']*gas_constant # HT limit of kinetic contribution (just after equation 29.)
        C_e = A*(1. - (Theta*Theta)/(Theta*Theta + temperature*temperature)) + B*np.power(temperature, 0.6) # Equation A1
        C_pot = (lmda*temperature + xi*params['theta']) / (params['theta'] + temperature) # Equation A15
        
        return C_kin + C_e + C_pot

    
    def molar_heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]` 
        """
        
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.molar_heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v*(1. + gr * alpha * temperature)
        
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        Currently found by numerical differentiation (1/V * dV/dT) 
        """

        delta_T = 1.
        V0 = self.volume(pressure, temperature-0.5*delta_T, params)
        V1 = self.volume(pressure, temperature+0.5*delta_T, params)
            
        return (1./volume)*(V1 - V0)/delta_T

    def gibbs_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        E + PV
        """
        
        return self.helmholtz_free_energy( pressure, temperature, volume, params) + \
            pressure * self.volume( pressure, temperature, params)

    def molar_internal_energy( self, pressure, temperature, volume, params):
        """
        Returns the internal energy at the pressure and temperature of the mineral [J/mol]
        """
        Ts = self._isentropic_temperature(volume, params)
        E = (params['E_0'] + self._isentropic_energy_change(volume, params)
                + self._isochoric_energy_change(Ts, temperature, volume, params))
            
        return E

    def entropy( self, pressure, temperature, volume, params):
        """
        Returns the entropy at the pressure and temperature of the mineral [J/K/mol]
        """
        T = temperature
        Ts = self._isentropic_temperature(volume, params)

        if np.abs(T- Ts) < 1.e-10:
            Delta_S = 0.
        else:
            A, B, Theta = self._ABTheta(volume, params)
            lmda, xi = self._lambdaxi(volume, params)
            S_kin = 1.5*params['n']*gas_constant*(np.log(T) - np.log(Ts))
            S_el = (A*(np.log(T/Ts) - 0.5*np.log(T*T*(Theta*Theta + Ts*Ts)/(Ts*Ts*(Theta*Theta + T*T)))) + 5./3.*B*(np.power(T, 0.6) - np.power(Ts, 0.6))) # A6
            S_pot = (lmda*np.log((params['theta'] + T)/(params['theta'] + Ts)) + xi*np.log((T*(params['theta'] + Ts))/(Ts*(params['theta'] + T)))) # A20
            Delta_S = S_kin + S_el + S_pot
        
        S = params['S_0'] + Delta_S
        return S 

    def enthalpy( self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        E + PV
        """
        
        return self.molar_internal_energy(pressure, temperature, volume, params) + \
            pressure * self.volume( pressure, temperature, params)

    def helmholtz_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        E - TS
        """
        return self.molar_internal_energy(pressure, temperature, volume, params) - temperature*self.entropy(pressure, temperature, volume, params)


    def validate_parameters(self, params):
        """
        Check for existence and validity of the parameters
        """

        # Now check all the required keys for the 
        # thermal part of the EoS are in the dictionary
        expected_keys = ['P_0', 'T_0', 'S_0', 'molar_mass', 'grueneisen_0']

        
        for k in expected_keys:
            if k not in params:
                raise KeyError('params object missing parameter : ' + k)
        
        # Finally, check that the values are reasonable.
        if params['T_0'] < 0.:
            warnings.warn( 'Unusual value for T_0', stacklevel=2 )
        if params['molar_mass'] < 0.001 or params['molar_mass'] > 10.:
            warnings.warn( 'Unusual value for molar_mass', stacklevel=2 )
        if params['n'] < 1. or params['n'] > 1000.:
            warnings.warn( 'Unusual value for n', stacklevel=2 )
