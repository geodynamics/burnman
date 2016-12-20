# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.

from __future__ import absolute_import

import numpy as np
from scipy.optimize import fsolve
import warnings

from . import equation_of_state as eos
from ..constants import gas_constant

class AA(eos.EquationOfState):
    """
    Base class for the liquid metal EOS detailed in Anderson and Ahrens (1994).

    This equation of state is described by a fourth order BM isentropic EoS
    (V_0, KS, KS', KS''), a description of the volumetric heat capacity
    (which gives energy along an isochore as a function of temperature), and
    a description of the gruneisen parameter as a function of volume and energy
    (which defines the thermal pressure).
    """

    # Electronic heat capacity functions
    def _ABTheta(self, V, params):
        Vfrac = V/params['V_0']
        
        A = params['a'][0] + params['a'][1]*Vfrac
        B = params['b'][0] + params['b'][1]*Vfrac*Vfrac
        Theta = params['Theta'][0]*np.power(Vfrac, -params['Theta'][1])

        return A, B, Theta

    # Potential heat capacity functions
    def _lambdaxi(self, V, params):
        rhofrac = params['V_0']/V
        xi = params['xi_0']*np.power(rhofrac, -0.6)
        F = 1./(1. + np.exp((rhofrac - params['F'][0])/params['F'][1]))
        lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac) + params['lmda'][2])*np.power(rhofrac, 0.4)
        #lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac + params['lmda'][2]))*np.power(rhofrac, 0.4) # this incorrect expression for lmda seems to provide a very close fit to figure 5

        return lmda, xi
    
    # Fourth order BM functions
    def _rhofracxksis(self, V, params):
        rhofrac = params['V_0']/V # rho/rho0 = V0/V
        x = np.power(rhofrac, 1./3.) # equation 18
        ksi1 = 0.75*(4. - params['Kprime_S']) # equation 19
        ksi2 = 0.375*(params['K_S']*params['Kprime_prime_S'] + params['Kprime_S']*(params['Kprime_S'] - 7.)) + 143./24. # equation 20
        return rhofrac, x, ksi1, ksi2

    '''     
    Contributions to the heat capacity
    '''
    
    # High temperature limit of the kinetic contribution to the heat capacity
    # Anderson and Ahrens (1994), just after equation 29.
    def _C_v_kin(self, V, T, params):
        return 1.5*params['n']*gas_constant
    
    # Equation A1
    def _C_v_el(self, V, T, params):
        A, B, Theta = self._ABTheta(V, params)
        C_e = A*(1. - (Theta*Theta)/(Theta*Theta + T*T)) + B*np.power(T, 0.6)
        return C_e

    # Equation A15
    def _C_v_pot(self, V, T, params):
        lmda, xi = self._lambdaxi(V, params)
        C_pot = (lmda*T + xi*params['theta']) / (params['theta'] + T)
        return C_pot

    '''
    Contributions to the internal energy
    '''
        
    def _internal_energy_kin(self, Ts, T, V, params):
        E_kin = 1.5*params['n']*gas_constant*(T - Ts)
        return E_kin
    
    def _internal_energy_el(self, Ts, T, V, params):
        A, B, Theta = self._ABTheta(V, params)
        E_el = A*(T - Ts - Theta*(np.arctan(T/Theta) - np.arctan(Ts/Theta))) + 0.625*B*(np.power(T, 1.6) - np.power(Ts, 1.6))
        return E_el
    
    def _internal_energy_pot(self, Ts, T, V, params):
        lmda, xi = self._lambdaxi(V, params)
        E_pot = (lmda*(T - Ts) + params['theta']*(xi - lmda)*np.log((params['theta'] + T)/(params['theta'] + Ts)))
        return E_pot
    
    '''
    Contributions to entropy
    '''
    
    def _entropy_kin(self, Ts, T, V, params):
        if np.abs(T- Ts) > 1.e-10:
            S_kin = 1.5*params['n']*gas_constant*(np.log(T) - np.log(Ts))
        else:
            S_kin = 0.
        return S_kin
        
    def _entropy_el(self, Ts, T, V, params):
        if np.abs(T- Ts) > 1.e-10:
            A, B, Theta = self._ABTheta(V, params)
            S_el = (A*(np.log(T/Ts) - 0.5*np.log(T*T*(Theta*Theta + Ts*Ts)/(Ts*Ts*(Theta*Theta + T*T)))) + 5./3.*B*(np.power(T, 0.6) - np.power(Ts, 0.6)))
        else:
            S_el = 0.
        return S_el
    
    def _entropy_pot(self, Ts, T, V, params):
        if np.abs(T- Ts) > 1.e-10:
            lmda, xi = self._lambdaxi(V, params)
            S_pot = (lmda*np.log((params['theta'] + T)/(params['theta'] + Ts)) + xi*np.log((T*(params['theta'] + Ts))/(Ts*(params['theta'] + T))))
        else:
            S_pot = 0.
        return S_pot
            
    '''
    Isentropic and isochoric calculations
    '''
    
    # Temperature along an isentrope (Anderson and Ahrens; Equation B5)
    def _isentropic_temperature(self, V, params):
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



    # Pressure along the reference isentrope
    def _isentropic_pressure(self, V, params):
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(V, params)
        x2 = x*x
        x3 = x*x*x
        x5 = x3*x2
        x7 = x5*x2
    
        Ps = 1.5*params['K_S'] * (x7 - x5) * (1. + ksi1 - ksi1*x2 + ksi2 * (x2 - 1.) * (x2 - 1.))
    
        return Ps

    # Birch Murnaghan equation of state expression for the energy change along an isentrope
    # Anderson and Ahrens, 1994 (Equation 21)
    def _isentropic_energy_change(self, V, params):
        rhofrac, x, ksi1, ksi2 = self._rhofracxksis(V, params)
        x2 = x*x
        x4 = x2*x2
        x6 = x4*x2
        x8 = x4*x4
        
        E_S = 4.5*params['V_0']*params['K_S'] * ((ksi1 + 1.) * (x4/4. - x2/2. + 0.25) - ksi1*(x6/6. - x4/4. + 1./12.)
                                                + ksi2*(x8/8. - x6/2. + 0.75*x4 - x2/2. + 0.125))
        return E_S

    # int Cv dT 
    def _isochoric_energy_change(self, Ts, T, V, params):
        return (self._internal_energy_kin(Ts, T, V, params)
                + self._internal_energy_el(Ts, T, V, params)
                + self._internal_energy_pot(Ts, T, V, params))
    
    
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


    def _volume(self, volume, pressure, temperature, params):
        return pressure - self.pressure(temperature, volume, params)

    def volume(self, pressure, temperature, params):
        """
        Returns molar volume. :math:`[m^3]`
        """
        return fsolve(self._volume, params['V_0']*0.1, args=(pressure, temperature, params))[0]

    def pressure( self, temperature, volume, params):
        """
        Returns the pressure of the mineral at a given temperature and volume [Pa]
        """
        
        Ts = self._isentropic_temperature(volume, params)
        dE = self._isochoric_energy_change(Ts, temperature, volume, params)
        E1 = self._isentropic_energy_change(volume, params) # should also include params['E_0'] given the expression in Anderson and Ahrens. Here, we take the energy change relative to the reference isentrope (effective E_0 = 0). The energy at standard state is *only* used to calculate the final energies, not the physical properties.
        E2 = E1 + dE

        # Integrate along isochore (\int dP = \int gr/V dE)
        dP = (params['grueneisen_0']*dE
              + 0.5*params['grueneisen_prime']*np.power(params['V_0']/volume,
                                                        params['grueneisen_n'])*(E2*E2 - E1*E1))/volume
        P = self._isentropic_pressure(volume, params) + dP
    
        return P

        

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]` 
        """
        gr = (params['grueneisen_0'] +
              params['grueneisen_prime'] *
              (np.power(params['V_0']/volume,
                        params['grueneisen_n']) *
               self.internal_energy(pressure, temperature, volume, params)))
        return gr

    def isothermal_bulk_modulus(self, pressure,temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]` 
        """
        # K_T = -V * dP/dV
        delta_V = params['V_0']*1.e-5
        delta_P = self.pressure(temperature, volume+delta_V, params) - pressure

        K_T = -volume*delta_P/delta_V
        
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

    def heat_capacity_v(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant volume. :math:`[J/K/mol]` 
        """
        C_v = (self._C_v_kin(volume, temperature, params)
                + self._C_v_pot(volume, temperature, params)
                + self._C_v_el(volume, temperature, params))
        return C_v
    
    def heat_capacity_p(self, pressure, temperature, volume, params):
        """
        Returns heat capacity at constant pressure. :math:`[J/K/mol]` 
        """
        
        alpha = self.thermal_expansivity(pressure, temperature, volume, params)
        gr = self.grueneisen_parameter(pressure, temperature, volume, params)
        C_v = self.heat_capacity_v(pressure, temperature, volume, params)
        C_p = C_v*(1. + gr * alpha * temperature)
        
        return C_p

    def thermal_expansivity(self, pressure, temperature, volume, params):
        """
        Returns thermal expansivity. :math:`[1/K]`
        Currently found by numerical differentiation (1/V * dV/dT) 
        """

        delta_T = 1.
        delta_V = ( self.volume(pressure, temperature+delta_T, params) - volume )
            
        return (1./volume)*delta_V/delta_T

    def gibbs_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Gibbs free energy at the pressure and temperature of the mineral [J/mol]
        E + PV
        """
        
        return self.helmholtz_free_energy( pressure, temperature, volume, params) + \
            pressure * self.volume( pressure, temperature, params)

    def internal_energy( self, pressure, temperature, volume, params):
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
        Ts = self._isentropic_temperature(volume, params)
        S = (params['S_0'] + self._entropy_kin(Ts, temperature, volume, params)
             + self._entropy_el(Ts, temperature, volume, params)
             + self._entropy_pot(Ts, temperature, volume, params))

        return S 

    def enthalpy( self, pressure, temperature, volume, params):
        """
        Returns the enthalpy at the pressure and temperature of the mineral [J/mol]
        E + PV
        """
        
        return self.internal_energy(pressure, temperature, volume, params) + \
            pressure * self.volume( pressure, temperature, params)

    def helmholtz_free_energy( self, pressure, temperature, volume, params):
        """
        Returns the Helmholtz free energy at the pressure and temperature of the mineral [J/mol]
        E - TS
        """
        
        return self.internal_energy(pressure, temperature, volume, params) - temperature*self.entropy(pressure, temperature, volume, params)


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
