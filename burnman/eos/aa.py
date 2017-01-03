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

    This is an EVS equation of state. Internal energy (E) is first calculated 
    along a reference isentrope using a fourth order BM EoS
    (V_0, KS, KS', KS''), which gives volume as a function of pressure, 
    and the thermodynamic identity:

    -dE/dV | S = P.

    The temperature along the isentrope is calculated via

    d(ln T)/d(ln rho) | S = grueneisen

    which gives:
    Ts/T0 = exp(int(grueneisen/rho drho))
    
    The thermal effect on internal energy is calculated at constant volume 
    using expressions for the kinetic, electronic and potential contributions 
    to the volumetric heat capacity, which can then be integrated with respect 
    to temperature:

    dE/dT | V = Cv
    dE/dS | V = T

    We note that Anderson and Ahrens (1994) also include a detailed description
    of the gruneisen parameter as a function of volume and energy, and incorporate
    their preferred expression into the table which provides all of their 
    parameters for the equation of state. However, this expression is not 
    required to formulate the equation of state 
    (other than finding the temperature along the principal isentrope), 
    and indeed is generally inconsistent with it. This can be seen 
    most simply by considering the following:

    1) As energy and entropy are defined by the equation of state at any 
    temperature and volume, the grueneisen parameter is also 
    (implicitly) defined via the expressions:

    dE = TdS - PdV (and so dE/dV | S = P)
    grueneisen = V dP/dE | V 

    Away from the reference isentrope, the grueneisen parameter 
    calculated using these expressions is not equal to the
    analytical expression given by Anderson and Ahrens (1994).
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
        lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac) + params['lmda'][2])*np.power(rhofrac, 0.4) # A17
        #lmda = (F*(params['lmda'][0] + params['lmda'][1]*rhofrac + params['lmda'][2]))*np.power(rhofrac, 0.4) # this incorrect expression for lmda seems to provide a very close fit to figure 5

        return lmda, xi
    
    def _rhofracxksis(self, V, params):
        """
        Functions for the fourth order Birch-Murnaghan equation of state
        """
        rhofrac = params['V_0']/V # rho/rho0 = V0/V
        x = np.power(rhofrac, 1./3.) # equation 18
        ksi1 = 0.75*(4. - params['Kprime_S']) # equation 19
        ksi2 = 0.375*(params['K_S']*params['Kprime_prime_S'] + params['Kprime_S']*(params['Kprime_S'] - 7.)) + 143./24. # equation 20
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
    
        Ps = 1.5*params['K_S'] * (x7 - x5) * (1. + ksi1 - ksi1*x2 + ksi2 * (x2 - 1.) * (x2 - 1.)) # Eq. 17
    
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
        
        E_S = 4.5*params['V_0']*params['K_S'] * ((ksi1 + 1.) * (x4/4. - x2/2. + 0.25) - ksi1*(x6/6. - x4/4. + 1./12.)
                                                + ksi2*(x8/8. - x6/2. + 0.75*x4 - x2/2. + 0.125)) # Eq. 21
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

        '''
        dE = self._isochoric_energy_change(Ts, temperature, volume, params)
        E1 = self._isentropic_energy_change(volume, params) # should also include params['E_0'] given the expression in Anderson and Ahrens. Here, we take the energy change relative to the reference isentrope (effective E_0 = 0). The energy at standard state is *only* used to calculate the final energies, not the physical properties.
        E2 = E1 + dE

        # Integrate at constant volume (V \int dP = \int gr dE)
        dP = (params['grueneisen_0']*(E2 - E1) +
              (0.5*params['grueneisen_prime'] *
               np.power(params['V_0']/volume, params['grueneisen_n']) *
               (E2*E2 - E1*E1))) / volume # eq. 23
              
        P = self._isentropic_pressure(volume, params) + dP
        '''

        dV = volume*1.e-5
        S0 = self.entropy(0., temperature, volume, params)
        E0 = self.internal_energy(0., temperature, volume, params)

        def Sdiff(args, S, dV):
            T = args[0]
            S1 = self.entropy(0., T, volume+dV, params)
            return S1 - S0

        T1 = fsolve(Sdiff, [temperature], args=(S0, dV))[0]
                  
        S1 = self.entropy(0., T1, volume+dV, params)
        E1 = self.internal_energy(0., T1, volume+dV, params)

        P = -(E1 - E0)/dV
        return P

        

    def grueneisen_parameter(self, pressure, temperature, volume, params):
        """
        Returns grueneisen parameter :math:`[unitless]` 
        """

        '''
        gr = (params['grueneisen_0'] +
              params['grueneisen_prime'] *
              (np.power(params['V_0']/volume, params['grueneisen_n']) *
               self.internal_energy(pressure, temperature, volume, params)))
        '''
        dT = 1.
        dE = (self.internal_energy(0., temperature, volume, params) -
                                   self.internal_energy(0., temperature+dT, volume, params))
        dP = (self.pressure(temperature, volume, params) -
              self.pressure(temperature+dT, volume, params))
        gr = volume*dP/dE
        return gr

    def isothermal_bulk_modulus(self, pressure,temperature, volume, params):
        """
        Returns isothermal bulk modulus :math:`[Pa]` 
        """
        # K_T = -V * dP/dV
        delta_V = params['V_0']*1.e-4
        delta_P = self.pressure(temperature, volume+delta_V, params) - pressure

        K_T = -(volume + 0.5*delta_V)*delta_P/delta_V
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

        A, B, Theta = self._ABTheta(volume, params)
        lmda, xi = self._lambdaxi(volume, params)
        
        C_kin = 1.5*params['n']*gas_constant # HT limit of kinetic contribution (just after equation 29.)
        C_e = A*(1. - (Theta*Theta)/(Theta*Theta + temperature*temperature)) + B*np.power(temperature, 0.6) # Equation A1
        C_pot = (lmda*temperature + xi*params['theta']) / (params['theta'] + temperature) # Equation A15
        
        return C_kin + C_e + C_pot

    
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
