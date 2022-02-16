# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import scipy.optimize as opt
from ..utils.math import bracket


class Calibrant(object):
    """
    The base class for a pressure calibrant material.

    Initialization of a Calibrant object requires the following parameters:

    calibrant_function : Function
        A function that takes either pressure, temperature and a params
        object as arguments, returning the volume, or takes volume, temperature
        and a params object, returning the pressure.

    calibrant_function_return_type : 'pressure' or 'volume'
        The return type of the calibrant function.

    params : dictionary
        A dictionary containing the parameters required by the
        calibrant function.
    """
    def __init__(self, calibrant_function,
                 calibrant_function_return_type,
                 params):

        if calibrant_function_return_type == 'pressure':
            self.pressure_function = calibrant_function
            self.volume_function = self._volume_using_pressure_function
        elif calibrant_function_return_type == 'volume':
            self.volume_function = calibrant_function
            self.pressure_function = self._pressure_using_volume_function
        else:
            raise Exception('calibrant function return type must either '
                            'be pressure or volume')

        self.calibrant_function = calibrant_function

        self.params = params

    def _volume_using_pressure_function(self, pressure, temperature, params):
        """
        Helper function to compute volume iteratively by Brent's method using
        a function of the form pressure(volume, temperature).
        """
        def func(x):
            return (self.pressure_function(x, temperature, params) - pressure)
        try:
            sol = bracket(func, params['V_0'], 1.e-2 * params['V_0'])
        except ValueError:
            raise ValueError('Cannot find a volume, perhaps you are outside '
                             'of the range of validity for the equation '
                             'of state?')
        return opt.brentq(func, sol[0], sol[1])

    def _pressure_using_volume_function(self, volume, temperature, params):
        """
        Helper function to compute pressure iteratively by Brent's method using
        a function of the form volume(pressure, temperature).
        """
        def func(x):
            (self.volume_function(x, temperature, params) - volume)
        try:
            sol = bracket(func, 0., 300.e9)
        except ValueError:
            raise ValueError('Cannot find a pressure, perhaps you are outside '
                             'of the range of validity for the equation '
                             'of state?')
        return opt.brentq(func, sol[0], sol[1])

    def pressure(self, volume, temperature, VT_covariance=None):
        """
        Returns the pressure of the calibrant as a function of
        volume, temperature and (optionally) a volume-temperature
        variance-covariance matrix.

        Parameters
        ----------
        volume : float
            The volume of the calibrant [m^3/mol]
        temperature : float
            The temperature of the calibrant [K]
        VT_covariance : 2x2 numpy array [optional]
            The volume-temperature
            variance-covariance matrix

        Returns
        -------
        pressure : float
            The pressure of the calibrant [Pa]

        PVT_covariance : 3x3 numpy array (if PT_covariance is provided)
            The pressure-volume-temperature variance-covariance matrix.
        """
        if VT_covariance is None:
            return self.pressure_function(volume, temperature, self.params)
        else:
            # Here we take the centered differences
            # We could alternatively use thermodynamic properties
            # but these have not yet been implemented.
            dV = volume / 1.e7
            dT = 0.01
            PdV0 = self.pressure_function(volume - dV/2., temperature,
                                          self.params)
            PdV1 = self.pressure_function(volume + dV/2., temperature,
                                          self.params)
            PdT0 = self.pressure_function(volume, temperature - dT/2.,
                                          self.params)
            PdT1 = self.pressure_function(volume, temperature + dT/2.,
                                          self.params)
            pressure = (PdV0 + PdV1 + PdT0 + PdT1)/4.

            gradPVT = np.zeros((2, 3))
            gradPVT[:, 1:] = np.eye(2)
            gradPVT[:, 0] = [(PdV1 - PdV0)/dV, (PdT1 - PdT0)/dT]
            PVT_covariance = gradPVT.T.dot(VT_covariance).dot(gradPVT)
            return pressure, PVT_covariance

    def volume(self, pressure, temperature, PT_covariance=None):
        """
        Returns the volume of the calibrant as a function of
        pressure, temperature and (optionally) a pressure-temperature
        variance-covariance matrix.

        Parameters
        ----------
        pressure : float
            The pressure of the calibrant [Pa]
        temperature : float
            The temperature of the calibrant [K]
        PT_covariance : 2x2 numpy array [optional]
            The pressure-temperature
            variance-covariance matrix

        Returns
        -------
        volume : float
            The volume of the calibrant [m^3/mol]

        VPT_covariance : 3x3 numpy array (if VT_covariance is provided)
            The volume-pressure-temperature variance-covariance matrix.
        """
        if PT_covariance is None:
            return self.volume_function(pressure, temperature, self.params)
        else:
            # Here we take the centered differences
            # We could alternatively use thermodynamic properties
            # but these have not yet been implemented.
            dP = 100.
            dT = 0.01
            VdP0 = self.volume_function(pressure - dP/2., temperature,
                                        self.params)
            VdP1 = self.volume_function(pressure + dP/2., temperature,
                                        self.params)
            VdT0 = self.volume_function(pressure, temperature - dT/2.,
                                        self.params)
            VdT1 = self.volume_function(pressure, temperature + dT/2.,
                                        self.params)
            volume = (VdP0 + VdP1 + VdT0 + VdT1)/4.

            gradVPT = np.zeros((2, 3))
            gradVPT[:, 1:] = np.eye(2)
            gradVPT[:, 0] = [(VdP1 - VdP0)/dP, (VdT1 - VdT0)/dT]
            VPT_covariance = gradVPT.T.dot(PT_covariance).dot(gradVPT)
            return volume, VPT_covariance
