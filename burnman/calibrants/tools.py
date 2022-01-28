import numpy as np


def pressure_to_pressure(old_calibrant, new_calibrant,
                         pressure, temperature, PT_covariance=None):
    """
    Convert from pressure defined by one calibrated equation of
    state of a material to pressure defined by an independent
    calibration of the same material.


    Arguments
    ---------
    old_calibrant : Calibrant object
        The original calibration used to estimate the pressure
    new_calibrant : Calibrant object
        The new calibration from which the pressure is desired
    pressure : float
        The pressure calculated using the old calibration
    temperature : float
        The temperature of the material

    PT_covariance : 2x2 numpy array [optional]
            The pressure-temperature
            variance-covariance matrix


    Returns
    -------
    pressure : float
        The pressure of the calibrant [Pa]

    PT_covariance : 2x2 numpy array (if PT_covariance is provided)
        The pressure-temperature variance-covariance matrix.
    """

    if PT_covariance is None:
        V = old_calibrant.volume(pressure, temperature)
        P = new_calibrant.pressure(V, temperature)
        return P
    else:
        V, var_VPT = old_calibrant.volume(pressure, temperature, PT_covariance)
        VT_covariance = var_VPT[np.ix_([0, 2], [0, 2])]
        P, var_PVT = new_calibrant.pressure(V, temperature, VT_covariance)
        PT_covariance = var_PVT[np.ix_([0, 2], [0, 2])]
        return P, PT_covariance
