import numpy as np


def pressure_to_pressure(
    old_calibrant, new_calibrant, pressure, temperature, PT_covariance=None
):
    """
    Convert from pressure defined by one calibrated equation of
    state of a material to pressure defined by an independent
    calibration of the same material.

    :param old_calibrant: The original calibration used to estimate the pressure
    :type old_calibrant: :class:`burnman.Calibrant`

    :param new_calibrant: The new calibration from which the pressure is desired
    :type new_calibrant: :class:`burnman.Calibrant`

    :param pressure: The pressure calculated using the old calibration
    :type pressure: float

    :param temperature: The temperature of the material
    :type temperature: float

    :param PT_covariance: The pressure-temperature variance-covariance matrix
    :type PT_covariance: 2x2 numpy.array [optional]

    :returns: The pressure of the calibrant [Pa] (float) and a 2x2 numpy array
        (if the PT_covariance is provided) containing the
        pressure-temperature variance-covariance matrix.

    :rtype: tuple
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
