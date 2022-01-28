# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2022 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_calibrants
------------------

This example demonstrates the use of BurnMan's library of pressure calibrants.
These calibrants are stripped-down versions of BurnMan's minerals, in that
they are only designed to return pressure as a function of
volume and temperature or volume as a function of pressure and temperature.

*Uses:*

* :class:`burnman.classes.calibrant.Calibrant`
* :func:`burnman.calibrants.tools.pressure_to_pressure`
* Decker (1971) calibration for NaCl

*Demonstrates:*

* Use of the Calibrant class
* Conversion between pressures given by two different calibrations.
"""

import numpy as np
from burnman.calibrants import decker_1971
from burnman.calibrants.tools import pressure_to_pressure
from burnman.tools.unitcell import molar_volume_from_unit_cell_volume


# In this example, we'll be using the Decker (1971)
# equation of state for NaCl.
# The volume at any pressure and temperature can be calculated like so:
V0 = decker_1971.B1_NaCl.volume(pressure=1.e5,
                                temperature=293.15)

# Bear in mind that volumes in BurnMan are usually calculated in molar units.
# To convert from unit cell volumes in A^3/unit cell to molar volumes
# in m^3/mol, you can use the BurnMan tool
# "molar_volume_from_unit_cell_volume", e.g.:
V0_std = molar_volume_from_unit_cell_volume(5.6404**3.,
                                            decker_1971.B1_NaCl.params['Z'])


# For the next part of the example,
# we'll use the volume and temperature given in the last line
# of Table II in Decker (1971):
volume = V0*(1.-0.2950)  # m^3/mol
temperature = 1073.15  # K

# The Calibrant class in BurnMan also allows users to propagate
# uncertainties in volume and temperature into a full
# variance-covariance matrix.

# Let's assume the uncertainty in volume is 1%,
# the uncertainty in temperature is 1 K,
# and the uncertainties are uncorrelated.
sigma_volume = volume / 100.
sigma_temperature = 10.

print('Experimental observations:')
print(f'Volume: {volume*1.e6:.4f} +/- {sigma_volume*1.e6:.4f} m^3/mol')
print(f'Temperature: {temperature} +/- {sigma_temperature} K')
print()

# The calculated pressure can be calculated like so:
pressure = decker_1971.B1_NaCl.pressure(volume, temperature)

print('Pressure estimated using the Decker equation of state:')
print(f'{pressure/1.e9:.4f} GPa')
print()

# The above command didn't do any error propagation.
# Let's utilise the reported uncertainties in the volume and
# temperature data to output a covariance matrix
var_VT = np.array([[sigma_volume**2., 0.],
                   [0., sigma_temperature**2.]])


pressure, var_PVT = decker_1971.B1_NaCl.pressure(volume, temperature, var_VT)

# Standard conversion from the covariance matrix to the
# standard deviations and correlation matrix
sigma_PVT = np.sqrt(np.diag(var_PVT))
corr_PVT = np.einsum('i, ij, j -> ij', 1./sigma_PVT, var_PVT, 1./sigma_PVT)

print('Pressure with propagated uncertainty\n'
      'estimated using the Decker equation of state:')
print(f'{pressure/1.e9:.4f} +/- {sigma_PVT[0]/1.e9:.4f} GPa')
print('PVT correlation matrix:')
print(f'{corr_PVT}')
print('The uncertainties and correlation matrix only take into account ')
print('the uncertainties in the measured volume and temperature, ')
print('not the uncertainties in the calibrated equation of state.')
print()


# We can also check that the inverse operation (pressure -> volume)
# produces the original covariance matrix is consistent
var_PT = var_PVT[np.ix_([0, 2], [0, 2])]
volume, var_VPT = decker_1971.B1_NaCl.volume(pressure, temperature, var_PT)

sigma_VPT = np.sqrt(np.diag(var_VPT))
corr_VPT = np.einsum('i, ij, j -> ij', 1./sigma_VPT, var_VPT, 1./sigma_VPT)

print('Consistency checks:')
print(f'Volume: {volume*1.e6:.4f} +/- {sigma_VPT[0]*1.e6:.4f} m^3/mol')
print(f'Temperature: {temperature} +/- {sigma_VPT[2]:.4f} K')
print(f'V-T correlation: {corr_VPT[0, 2]:.6f}')
print()

# Finally, let's convert from one calibration to another
# Here, we convert from the Decker NaCl equation of state
# to the same equation of state, so we should see no change
# in the estimated pressure and covariance matrix.
new_P, new_var_PT = pressure_to_pressure(decker_1971.B1_NaCl,
                                         decker_1971.B1_NaCl,
                                         pressure, temperature, var_PT)

print('Result of conversion from Decker calibration to Decker calibration:')
print(f'{new_P/1.e9:.4f} +/- {np.sqrt(new_var_PT[0, 0])/1.e9:.4f} GPa')
print('(This should be the same as the pressures above)')
