* BurnMan now has a Calibrant class. This class is a stripped-down
  version of a Material that is initialised with a params dictionary
  and a function to compute the pressure or volume at given conditions.

  Objects derived from this class have the ability to output
  pressure/volume as a function of volume/pressure and temperature.
  The user can pass a V-T or P-T covariance matrix as an additional
  optional argument, in which case BurnMan will propagate the errors
  and output a full P-V-T or V-P-T covariance matrix.

  A helper function is provided to convert from pressure calculated
  using one calibration into pressure calculated using another
  calibration of the same material.

  This class is expected to be used by experimental petrologists and
  those interpreting high pressure experimental data.

  *Bob Myhill, 2022/01/28*
