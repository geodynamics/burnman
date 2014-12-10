Averaging Schemes
=================

Given a set of mineral physics parameters and an equation of state we can calculate the density, bulk, and shear modulus for a given phase. 
However, as soon as we have a composite material (e.g., a rock), the determination of elastic properties become more complicated. 
The bulk and shear modulus of a rock are dependent on the specific geometry of the grains in the rock, so there is no general 
formula for its averaged elastic properties.  Instead, we must choose from a number of averaging schemes if we want a single value, 
or use bounding methods to get a range of possible values.  The module :mod:`burnman.averaging_schemes` provides a number of different 
average and bounding schemes for determining a composite rock's physical parameters.

Base class
-----------
.. autoclass:: burnman.averaging_schemes.AveragingScheme

Voigt bound
-----------
.. autoclass:: burnman.averaging_schemes.Voigt

Reuss bound
-----------
.. autoclass:: burnman.averaging_schemes.Reuss

Voigt-Reuss-Hill average
------------------------
.. autoclass:: burnman.averaging_schemes.VoigtReussHill

Hashin-Shtrikman upper bound
----------------------------
.. autoclass:: burnman.averaging_schemes.HashinShtrikmanUpper

Hashin-Shtrikman lower bound
----------------------------
.. autoclass:: burnman.averaging_schemes.HashinShtrikmanLower

Hashin-Shtrikman arithmetic average
-----------------------------------
.. autoclass:: burnman.averaging_schemes.HashinShtrikmanAverage

