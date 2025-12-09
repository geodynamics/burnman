.. _ref-geotherms:

Geotherms
=========

Unlike the pressure, the temperature of the lower mantle is relatively unconstrained.
BurnMan provides lightweight :class:`burnman.classes.geotherm.Geotherm` and
:class:`burnman.classes.geotherm.GeothermFromPressures` classes
that store a geotherm defined by arrays of depths and temperatures or
pressures and temperatures respectively. Both classes have a `temperatures()` method
that takes an array or list of depths as input and outputs the corresponding
temperatures along the geotherm. In addition, the
:class:`burnman.classes.geotherm.GeothermFromPressures` class has a
`temperatures_from_pressures()` method that takes an array or list of pressures as
input. A number of geotherms from the literature are provided.
See :ref:`ref-api-geotherms` for a list.

Also implemented in BurnMan is an :class:`burnman.classes.geotherm.AdiabaticGeotherm`
which allows for the computation of adiabatic temperature profiles given a pressure to
depth conversion, a material, and an anchor temperature and pressure.
Again, see :ref:`ref-api-geotherms` for more details.
