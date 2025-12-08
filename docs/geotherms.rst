.. _ref-geotherms:

Geotherms
=========

Unlike the pressure, the temperature of the lower mantle is relatively unconstrained.
BurnMan provides a lightweight :class:`burnman.classes.geotherm.Geotherm` class to represent
temperature profiles as a function of depth. A number of geotherms from the literature are provided.
See :ref:`ref-api-geotherms` for a list.

Also implemented in BurnMan is an :class:`burnman.classes.geotherm.AdiabaticGeotherm`
which allows for the computation of adiabatic temperature profiles given a pressure to
depth conversion, a material, and an anchor temperature and pressure.
Again, see :ref:`ref-api-geotherms` for more details.