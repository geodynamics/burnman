.. _ref-calibrants:

Calibrants
**********

Sometimes, it is unnecessary or inconvenient to define a full material
model to represent a substance. One of these cases is when dealing with
experimental calibrants, where we are only interested in finding the
pressure and/or temperature conditions based on measured unit cell volumes.
BurnMan includes a Calibrant class for this purpose. A Calibrant is
defined by a volumetric equation of state (EOS) and accompanying parameters.
It may or may not be temperature dependent.

See :ref:`ref-api-calibrant-database` for a list of built-in calibrants
available in BurnMan.
