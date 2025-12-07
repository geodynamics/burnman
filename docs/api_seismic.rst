.. _ref-api-seismic:

Seismic Models
==============

.. contents::
   :local:
   :class: this-will-duplicate-information-and-it-is-still-useful-here
   :depth: 2

Base class for all seismic models
---------------------------------
.. autoclass:: burnman.classes.seismic.Seismic1DModel


Class for 1D Models
-------------------
.. autoclass:: burnman.classes.seismic.SeismicTable



Models currently implemented
----------------------------

.. autoclass:: burnman.classes.seismic.PREM

.. autoclass:: burnman.classes.seismic.Slow

.. autoclass:: burnman.classes.seismic.Fast

.. autoclass:: burnman.classes.seismic.STW105

.. autoclass:: burnman.classes.seismic.IASP91

.. autoclass:: burnman.classes.seismic.AK135


Attenuation Correction
-----------------------

.. autofunction:: burnman.tools.seismic.attenuation_correction
