.. _ref-tools-thermobarometry:

Thermobarometry
---------------

Overview
^^^^^^^^

One of the most common applications of thermodynamic modeling in
geosciences is to estimate the pressure and temperature conditions at
which a given mineral assemblage is stable. This process is known as
thermobarometry. Burnman provides tools to perform thermobarometric
calculations using a least-squares inversion approach, following
and extending the Optimal Thermobarometry
method described by :cite:`Powell1994`.

The essential idea behind the method is to minimize a misfit scalar that
is a function of reaction affinities, weighted by a covariance matrix that
is a function of both measurement uncertainties and thermodynamic model
uncertainties. The misfit scalar is defined as:

.. math::
    \Phi = \sum_i \sum_j a_i W_{ij} a_j \\
    W_{ij} = (C_{ij} + M_{ij})^{-1}


where :math:`a_i` is the affinity of reaction :math:`i`, :math:`W_{ij}` is the
weighting matrix, :math:`C_{ij}` is the covariance matrix of measurement
uncertainties, and :math:`M_{ij}` is the covariance matrix of thermodynamic
model uncertainties.

The full documentation for the thermobarometry tools can be found at
:ref:`ref-api-tools-thermobarometry`. The documentation for the
Optimal Thermobarometry function is given below.


Implemented function
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: burnman.tools.thermobarometry.estimate_conditions
    :no-index: