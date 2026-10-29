.. _ref-fitting-solutions:

Fitting Solution Model Parameters to Data
-----------------------------------------

Overview
^^^^^^^^

BurnMan implements completely general nonlinear least squares fitting routines
in the function :func:`burnman.optimize.eos_fitting.nonlinear_least_squares_fit`.
This function takes an instance of a :class:`burnman.optimize.nonlinear_fitting.NonLinearModel`
as input, along with data to be fit, and returns best-fit parameters along with
uncertainties.

For many users, the process of constructing the `NonLinearModel` class may be
unfamiliar and tedious. To facilitate fitting solution parameters to data, BurnMan
provides a helper function that constructs `NonLinearModel` instances
so that users can easily fit solution parameters to data.

Available Fitting Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: burnman.optimize.eos_fitting.fit_XPTp_data
    :no-index:
