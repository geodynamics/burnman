.. _ref-fitting-equations-of-state:

Fitting Equation of State Parameters to Data
--------------------------------------------

Overview
^^^^^^^^

BurnMan implements completely general nonlinear least squares fitting routines
in the function :func:`burnman.optimize.eos_fitting.nonlinear_least_squares_fit`.
This function takes an instance of a :class:`burnman.optimize.nonlinear_fitting.NonLinearModel`
as input, along with data to be fit, and returns best-fit parameters along with
uncertainties.

For many users, the process of constructing the `NonLinearModel` class may be
unfamiliar and tedious. To facilitate fitting equations of state to data, BurnMan
provides a number of helper functions that construct `NonLinearModel` instances
so that users can easily fit equations of state to data.

Available Fitting Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: burnman.optimize.eos_fitting.fit_PTV_data
    :no-index:

.. autofunction:: burnman.optimize.eos_fitting.fit_PTp_data
    :no-index:

.. autofunction:: burnman.optimize.eos_fitting.fit_VTP_data
    :no-index:

.. autofunction:: burnman.optimize.eos_fitting.fit_VTp_data
    :no-index: