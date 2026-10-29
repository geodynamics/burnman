2022-02-05
^^^^^^^^^^

* New: BurnMan now has a new fitting function,
  :func:`burnman.optimize.eos_fitting.fit_XPTp_data`.
  It has the same arguments as
  :func:`burnman.optimize.eos_fitting.fit_PTp_data`, but instead of
  fitting data to a mineral of fixed composition,
  it is able to simultaneously fit parameters of solution model
  parameters, including equation of state parameters for each
  endmember.

  *Bob Myhill, 2022/02/05*
