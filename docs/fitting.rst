.. _ref-fitting:

Fitting Model Parameters to Data
================================

BurnMan provides several tools for fitting model parameters to experimental or observational data.
These include:

- calculation of phase proportions in an assemblage given a bulk composition and phase compositions
- calculation of independent endmember proportions in a solution given its composition
- fitting equations of state parameters to P-T data including volumes,
    heat capacities and seismic velocities.
- fitting solution model parameters to P-T-compositional data, also including volumes,
    heat capacities and seismic velocities.

In all cases, the optimal solution is found in a weighted least-squares sense, taking into account
uncertainties in the data where provided, and propogating uncertainties into the fitted parameters.
Data uncertainties can include full covariance matrices in P, T, V, composition space, etc.
Prior estimates of parameters can also be provided.

.. toctree::
   :maxdepth: 2

   fitting_01_compositions
   fitting_02_equations_of_state
   fitting_03_solutions