.. _ref-geotherms:

Geotherms
=========

Unlike the pressure, the temperature of the lower mantle is relatively unconstrained.
As elsewhere, BurnMan provides a number of built-in geotherms, as well as the ability to use user-defined temperature-depth relationships.
A geotherm in BurnMan is an object that returns temperature as a function of pressure.
Alternatively, the user could ignore the geothermal and compute elastic velocities for a range of temperatures at any give pressure.

Currently, we include geotherms published by :cite:`Brown1981` and :cite:`anderson1982earth`.
Alternatively one can use an adiabatic gradient defined by the thermoelastic properties of a given mineralogical model.
For a homogeneous material, the adiabatic temperature profile is given by integrating the ordinary differential equation (ODE)

.. math::
    \left(\frac{\text{d}T}{\text{d}P}\right)_S = \frac{\gamma T}{K_S}.
    :label: geoth

This equation can be extended to multiphase composite using the first law of thermodynamics to arrive at

.. math::
    \left(\frac{\text{d}T}{\text{d}P}\right)_S = \frac{ T \displaystyle\sum_{i} \frac{ n_i C_{Pi} \gamma_i }{K_{Si}}}{ \displaystyle\sum_{i} n_i C_{Pi} },
    :label: geoth2

where the subscripts correspond to the :math:`i` th phase, :math:`C_P` is the heat capacity at constant pressure of a phase, and the other symbols are as defined above.
Integrating this ODE requires a choice in anchor temperature (:math:`T_0`) at the top of the lower mantle (or including this as a parameter in an inversion).
As the adiabatic geotherm is dependent on the thermoelastic parameters at high pressures and temperatures, it is dependent on the equation of state used.
