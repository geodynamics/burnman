.. _ref-materials-base-class:

Material Base Class
-------------------


Overview
~~~~~~~~

.. currentmodule:: burnman

The :class:`Material` class serves as the **base class** for all materials
in BurnMan, including :ref:`ref-materials-mineral` (:class:`Mineral`),
:ref:`ref-materials-solution` (:class:`Solution`),
and :ref:`ref-materials-composite` (:class:`Composite`).
It defines the common interface through
which all thermodynamic and elastic properties are computed, as well as the
workflow for setting state variables and selecting thermodynamic models.

Subclasses implement the specific physical behaviour (equations of state,
solution models, mechanical mixing laws), while :class:`Material` ensures
that all materials (simple or composite) expose a consistent and unified interface.

Regardless of type, materials support the same set of state variables and
property queries.

Typical usage follows this pattern:

1. **Instantiate** a material object.
2. **Set the state** using one or more of:
   
   * :meth:`Material.set_state` — set pressure and temperature,
   * :attr:`Material.number_of_moles` — specify the material amount.

3. **Query thermodynamic properties**, such as density, heat capacities,
   expansivity, elastic parameters, or seismic velocities.

Thermodynamic properties are computed on demand, with caching to optimise
performance. Changing any state variable automatically resets cached values
to ensure consistency.

Property Evaluation
~~~~~~~~~~~~~~~~~~~

A large number of physical properties can be evaluated. Examples include:

* Pressure (:attr:`pressure`)
* Temperature (:attr:`temperature`)
* Number of moles (:attr:`number_of_moles`)
* Mass (:attr:`mass`)
* Gibbs energy (:attr:`gibbs`)
* Helmholtz energy (:attr:`helmholtz`)
* Enthalpy (:attr:`enthalpy`)
* Internal energy (:attr:`internal_energy`)
* Entropy (:attr:`entropy`)
* Density (:attr:`density`)
* Volume (:attr:`volume`) and molar volume (:attr:`molar_volume`)
* Heat capacities (:attr:`molar_heat_capacity_p` and :attr:`molar_heat_capacity_v`)
* Thermal expansivity (:attr:`thermal_expansivity`)
* Gruneisen parameter (:attr:`gruneisen_parameter`) 
* Isothermal bulk modulus (Reuss) (:attr:`isothermal_bulk_modulus_reuss`)
* Isentropic bulk modulus (Reuss) (:attr:`isentropic_bulk_modulus_reuss`)
* Corresponding compressibilities (:attr:`isothermal_compressibility_reuss`,
  :attr:`isentropic_compressibility_reuss`)
* Shear modulus (:attr:`shear_modulus`)
* Effective seismic bulk and shear moduli (:attr:`effective_isentropic_bulk_modulus`,
  :attr:`effective_shear_modulus`), see :ref:`ref-averaging-schemes`.
* Seismic velocities (:attr:`p_wave_velocity`, :attr:`shear_wave_velocity`, :attr:`bulk_sound_velocity`), see :ref:`ref-averaging-schemes`.

Molar quantities are available via ``molar_``-prefixed properties. Aliases
exist for common properties (e.g., :attr:`C_p` for :attr:`heat_capacity_p`).

Properties generally require that the state (pressure and temperature) has
been set before computation.


State Management
~~~~~~~~~~~~~~~~

.. method:: Material.set_state(pressure, temperature)
    :noindex:

    Sets the thermodynamic state of the material by updating its **pressure** (Pa) and
    **temperature** (K). 


.. method:: Material.set_state_with_volume(volume, temperature)
    :noindex:

    Sets the state of the material using **volume** (:math:`\mathrm{m}^3`)
    rather than **pressure** (Pa).
    The method internally performs a pressure inversion: it solves for the
    pressure that produces the specified volume at the given temperature.


Introspection and Debugging
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. method:: Material.to_string()
    :noindex:

    Returns a human-readable string describing the material. Subclasses typically
    override this to include formula, composition, or model information.


.. method:: Material.debug_print()
    :noindex:

    Prints detailed diagnostic information about the material, including its
    current state, cached properties, and any internal solver or model details
    relevant to the subclass. Intended primarily for debugging complex materials
    such as solutions or composites.


.. method:: Material.unroll()
    :noindex:

    Expands the material into a flat list of :class:`burnman.Mineral` objects and
    their molar fractions.


.. method:: Material.print_minerals_of_current_state()
    :noindex:

    Prints a formatted list of all constituent minerals after applying the
    current pressure, temperature, and composition. This method invokes the same
    internal expansion logic as :meth:`Material.unroll` and is useful for
    inspecting the active mineral assemblage within solutions or composites.

Evaluation Routines
~~~~~~~~~~~~~~~~~~~

.. method:: Material.evaluate(property_list, pressures, temperatures, molar_fractions=None)
    :noindex:

    Computes and returns a thermodynamic property or properties at the requested conditions.
    The list of properties should be provided as a list of strings in ``property_list``.
    The conditions are specified via arrays which may be of arbitrary shape, with broadcasting
    applied as needed. Pressure and temperature arrays must have the same shape.


.. method:: Material.evaluate_with_volumes(property_name, volumes, temperatures, molar_fractions=None)
    :noindex:
   
    Evaluates a specified property using externally supplied **volume** values
    rather than pressures. This method performs a pressure inversion internally to
    determine the corresponding pressures for the given volumes at the specified temperatures.
