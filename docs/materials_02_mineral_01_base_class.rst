.. _ref-mineral:

Mineral Base Class
~~~~~~~~~~~~~~~~~~

.. currentmodule:: burnman

The :class:`Mineral` class is the **base class** for all mineral phases
in BurnMan. It derives from the
:ref:`ref-materials-base-class` (:class:`Material`), using an
equation of state to implement functionality specific to
pure mineral phases with fixed chemical compositions
and crystal structures.

Key Features
^^^^^^^^^^^^

* **Equations of State**: Minerals use a variety of implemented
  equations of state (EOS) to model their thermodynamic behaviour
  under varying pressure and temperature conditions. The EOS is
  selected at instantiation and governs property calculations.

* **Property Modifiers**: Minerals can apply property modifiers
  to adjust baseline properties from the EOS. This allows
  incorporation of experimental corrections or additional
  physical effects.

* **Composition**: Each :class:`Mineral` instance has a fixed
  chemical composition defined at creation. This composition
  influences density, molar volume, and other derived properties.

Mineral Datasets
^^^^^^^^^^^^^^^^

BurnMan includes a database of predefined mineral phases
with commonly used equations of state and compositions.
See :ref:`ref-mineral-datasets` for a list of the available
mineral models.

Instantiation
^^^^^^^^^^^^^

Mineral objects are instantiated via a constructor that takes
a dictionary of parameters as argument. The parameters required
depend on the selected equation of state. Instantiation looks
like this:

.. code-block:: python

    from burnman import Mineral

    params = {
        'equation_of_state': 'slb3',  # using the Stixrude and Lithgow-Bertelloni EoS
        # EOS-specific parameters...
    }

    property_modifiers = [
        ['linear', {'delta_E': 1.e3, 'delta_S': 0., 'delta_V': 0.}],
        # Additional modifiers...
    ]

    mineral = Mineral(params, property_modifiers)

An important parameter in the `params` dictionary is the
`equation_of_state` key, which specifies the EOS to use.
This can either be a string referring to a built-in EOS,
or a custom EOS object implementing the required interface.
Available equations of state, and their string identifiers,
are described in :ref:`ref-eos`.

Property modifiers are optional and can be omitted. See
:ref:`ref-property-modifiers` for details.