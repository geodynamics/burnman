.. _ref-materials-mineral:

Mineral Classes
---------------

.. currentmodule:: burnman

The Mineral classes in BurnMan are intended to represent pure
phases with fixed chemical compositions and structures.
They form the foundation for more complex materials such as solutions
and multi-phase composites.

In this section, we describe the base :ref:`ref-mineral`
(:class:`Mineral`) class and its key features,
the :ref:`ref-combined-mineral` (:class:`CombinedMineral`) subclass,
which allows combining multiple minerals into a single effective phase,
the various implemented :ref:`ref-eos`, and :ref:`ref-property-modifiers`
that can be applied to adjust mineral properties.

.. toctree::
  :maxdepth: 1

  materials_02_mineral_01_base_class
  materials_02_mineral_02_combined_mineral
  materials_02_mineral_03_equations_of_state
  materials_02_mineral_04_property_modifiers
