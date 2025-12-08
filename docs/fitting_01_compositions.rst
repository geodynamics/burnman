.. _ref-fitting-compositions:

Fitting Model Compositions to Data
----------------------------------

Overview
^^^^^^^^

In experimental petrology and geochemistry, it is often useful to:

- fit proportions of mineral phases to bulk rock compositions.
    This can be used, for example, to perform mass balance calculations
    to check if an experimental run product is consistent with the starting composition.
- fit molar fractions of solution endmembers to analysed elemental or oxide compositions.
    This can form a starting point for :ref:`ref-tools-thermobarometry`,
    or for fitting parameters of thermodynamic datasets. 

BurnMan provides functions to perform these types of composition fitting.

Available Fitting Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: burnman.composition_fitting.fit_phase_proportions_to_bulk_composition
    :no-index:

.. autofunction:: burnman.composition_fitting.fit_composition_to_solution
    :no-index:
