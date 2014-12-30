
Materials
=========

Burnman operates on materials (type :class:`~burnman.material.Material`) most prominently in form of minerals (:class:`~burnman.mineral.Mineral`) and composites (:class:`~burnman.composite.Composite`).


.. inheritance-diagram:: burnman.Material burnman.Composite burnman.Mineral burnman.solidsolution.SolidSolution burnman.mineral_helpers.HelperSpinTransition


.. inheritance-diagram:: burnman.Material burnman.Composite burnman.Mineral burnman.solutionmodel.SolutionModel burnman.solidsolution.SolidSolution burnman.mineral_helpers.HelperSpinTransition


Base Classes
------------

Materials
^^^^^^^^^

.. autoclass:: burnman.material.Material

Solution models
^^^^^^^^^^^^^^^

.. autoclass:: burnman.solutionmodel.SolutionModel


Minerals
--------

Endmembers
^^^^^^^^^^

.. autoclass:: burnman.mineral.Mineral

Solid solutions
^^^^^^^^^^^^^^^

.. autoclass:: burnman.solidsolution.SolidSolution

Solid solution helpers (deprecated)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: burnman.mineral_helpers.HelperSpinTransition

.. autoclass:: burnman.mineral_helpers.HelperSolidSolution

.. autoclass:: burnman.mineral_helpers.HelperFeDependent


Composites
----------

.. autoclass:: burnman.composite.Composite

