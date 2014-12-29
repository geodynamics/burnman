
Materials
=========

Burnman operates on materials (type :class:`~burnman.material.Material`) most prominently in form of minerals (:class:`~burnman.mineral.Mineral`) and composites (:class:`~burnman.composite.Composite`).


.. inheritance-diagram:: burnman.Material burnman.Composite burnman.Mineral burnman.solidsolution.SolidSolution burnman.mineral_helpers.HelperSpinTransition


Base Material
-------------

.. autoclass:: burnman.material.Material


Minerals
----------------------------

.. autoclass:: burnman.mineral.Mineral

.. autoclass:: burnman.SolidSolution

.. autoclass:: burnman.mineral_helpers.HelperSpinTransition

.. autoclass:: burnman.mineral_helpers.HelperSolidSolution

.. autoclass:: burnman.mineral_helpers.HelperFeDependent


Composite
---------

.. autoclass:: burnman.composite.Composite
