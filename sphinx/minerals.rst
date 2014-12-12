
Minerals
========

Burnman operates on rocks (type :class:`~burnman.material.Material`) most prominently in form of minerals (:class:`~burnman.mineral.Mineral`) and composites (:class:`~burnman.composite.Composite`).


.. inheritance-diagram:: burnman.Material burnman.Composite burnman.Mineral burnman.mineral_helpers.HelperSolidSolution burnman.mineral_helpers.HelperSpinTransition burnman.mineral_helpers.HelperFeDependent


Base Class
----------

.. autoclass:: burnman.material.Material


Base for individual Minerals
----------------------------


.. autoclass:: burnman.mineral.Mineral

Composite
---------

.. autoclass:: burnman.composite.Composite





Mineral helpers
---------------

.. autoclass:: burnman.SolidSolution
.. autoclass:: burnman.SolutionModel
.. autoclass:: burnman.mineral_helpers.HelperSpinTransition
.. autoclass:: burnman.mineral_helpers.HelperSolidSolution
.. autoclass:: burnman.mineral_helpers.HelperFeDependent
