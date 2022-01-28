
Materials
=========

Burnman operates on materials (type :class:`~burnman.Material`)
most prominently in the form of minerals
(:class:`~burnman.Mineral`) and composites (:class:`~burnman.Composite`).


.. inheritance-diagram:: burnman.Material burnman.Mineral burnman.PerplexMaterial burnman.SolidSolution burnman.Composite burnman.CombinedMineral burnman.AnisotropicMaterial burnman.AnisotropicMineral


Material Base Class
-------------------

.. autoclass:: burnman.Material

Perple_X Class
--------------

.. autoclass:: burnman.PerplexMaterial

Minerals
--------

Endmembers
^^^^^^^^^^

.. autoclass:: burnman.Mineral

Solid solutions
^^^^^^^^^^^^^^^

.. autoclass:: burnman.SolidSolution

Mineral helpers
^^^^^^^^^^^^^^^

.. autoclass:: burnman.classes.mineral_helpers.HelperSpinTransition

Anisotropic materials
^^^^^^^^^^^^^^^^^^^^^

.. autoclass:: burnman.AnisotropicMaterial

.. autoclass:: burnman.AnisotropicMineral

.. autofunction:: burnman.cell_parameters_to_vectors

.. autofunction:: burnman.cell_vectors_to_parameters

Composites
----------

.. autoclass:: burnman.Composite


Calibrants
----------

.. autoclass:: burnman.Calibrant
