
Materials
=========

Burnman operates on materials (type :class:`~burnman.Material`)
most prominently in the form of minerals
(:class:`~burnman.Mineral`) and composites (:class:`~burnman.Composite`).


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

Solutions
^^^^^^^^^

.. autoclass:: burnman.Solution

.. autoclass:: burnman.SolidSolution
  :no-inherited-members:

.. autoclass:: burnman.ElasticSolution

.. autoclass:: burnman.ElasticSolidSolution
  :no-inherited-members:

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
