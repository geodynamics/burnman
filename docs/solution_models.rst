Solution models
===============

SolidSolution objects in Burnman (type :class:`~burnman.SolidSolution`)
take one of several methods which define the properties of the solution.

.. inheritance-diagram::
  burnman.material_classes.solutionmodel.MechanicalSolution
  burnman.material_classes.solutionmodel.IdealSolution
  burnman.material_classes.solutionmodel.AsymmetricRegularSolution
  burnman.material_classes.solutionmodel.SymmetricRegularSolution
  burnman.material_classes.solutionmodel.SubregularSolution
  
Base class
----------
.. autoclass:: burnman.SolidSolution
  :noindex:

.. autoclass:: burnman.SolutionModel

Mechanical solution
-------------------

.. autoclass:: burnman.material_classes.solutionmodel.MechanicalSolution

Ideal solution
--------------

.. autoclass:: burnman.material_classes.solutionmodel.IdealSolution

Asymmetric regular solution
---------------------------

.. autoclass:: burnman.material_classes.solutionmodel.AsymmetricRegularSolution

Symmetric regular solution
--------------------------

.. autoclass:: burnman.material_classes.solutionmodel.SymmetricRegularSolution

Subregular solution
-------------------

.. autoclass:: burnman.material_classes.solutionmodel.SubregularSolution


Solution tools
==============

.. autofunction:: burnman.solutiontools.transform_solution_to_new_basis
