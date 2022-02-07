Solution models
===============

Solution objects in Burnman (type :class:`~burnman.Solution`, alias :class:`~burnman.SolidSolution`)
take one of several methods which define the properties of the solution.

.. inheritance-diagram::
  burnman.classes.solutionmodel.MechanicalSolution
  burnman.classes.solutionmodel.IdealSolution
  burnman.classes.solutionmodel.AsymmetricRegularSolution
  burnman.classes.solutionmodel.SymmetricRegularSolution
  burnman.classes.solutionmodel.SubregularSolution

Base class
----------
.. autoclass:: burnman.Solution
  :noindex:

.. autoclass:: burnman.SolutionModel

Mechanical solution
-------------------

.. autoclass:: burnman.classes.solutionmodel.MechanicalSolution

Ideal solution
--------------

.. autoclass:: burnman.classes.solutionmodel.IdealSolution

Asymmetric regular solution
---------------------------

.. autoclass:: burnman.classes.solutionmodel.AsymmetricRegularSolution

Symmetric regular solution
--------------------------

.. autoclass:: burnman.classes.solutionmodel.SymmetricRegularSolution

Subregular solution
-------------------

.. autoclass:: burnman.classes.solutionmodel.SubregularSolution


Solution tools
==============

.. autofunction:: burnman.tools.solution.transform_solution_to_new_basis
