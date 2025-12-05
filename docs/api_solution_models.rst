Solution models
===============

Solution objects in BurnMan are instances of one of three classes:
type :class:`~burnman.Solution` (alias :class:`~burnman.SolidSolution`),
type :class:`~burnman.ElasticSolution` and
type :class:`~burnman.AnisotropicSolution`.

The :class:`~burnman.Solution` class implements commonly used models
(in petrology). Excess properties are defined relative to the endmember
properties at fixed pressure and temperature.
The formulations are defined with interaction parameters
such as excess energies, volumes and entropies.

The :class:`~burnman.AnisotropicSolution` class extends the
:class:`~burnman.Solution` class to include anisotropic elastic properties.

The standard and anisotropic solution classes can also be "relaxed" via
two derived classes: :class:`~burnman.RelaxedSolution` and
:class:`~burnman.RelaxedAnisotropicSolution`. In these classes, specific
linear combinations of molar fractions can be designated as relaxation vectors.
These are allowed to vary freely to minimize the Gibbs free energy of the solution.

An independent class of solution models is provided by the
:class:`~burnman.ElasticSolution` class. This class instead defines excess
properties are relative to the endmember properties at fixed volume and
temperature. Such models have their roots in atom-scale considerations;
mixing of two instances of the same lattice type requires deformation
(local lattice distortions), that can be considered to induce
a local chemical stress. Therefore, volume may be a more useful
independent variable than pressure. For more details, see :cite:`Myhill2018`.

.. rubric:: Available solution classes
.. inheritance-diagram::
  burnman.Solution
  burnman.RelaxedSolution
  burnman.ElasticSolution
  burnman.AnisotropicSolution
  burnman.RelaxedAnisotropicSolution

All of the above classes accept several models which define the properties
of the solution. The standard and anisotropic solution classes share the same
solution models, while the elastic solution class has its own set of
closely related solution models.

.. rubric:: Available solution model classes
.. inheritance-diagram::
  burnman.classes.solutionmodel.MechanicalSolution
  burnman.classes.solutionmodel.IdealSolution
  burnman.classes.solutionmodel.AsymmetricRegularSolution
  burnman.classes.solutionmodel.SymmetricRegularSolution
  burnman.classes.solutionmodel.SubregularSolution
  burnman.classes.solutionmodel.FunctionSolution

.. rubric:: Available elastic solution model classes
.. inheritance-diagram::
  burnman.classes.elasticsolutionmodel.ElasticMechanicalSolution
  burnman.classes.elasticsolutionmodel.ElasticIdealSolution
  burnman.classes.elasticsolutionmodel.ElasticAsymmetricRegularSolution
  burnman.classes.elasticsolutionmodel.ElasticSymmetricRegularSolution
  burnman.classes.elasticsolutionmodel.ElasticSubregularSolution
  burnman.classes.elasticsolutionmodel.ElasticFunctionSolution

Base classes
------------
.. autoclass:: burnman.Solution
  :noindex:

.. autoclass:: burnman.ElasticSolution
  :noindex:

.. autoclass:: burnman.AnisotropicSolution

.. autoclass:: burnman.SolutionModel

.. autoclass:: burnman.ElasticSolutionModel

Relaxed solutions
-----------------

.. autoclass:: burnman.RelaxedSolution

.. autoclass:: burnman.RelaxedAnisotropicSolution

Mechanical solution
-------------------

.. autoclass:: burnman.classes.solutionmodel.MechanicalSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticMechanicalSolution

Ideal solution
--------------

.. autoclass:: burnman.classes.solutionmodel.IdealSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticIdealSolution

Asymmetric regular solution
---------------------------

.. autoclass:: burnman.classes.solutionmodel.AsymmetricRegularSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticAsymmetricRegularSolution

Symmetric regular solution
--------------------------

.. autoclass:: burnman.classes.solutionmodel.SymmetricRegularSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticSymmetricRegularSolution

Subregular solution
-------------------

.. autoclass:: burnman.classes.solutionmodel.SubregularSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticSubregularSolution

Function solution
-----------------

.. autoclass:: burnman.classes.solutionmodel.FunctionSolution

.. autoclass:: burnman.classes.elasticsolutionmodel.ElasticFunctionSolution


Solution tools
==============

.. autofunction:: burnman.tools.solution.transform_solution_to_new_basis
