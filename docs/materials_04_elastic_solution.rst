.. _ref-materials-elastic-solution:

Elastic Solution Class
----------------------

Elastic solution models in BurnMan are much like standard
:ref:`ref-materials-solution-models`. The only difference is
that instead of considering excess **Gibbs energy** of mixing when
all of the endmembers are evaluated at the same **pressure** and temperature,
elastic solution models consider excess **Helmholtz energy** of mixing when
all of the endmembers are evaluated at the same **volume** and temperature.

As a consequence of the different thermodynamic potential being used,
elastic solution models can include **excess contributions to pressure**,
rather than **excess contributions to volume**.

In BurnMan, we define :class:`burnman.ElasticSolution` to represent a solution phase. Each instance of this class
contains a :class:`burnman.ElasticSolutionModel` which defines how the Helmholtz energy of the solution is calculated
from its components. The following sections describe these classes in more detail.

.. toctree::
  :maxdepth: 2

  materials_04_elastic_solution_01_base_class
  materials_04_elastic_solution_02_solution_models
