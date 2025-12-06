.. _ref-materials-solution:

Solution Classes
----------------

Many phases (whether minerals, melts, fluids or gases) can exist over a finite region of composition space.
These spaces are bounded by endmembers (which may themselves not be stable), and each phase can then be described
as a solution of those endmembers. In a solid solution, different elements substitute for one another on distinct
crystallographic sites in the structure. For example, low pressure silicate garnets have two distinct sites on
which mixing takes place; a dodecahedral site (of which there are three per unit cell on an eight-cation basis)
and octahedral site (of which there are two per unit cell). A third tetrahedral cation site (three per unit cell)
is usually assumed to be occupied solely by silicon, and therefore can be ignored in solution calculations.
The chemical formula of many low pressure garnets exist within the solution:

.. math::
    \textrm{[Mg,Fe,Mn,Ca]}_3\textrm{[Al,Fe,Cr]}_2\textrm{Si}_3\textrm{O}_{12}


We typically calculate solution properties by appropriate differentiation of the Gibbs energy, where

.. math::
    \mathcal{G} = \sum_i n_i \left( \mathcal{G}_i + RT \ln \alpha_i \right)\\
    \alpha_i = \gamma_i \alpha_{\textrm{ideal}, i}


.. toctree::
  :maxdepth: 2

  materials_03_solution_01_base_class
  materials_03_solution_02_solution_models
