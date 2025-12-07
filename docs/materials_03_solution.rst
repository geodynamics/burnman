.. _ref-materials-solution:

Solutions
---------

Many phases (whether minerals, melts, fluids or gases) can exist over a finite region of composition space.
These spaces are bounded by endmembers (which may themselves not be stable), and each phase can then be described
as a solution of those endmembers. 

In a mechanically and thermally equilibrated mechanical mixture (i.e., not a solution),
the Gibbs energy of the mixture is simply the sum of the Gibbs energies of the individual components.
However, in a solution, the Gibbs energy also includes an entropic contribution from the mixing of
the different components. This entropic contribution generally lowers the overall Gibbs energy of the phase,
making the solution more stable than a mechanical mixture of its endmembers. In addition, there may be
enthalpic contributions to the Gibbs energy of mixing, which can either further stabilize or destabilize the
solution relative to the mechanical mixture, such that:

.. math::
    \mathcal{G}_{\textrm{solution}} = \sum_i n_i \mathcal{G}_i + \Delta \mathcal{G}_{\textrm{mixing}} \\
    \Delta \mathcal{G}_{\textrm{mixing}} = \Delta \mathcal{H}_{\textrm{mixing}} - T \Delta S_{\textrm{mixing}}

where :math:`n_i` is the number of moles of component :math:`i`,
:math:`\mathcal{G}_i` is the molar Gibbs energy of component :math:`i`,
:math:`\Delta \mathcal{H}_{\textrm{mixing}}` is the enthalpy of mixing,
and :math:`\Delta S_{\textrm{mixing}}` is the entropy of mixing.

It is important to recognise that the components of a solution must be compatible with the structure of the phase;
that is, they must all be contained within the same crystallographic framework. They do not need to have unique
compositions. For example, we can model a displacive phase transition through a solution of two identical components
which differ only in the element of symmetry lost during the transition, and we can model order-disorder transitions
through solutions of components which differ only in the distribution of cations over crystallographic sites.

Derivatives of the Gibbs energy of the solution can now be defined with respect to amounts of components,
in addition to pressure and temperature. The first derivatives give rise to so-called chemical potentials,
which are extremely useful in defining phase equilibria.

In BurnMan, we define :class:`burnman.Solution` to represent a solution phase. Each instance of this class
contains a :class:`burnman.SolutionModel` which defines how the Gibbs energy of the solution is calculated
from its components. The following sections describe these classes in more detail.


.. toctree::
  :maxdepth: 2

  materials_03_solution_01_base_class
  materials_03_solution_02_solution_models
