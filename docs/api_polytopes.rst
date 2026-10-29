.. _ref-api-polytopes:

Polytopes
=========

Often in mineral physics, solutions are subject to a set of linear constraints.
For example, the set of valid site-occupancies in solution models are
constrained by positivity and fixed sum constraints
(the amount of each chemical species must be greater than or equal to zero,
the sites in the structure are present in fixed ratios).
Similarly, the phase amounts in a composite must be more than or
equal to zero, and the compositions of each phase must sum to the
bulk composition of the composite.

Geometrically, linear equality and inequality constraints can be visualised
as a polytope (an n-dimensional) polyhedron. There are several situations
where it is convenient to be able to interrogate such objects to understand
the space of validity. In BurnMan, we make use of the module pycddlib to
create polytope objects. We also provide a number of tools for common
chemically-relevant operations.

Base class
----------

.. autoclass:: burnman.MaterialPolytope

Polytope tools
--------------

.. automodule:: burnman.tools.polytope
    :noindex:
