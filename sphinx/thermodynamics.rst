
Thermodynamics
==============

Burnman has a number of functions and classes which deal with the thermodynamics of single phases and aggregates.

Lattice Vibrations
------------------

Debye model
^^^^^^^^^^^
.. automodule:: burnman.debye

Einstein model
^^^^^^^^^^^^^^
.. automodule:: burnman.eos.einstein

Solution models
---------------

.. automodule:: burnman.solutionmodel


Chemistry parsing
-----------------

.. autofunction:: burnman.processchemistry.read_masses
.. autofunction:: burnman.processchemistry.dictionarize_formula
.. autofunction:: burnman.processchemistry.formula_mass
.. autofunction:: burnman.processchemistry.dictionarize_site_formula
.. autofunction:: burnman.processchemistry.process_solution_chemistry
.. autofunction:: burnman.processchemistry.compositional_array
.. autofunction:: burnman.processchemistry.ordered_compositional_array


Chemical potentials
-------------------

.. autofunction:: burnman.chemicalpotentials.chemical_potentials
.. autofunction:: burnman.chemicalpotentials.fugacity
.. autofunction:: burnman.chemicalpotentials.relative_fugacity
