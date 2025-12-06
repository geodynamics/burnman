2021-08-05
^^^^^^^^^^

* :class:`burnman.composite.Composite` now has new properties which include
  the *endmember_formulae*, a list of *elements* which make up those formulae,
  the *stoichiometric_matrix* (number of atoms of element j in formula i),
  and an independent *reaction_basis*. These properties are calculated once
  when they are first needed, and then cached.

  *Bob Myhill, 2021/08/05*
