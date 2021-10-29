* The :class:`burnman.Composite` class now has a `chemical_potential` method.
  This replaces the more limited function
  :func:`burnman.tools.chemistry.chemical_potentials`. The related functions
  :func:`burnman.tools.chemistry.fugacity` and
  :func:`burnman.tools.chemistry.relative_fugacity` have been updated
  to use this method, and to take Composites as arguments, rather than
  lists of Minerals.

  *Bob Myhill, 2021/10/30*
