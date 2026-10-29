2022-10-03
^^^^^^^^^^

* BurnMan now has two new helper functions:
  :func:`burnman.tools.chemistry.reactions_from_stoichiometric_matrix` 
  :func:`burnman.tools.chemistry.reactions_from_formulae`.
  These functions generate a complete list of reactions
  (forward and reverse) from either the stoichiometric matrix
  (a 2D numpy array containing the amount of component j in phase i),
  or from a list of formulae
  (as strings or dictionaries of elemental amounts).

  *Bob Myhill, 2022/10/03*
