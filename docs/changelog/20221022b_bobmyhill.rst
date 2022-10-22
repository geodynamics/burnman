* BurnMan now has a generalised PolynomialSolution class. The
  non-ideal excesses in this model are polynomial functions
  of composition. The coefficients of these polynomials are either
  constant internal energy, entropy and volume
  (i.e. Gxs = (Eijk... - T Sijk... + P Vijk...) xi xj xk...) or a
  linear sum of Mineral objects
  (i.e. Gxs = amijk mineralm.gibbs ... xi xj xk ...).
  Coefficients for both of these forms are passed as a list of lists
  to instantiate a PolynomialSolution object, which can then be passed
  as usual to instantiate a Solution object.

  Optionally, a transformation matrix can be passed that 
  allows users to define the coefficients above as functions of
  a modified set of endmember proportions: p'i = Aij pj. This
  may be useful when dealing with ordering, as expressing
  excess terms as a function of order parameter is often much
  more illuminating than expressing them in terms of endmember
  proportions.

  This class can deal with arbitrarily high powers in endmember
  proportions. However, because the class internally converts the
  list of lists to numpy arrays, high powers of solutions with a
  large number of endmembers will create very large
  arrays (with order n_endmembers^(highest power) elements). 
  This may significantly slow down calculations.

  For most purposes, using dense numpy arrays is much faster than
  using sparse arrays in COO format (a downside of using python).
  Should users need particularly complex solutions with high power
  terms, a modified solution model will be required. 

  *Bob Myhill, 2022/10/22*
