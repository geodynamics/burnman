* BurnMan now has a new anisotropic equation of state class,
  :class:`burnman.AnisotropicMineral`, which can be used to model materials of
  arbitrary symmetry under hydrostatic conditions.
  Users can set the state (pressure and temperature)
  of AnisotropicMineral objects and then retrieve their anisotropic properties.
  Details of the formulation can be found in :cite:`Myhill2022`.
  Examples are provided in the file examples/example\_anisotropic\_mineral.py.

  *Bob Myhill, 2021/10/03*
