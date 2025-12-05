.. _ref-anisotropic-materials:

Anisotropic Materials
*********************

In addition to BurnMan's standard materials, which have
properties that are a function of pressure, temperature,
and composition, BurnMan also includes classes for
anisotropic elastic materials. These materials have elastic
properties that depend on stress rather than pressure,
or strain rather than volumetric compression.

Because most natural materials in the Earth can only withstand
small deviatoric stresses before yielding or fracturing,
the anisotropic materials in BurnMan have properties that are
defined as a function of pressure, but these properties include
those that describe the response of the material to small
deviatoric stresses. Such properties include the full elastic
tensor (isothermal or isentropic), the thermal expansion tensor,
and the heat capacity at constant strain.


.. toctree::
  :maxdepth: 3

  anisotropic_materials_01_anisotropic_material
  anisotropic_materials_02_anisotropic_mineral
  anisotropic_materials_03_anisotropic_solution
