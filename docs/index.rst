.. BurnMan documentation master file, created by
   sphinx-quickstart on Sun Aug  1 20:20:43 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BurnMan |version|: a thermodynamic and geophysics toolkit for the Earth and planetary sciences
===============================================================================================


.. image:: burnjack-bar.png

BurnMan is a Python library for computing the thermodynamic and thermoelastic properties
of geological materials from simple mineral endmembers to complex multilayered planetary interiors.

BurnMan is released under the GNU GPL v2 or newer. It relies heavily on numpy, scipy, and matplotlib.

 - Homepage: https://geodynamics.github.io/burnman/
 - Documentation: http://burnman.readthedocs.io
 - Source code: https://github.com/geodynamics/burnman

If you haven't yet installed BurnMan, you can go straight to :ref:`ref-installation` for detailed
instructions. After that, you might want to try out :ref:`ref-example-tutorial` or the other
:ref:`ref-examples`. Finally, and most importantly, have fun!


Citing BurnMan
--------------

If you use BurnMan in your work, we ask that you cite the following publications:

 - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
   Dannberg, J. and Gassmoeller, R. (2023). BurnMan - a Python toolkit for
   planetary geophysics, geochemistry and thermodynamics. Journal of Open Source Software.
   `(https://doi.org/10.21105/joss.05389) <https://doi.org/10.21105/joss.05389>`_

 - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
   Dannberg, J., Gassmoeller, R. and Farla, R. (2024):
   BurnMan v2.1.0 [Software]. Computational Infrastructure for Geodynamics. Zenodo.
   `(https://doi.org/10.5281/zenodo.14238360) <https://doi.org/10.5281/zenodo.14238360>`_

 - Cottaar S., Heister, T., Rose, I., and Unterborn, C. (2014). BurnMan: A
   lower mantle mineral physics toolkit, Geochemistry, Geophysics, and
   Geosystems, 15(4), 1164-1179 `(https://doi.org/10.1002/2013GC005122) <https://doi.org/10.1002/2013GC005122>`_

Contributing to BurnMan
-----------------------

If you would like to contribute bug fixes, new functions or new modules
to the existing codebase, please contact us at info@burnman.org or make a
pull request at `https://github.com/geodynamics/burnman <https://github.com/geodynamics/burnman>`_.

BurnMan also includes a contrib directory that contains python and ipython
scripts used to reproduce published results. We welcome the submission of
new contributions to this directory. As with the contribution of code,
please contact us at info@burnman.org or make a pull request at
`https://github.com/geodynamics/burnman <https://github.com/geodynamics/burnman>`_.

Acknowledgement and Support
---------------------------

 - This project was initiated at, and follow-up research support was received
   through, Cooperative Institute of Deep Earth Research, CIDER (NSF FESD
   grant 1135452) -- see `www.deep-earth.org <http://www.deep-earth.org>`_

 - We thank all the members of the CIDER Mg/Si team for their input:
   Valentina Magni, Yu Huang, JiaChao Liu, Marc Hirschmann, and Barbara
   Romanowicz. We also thank Lars Stixrude for providing benchmarking calculations
   and Zack Geballe, Motohiko Murakami, Bill McDonough, Quentin Williams,
   Wendy Panero, and Wolfgang Bangerth for helpful discussions.

 - We thank CIG (`www.geodynamics.org <http://www.geodynamics.org>`_) for support
   and accepting our donation of BurnMan as an official project.


.. toctree::
  :maxdepth: 4

  introduction
  tutorial
  materials
  anisotropic_materials
  relaxed_materials
  averaging_schemes
  calibrants
  seismic_models
  geotherms
  datasets
  compositions
  layers_and_planets
  material_polytopes
  fitting
  tools
  utils
  examples
  api
  zreferences
  