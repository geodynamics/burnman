# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
Introducing BurnMan |version|
=============================

Overview
--------

BurnMan is an open source mineral physics and seismological toolkit written in Python
which can enable a user to calculate (or fit) the physical and chemical properties
of endmember minerals, fluids/melts, solutions, and composite assemblages.

Properties which BurnMan can calculate include:

  - the thermodynamic free energies, allowing phase equilibrium calculations,
    endmember activities, chemical potentials and oxygen (and other) fugacities.

  - entropy, enabling the user to calculate isentropes for a given assemblage.

  - volume, to allow the user to create density profiles.

  - seismic velocities, including Voigt-Reuss-Hill and Hashin-Strikman bounds
    and averages.

The toolkit itself comes with a large set of classes and functions which are
designed to allow the user to easily combine mineral physics with
geophysics, and geodynamics. The features of BurnMan include:

  - the full codebase, which includes implementations of many static and thermal equations of state
    (including Vinet, Birch Murnaghan, Mie-Debye-Grueneisen, Modified Tait),
    and solution models (ideal, symmetric, asymmetric, subregular).
  - popular endmember and solution datasets already coded into burnman-usable format
    (including :cite:`HP2011`, :cite:`Stixrude2005`, :cite:`Stixrude2011` and :cite:`Stixrude2022`)
  - Optimal least squares fitting routines for multivariate data with (potentially correlated) errors
    in pressure and temperature. As an example, such functions can be used to
    simultaneously fit volumes, seismic velocities and enthalpies.
  - a "Planet" class, which self-consistently calculates gravity profiles, mass, moment of
    inertia of planets given the chemical and temperature structure of a planet
  - published geotherms
  - a tutorial on the basic use of BurnMan
  - a large collection of annotated examples
  - a set of high-level functions which create files readable by seismological and geodynamic software,
    including: Mineos :cite:`Masters2011`, AxiSEM :cite:`NissenMeyer2014` and ASPECT
  - an extensive suite of unit tests to ensure code functions as intended
  - a series of benchmarks comparing BurnMan output with published data
  - a directory containing user-contributed code from published papers

BurnMan makes extensive use of `SciPy <http://www.scipy.org/>`_,
`NumPy <http://www.numpy.org/>`_ and `SymPy <http://www.sympy.org/>`_
which are widely used Python libraries for scientific computation.
`Matplotlib <http://matplotlib.org/>`_ is used
to display results and produce publication quality figures.
The computations are consistently formulated in terms of SI units.

The code documentation including class and function descriptions can be found online at
http://burnman.readthedocs.io.

This software has been designed to allow the end-user a great deal of freedom
to do whatever calculations they may wish and to add their own modules.
The underlying Python classes have been designed to make new endmember,
solution and composite models easy to read and create.
We have endeavoured to provide examples and benchmarks which cover the
most popular uses of the software, some of which are included in the figure below.
This list is certainly not exhaustive, and we will definitely have missed interesting
applications. We will be very happy to accept contributions in
form of corrections, examples, or new features.

Structure
---------
.. image:: figures/structure.png


.. _ref-installation:


Requirements
------------

  - Python 3.8+
  - Python modules: NumPy, SciPy, SymPy, Sparse, Matplotlib

Optional modules
^^^^^^^^^^^^^^^^

Needed for some functionality:

  - cvxpy: required for some least squares fitting routines
    and solution polytope calculations.
  - pycddlib: required for solution polytope calculations.
  - autograd: required for esoteric solution models defined using a single excess
    function. Not required for the vast majority of users.

Installation
------------

Installation of BurnMan is mostly platform independent.
As long as you know how to use a terminal, the process should be straightforward.
The following instructions should help, but let us know if you have any problems.

Dependencies
^^^^^^^^^^^^
First, make sure you have a sufficiently recent version of python installed
on your machine (see above for the latest requirements).
To check your version of python, type the following in a terminal:
    python --version
If your version is not recent enough, visit https://www.python.org/ to
find out how to install a newer version.

Once you have checked your version of python, you should make sure you have
installed the python module pip. We will use this module to install BurnMan.
If you don't have it already, you can install it by opening a
terminal window and typing:

    python -m ensurepip --upgrade

Mac users will also need to install Xcode, which can be found in the MacStore.

Stable version
^^^^^^^^^^^^^^

If you are only interested in using BurnMan
(rather than developing the software), and you aren't interested in any of the
latest changes, you can install the stable version by typing
the following into a terminal window:

    python -m pip install burnman

This method of installation does not give you easy access to all the examples,
or the test suite. These can be found in the latest release package which can
be downloaded from https://github.com/geodynamics/burnman/releases.

Development version
^^^^^^^^^^^^^^^^^^^
If you want to install the development version of BurnMan
(with all the latest features), you will first need to download the source code.
The best way to do this is by using git (a version control system).
To install git, follow the instructions at https://git-scm.com/downloads.

Then, using a terminal, navigate to the directory into which you want to
clone the BurnMan repository, and type

    git clone https://github.com/geodynamics/burnman.git

(If you don't want to use git, you can download the current master branch
from https://github.com/geodynamics/burnman/archive/master.zip.)

Once the repository is cloned, navigate to the top-level directory by typing
`cd burnman` in the terminal, and then install BurnMan, either in static mode:
`python -m pip install .` or in development mode
(if you want to develop or modify the code): `python -m pip install -e .`.

Checking that the installation worked
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To check that the installation has worked, you can run the test suite
(`./test.sh`). This takes a few minutes to run.

A more basic check that BurnMan is installed is to navigate to the Burnman
examples directory and type:

    python example_beginner.py

If figures show up, BurnMan has been installed.


Citing BurnMan
--------------

If you use BurnMan in your work, we ask that you cite the following publications:

  - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
    Dannberg, J. and Gassmoeller, R. (2023). BurnMan - a Python toolkit for
    planetary geophysics, geochemistry and thermodynamics. Journal of Open Source Software.
    `(https://doi.org/10.21105/joss.05389) <https://doi.org/10.21105/joss.05389>`_

  - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
    Dannberg, J., Gassmoeller, R. and Farla, R. (2024):
    BurnMan v2.0.0 [Software]. Computational Infrastructure for Geodynamics. Zenodo.
    `(https://doi.org/10.5281/zenodo.14165625) <https://doi.org/10.5281/zenodo.14165625>`_

  - Cottaar S., Heister, T., Rose, I., and Unterborn, C., (2014). BurnMan: A
    lower mantle mineral physics toolkit, Geochemistry, Geophysics, and
    Geosystems, 15(4), 1164-1179 `(https://doi.org/10.1002/2013GC005122)
    <https://doi.org/10.1002/2013GC005122>`_

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

"""
from __future__ import absolute_import
import importlib.metadata

# Low level utility functions
from . import utils
from .utils import geotherm

# Classes and associated functions for representing rocks and minerals:
from .classes.material import Material, material_property
from .classes.perplex import PerplexMaterial
from .classes.mineral import Mineral
from .classes.combinedmineral import CombinedMineral
from .classes.solution import Solution, SolidSolution, RelaxedSolution
from .classes.elasticsolutionmodel import ElasticSolutionModel
from .classes.elasticsolution import ElasticSolution, ElasticSolidSolution
from .classes.composite import Composite
from .classes.calibrant import Calibrant
from .classes.anisotropy import AnisotropicMaterial
from .classes.anisotropicmineral import AnisotropicMineral
from .classes.anisotropicmineral import cell_parameters_to_vectors
from .classes.anisotropicmineral import cell_vectors_to_parameters
from .classes.anisotropicsolution import AnisotropicSolution, RelaxedAnisotropicSolution
from .classes.mineral_helpers import HelperLowHighPressureRockTransition
from .classes.mineral_helpers import HelperSpinTransition
from .classes.mineral_helpers import HelperRockSwitcher

# Other classes
from .classes.composition import Composition
from .classes.layer import Layer, BoundaryLayerPerturbation
from .classes.planet import Planet
from .classes.polytope import MaterialPolytope
from .classes import seismic
from .classes import averaging_schemes

# Mineral library
from . import minerals

# Equations of state
from . import eos

# Calibrants
from . import calibrants

# High level tools
from . import tools
from .tools.equilibration import equilibrate
from .tools.partitioning import calculate_nakajima_fp_pv_partition_coefficient

# Optimization functions
from .optimize import composition_fitting
from .optimize import linear_fitting
from .optimize import nonlinear_fitting
from .optimize import nonlinear_solvers
from .optimize import eos_fitting

__version__ = importlib.metadata.version("burnman")
