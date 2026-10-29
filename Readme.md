[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14238360.svg)](https://doi.org/10.5281/zenodo.14238360)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.05389/status.svg)](https://doi.org/10.21105/joss.05389)
[![coverage](https://img.shields.io/badge/coverage-94%25-green")](https://github.com/geodynamics/burnman/actions/workflows/coverage.yml)

# BurnMan - a Python toolkit for planetary geophysics, geochemistry and thermodynamics

<img src="docs/burnjack.png" width="256">

## About

BurnMan is a free, open-source toolkit written in Python. It is designed to compute thermodynamic and thermoelastic properties of geological and planetary materials — from pure mineral ``endmembers'', through solid solutions, to composites of multiple phases, through to planetary layers and whole planets.

BurnMan is a modular, class-based toolkit, and as such it has a wide range of use cases. It is well-suited to both discipline-specific and interdisciplinary problems, including:
- Fitting experimental data on mineral/melt endmembers (the `Mineral` class) and solutions (the `Solution` class) to a range of equations of state and solution models
- Designing experimental starting mix compositions (via the `Composition` class) and performing mass balance calculations on post-run products (via the `burnman.optimize.composition_fitting` submodule)
- Predicting thermodynamic properties from user-generated or existing thermodynamic datasets (e.g. Stixrude + Lithgow-Bertelloni, Holland, Powell et al.; see the `minerals` submodule)
- Creation and property interrogation of fully anisotropic models of endmembers (the `AnisotropicMineral` class) and solutions (the `AnistropicSolution` class)
- Pressure-temperature estimation of experiments from XRD data (via the `Calibrant` class)
- Pressure-temperature estimation of natural samples (via the `optimize.thermobarometry.estimate_conditions` function)
- Calculating equilibrium compositions of known mineralogical assemblages given constraints (via the `equilibrate` function)
- Calculating assemblage boundaries or isopleths, including uncertainties (also via the `equilibrate` function)
- Calculating seismic velocities of multiphase assemblages according to different averaging schemes (e.g. Reuss, Voigt, VRH, Hashin-Strikman) and with or without isochemical reaction during passage of the wave (``seismic softening'', see the `RelaxedSolution` and `RelaxedAnisotropicSolution` class).
- Generating 1D profiles of planetary interiors (via the `Layer` and `Planet` classes).

A variety of tutorials, examples and tests accompany the BurnMan module, and can be used as a basis for your own work.


BurnMan began as a working-group at the 2012 CIDER workshop in Santa Barbara. It is released under the GNU GPL v2 or newer

* Homepage: https://geodynamics.github.io/burnman/
* Documentation: http://burnman.readthedocs.io
* Source code: https://github.com/geodynamics/burnman
* Forums: https://community.geodynamics.org/c/burnman

Authors (as of 2025):
* Bob (Robert) Myhill (main contributor)
* Sanne Cottaar
* Timo Heister
* Ian Rose
* Cayman Unterborn
* Benedict Heinen
* Robert Farla
* Juliane Dannberg
* Rene Gassmoeller

## Citing BurnMan

If you use BurnMan in your work, we ask that you cite the following publications:

  - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
    Dannberg, J. and Gassmoeller, R. (2023). BurnMan - a Python toolkit for
    planetary geophysics, geochemistry and thermodynamics. Journal of Open Source Software.
    https://doi.org/10.21105/joss.05389

  - Myhill, R., Cottaar, S., Heister, T., Rose, I., Unterborn, C.,
    Dannberg, J., Gassmoeller, R. and Farla, R. (2024):
    BurnMan v2.1.0 [Software]. Computational Infrastructure for Geodynamics. Zenodo.
    https://doi.org/10.5281/zenodo.14238360

  - Cottaar S., Heister, T., Rose, I., and Unterborn, C., (2014). BurnMan: A
    lower mantle mineral physics toolkit, Geochemistry, Geophysics, and
    Geosystems, 15(4), 1164-1179 https://doi.org/10.1002/2013GC005122
    
## Installation requirements

* Python 3.10+
* Python modules:
  NumPy, SciPy, SymPy, Matplotlib

## Optional packages needed for some functionality

* cvxpy - required for some least squares fitting routines and solution polytope calculations.
* pycddlib - required for solution polytope calculations.
* autograd - required for esoteric solution models defined using a single excess function. Not required for the vast majority of users.

## Installation process
Installation of BurnMan is mostly platform independent.
As long as you know how to use a terminal, the process should be straightforward.
The following instructions should help, but let us know if you have any problems.

### Environment management
We strongly recommend using a python environment manager like conda or pyenv to install
BurnMan and its dependencies. This is especially the case for installations on modern Mac systems.

For pyenv, we suggest you select the most recent version of python supported by BurnMan, install all of the dependencies into that environment and set the burnman root directory to use that environment automatically.

For conda, we suggest making a new environment for the burnman installation, install the most recent version of python supported by BurnMan into that environment, and install all of the dependencies into that environment. Remember to activate the environment before installing new dependencies or using BurnMan.

### Dependencies
First, make sure you have a sufficiently recent version of python installed on your machine (see above for the latest requirements).  To check your version of python, type the following in a terminal:
    python --version
If your version is not recent enough, visit https://www.python.org/ to find out how to install a newer version.

Once you have checked your version of python, you should make sure you have installed the python module pip. We will use this module to install BurnMan. If you don't have it already, you can install it by opening a terminal window and typing:

    python -m ensurepip --upgrade

Mac users will also need to install Xcode, which can be found in the MacStore, or can be installed with:

    xcode-select --install

### Stable version
If you are only interested in using BurnMan (rather than developing the software), and you aren't interested in any of the latest changes, you can install the stable version by typing the following into a terminal window:

    python -m pip install burnman

This method of installation does not give you easy access to all the examples, or the test suite. These can be found in the latest release package which can be downloaded from https://github.com/geodynamics/burnman/releases.

### Development version
If you want to install the development version of BurnMan (with all the latest features), you will first need to download the source code. The best way to do this is by using git (a version control system). To install git, follow the instructions at https://git-scm.com/downloads.

Then, using a terminal, navigate to the directory into which you want to clone the BurnMan repository, and type

    git clone https://github.com/geodynamics/burnman.git

(If you don't want to use git, you can download the current master branch from https://github.com/geodynamics/burnman/archive/master.zip.)

Once the repository is cloned, navigate to the top-level directory by typing `cd burnman` in the terminal, and then install BurnMan, either in static mode: `python -m pip install .` or in development mode (if you want to develop or modify the code): `python -m pip install -e .`.

### Checking that the installation worked

To check that the installation has worked, you can run the test suite (`./test.sh`). This takes a few minutes to run. You may find that you need to install some other dependencies (latex, pycddlib) if you don't already have them on your system.

A more basic check that BurnMan is installed is to navigate to the Burnman examples directory and type:

    python example_beginner.py

If figures show up, BurnMan has been installed.


## Getting started: The BurnMan Tutorial

An introduction to the tools available in BurnMan can be found in the BurnMan
Tutorial. These are ipython notebooks that can be run from inside jupyterlab.

If you want to run these notebooks without installing BurnMan, you can access
them on binder.org via the links below.

1. [tutorial_01_material_classes.ipynb](https://mybinder.org/v2/gh/geodynamics/burnman/HEAD?labpath=tutorial%2Ftutorial_01_material_classes.ipynb "Part 1: The Material Classes") (introduces the main classes in BurnMan, which are used to calculate material properties)
2. [tutorial_02_composition_class.ipynb](https://mybinder.org/v2/gh/geodynamics/burnman/HEAD?labpath=tutorial%2Ftutorial_02_composition_class.ipynb "Part 2: The Composition Class") (introduces the composition class, used for processing and converting chemical compositions)
3. [tutorial_03_layers_and_planets.ipynb](https://mybinder.org/v2/gh/geodynamics/burnman/HEAD?labpath=tutorial%2Ftutorial_03_layers_and_planets.ipynb "Part 3: Layers and Planets") (introduces the layer and planet classes used to create models of planetary interiors)
4. [tutorial_04_fitting.ipynb](https://mybinder.org/v2/gh/geodynamics/burnman/HEAD?labpath=tutorial%2Ftutorial_04_fitting.ipynb "Part 4: Fitting") (introduces the various functions used to fit model parameters to experimental and analytical data)
5. [tutorial_05_equilibrium.ipynb](https://mybinder.org/v2/gh/geodynamics/burnman/HEAD?labpath=tutorial%2Ftutorial_05_equilibrium.ipynb "Part 5: Equilibrium problems") (introduces the equilibrate function, used to equilibrate mineral assemblages)


## More detail: The Examples Suite

The BurnMan tutorials provide a basic but incomplete understanding of what
the module can do. To supplement the tutorials, BurnMan includes a
large suite of examples that provide more in-depth coverage of the
potential uses of the module.

For an up-to-date summary of all the examples, including the generated figures,
the user is referred to the BurnMan manual (http://burnman.readthedocs.io).

## About scripting in Python

Burnman has the advantage of being adaptable and extensible in easy scripts.
As BurnMan is a toolkit, a graphical user interface would be impractical
(although we do plan to add GUIs for some functionality in the future).
Nevertheless, we hope that we have succeeded in making BurnMan accessible to
coding novices. For those of you who have little experience with Python,
here are some specific features and pitfalls of the language:

* Python uses specific indentation. A script might fail if a code block is not indented correctly. We use four spaces and no tabs. Mixing spaces and tabs can cause trouble.
* Indices should be given inside square brackets and function or method call arguments inside parentheses (different from Matlab).
* The first index of an array or list is 0 (e.g. x[0]), not 1.
* Put dots after numbers to make them floats instead of integers.
