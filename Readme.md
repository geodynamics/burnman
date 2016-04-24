# BurnMan - a thermoelastic and thermodynamic toolkit for Earth and planetary sciences
<img src="sphinx/burnjack.png", width="256">
## About

BurnMan is a Python library for generating thermodynamic and thermoelastic models of planetary interiors.

It began as a working-group at the 2012 CIDER workshop in Santa Barbara.

BurnMan is released under the GNU GPL v2 or newer

Homepage: http://burnman.org

Documentation: http://burnman.org/doc-0.9.0

Source code: https://github.com/geodynamics/burnman

Authors (as of 2015, listed alphabetically by first name):
* Bob Myhill
* Cayman Unterborn
* Ian Rose
* Sanne Cottaar
* Timo Heister

Contact: Ian Rose 
ian.rose@berkeley.edu

## Requirements

* Python 2.7.x or Python 3.4+
* Python modules:
  NumPy, SciPy, Matplotlib


## Install under Ubuntu

1. Install using apt by opening a terminal window and entering
`sudo apt-get install python python-scipy python-numpy python-matplotlib` 
2. Go to the Burnman examples directory and type:
```python example_beginner.py```
Figures should show up, indicating that it is working.


## Install on a Mac

1. get Xcode
2. If you don't have Python yet, download it (for free) from
   python.org/download . Make sure to use either Python 2.7 or Python 3.4+.
   To check your version of python, type the following in a
   terminal: 
     python --version
3. Install the latest Numpy version: http://sourceforge.net/projects/numpy/files/NumPy/
4. Install the latest Scipy at http://sourceforge.net/projects/scipy/files/
5. Install the latest Matplotlib from http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.1.1/
6. Go to the Burnman examples directory and type:
	python example_beginner.py
    Figures should show up, indicating that it is working.

Problems you might run into:

* Installing numpy/scipy/matplotlib for a different python version than the
  one on your computer

* Having matplotlib for 32-bit instead of 64-bit (for me this got fixed by
  installing the very latest version). This will give you the error `no
  matching architecture in universal wrapper`. You can check if your python
  distribution is 32 or 64 bit with the following lines:
```
python 
>>> import platform
>>> print platform.architecture()
```

## Install under Windows

To get Python 2.7.x (for example) running under Windows:

1. Download Python from http://www.python.org/ and install the version at C:\Python27\; the 32-bit version is recommended
2. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy, download "numpy-MKL-1.6.2.win32-py2.7.exe" and install
3. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy, download "scipy-0.10.1.win32-py2.7.exe" and install
4. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib, download "matplotlib-1.1.1.win32-py2.7.exe" and install
5. Open Python Shell (IDLE Python GUI)
6. File -- Open -- find one of the example files
7. Run the module (or press F5)

## Start Here

To begin, the user may want to look at these examples to begin to understand
what tools are available in BurnMan and how values are calculated. Below is a
suggested order of examples that begin by introducing each of the user inputs
possible as well as each of the helpers involved with each example.

1. example_beginner.py (Creating a composite and computing seismic properties)
1. example_geotherms.py (Demonstrates built in geotherms and how to create
   your own).
2. example_seismic.py (Explains the various seismic models included in
   BurnMan)
3. example_composition.py (Explains how to create different mineralogical models) 
4. example_user_input_materials.py (Explains how to create user-defined
   minerals)
5. example_averaging.py (Explains how moduli and density are averaged to
   calculate seismic velocities)


## About scripting in Python

Burnman has the advantage of being adaptable and extendable in easy scripts. The downside might be that we do not
provide a graphical user interface. For those of you who are not familiar  with python, we suspect it will still be 
relatively easy to adapt the scripts for computations and plotting. 
Here are some specific features and pitfalls on Python:

* Python uses specific identation. A script might fail if a code block is not idented correctly. We use four spaces and no tabs, 
  mixing these can give trouble.
* Indices require square brackets and function or method calls parentheses (mainly different from Matlab).
* The first index of an array is 0 (e.g. x[0])
* Put dots after numbers to make them floats instead of integers (e.g. 5/3 will give 1 (Python 2.x rounds downward), while 5./3. will give 1.66666666667)


## Examples 

example_beginner
----------------

This example script is intended for absolute beginners to BurnMan.
We cover importing BurnMan modules, creating a composite material,
and calculating its seismic properties at lower mantle pressures and
temperatures.  Afterwards, we plot it against a 1D seismic model
for visual comparison.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.material.Material.evaluate`


*Demonstrates:*

* creating basic composites
* calculating thermoelastic properties
* seismic comparison


example_geotherms
-----------------

This example shows each of the geotherms currently possible with BurnMan.
These are:

1. Brown and Shankland, 1981 :cite:`Brown1981`
2. Anderson, 1982 :cite:`anderson1982earth`
3. Watson and Baxter, 2007 :cite:`Watson2007`
4. linear extrapolation
5. Read in from file from user
6. Adiabatic from potential temperature and choice of mineral

*Uses:*

* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.geotherm.anderson`
* input geotherm file *input_geotherm/example_geotherm.txt* (optional)
* :class:`burnman.composite.Composite` for adiabat

*Demonstrates:*

* the available geotherms



example_seismic
---------------

Shows the various ways to input seismic models (:math:`V_s, V_p, V_{\phi}, \rho`) as a
function of depth (or pressure) as well as different velocity model libraries
available within Burnman:

1. PREM :cite:`dziewonski1981`
2. STW105 :cite:`kustowski2008`
3. AK135 :cite:`kennett1995`
4. IASP91 :cite:`kennett1991`

This example will first calculate or read in a seismic model and plot the
model along the defined pressure range. The example also illustrates how to import a seismic model of your choice, here shown by importing AK135 :cite:`kennett1995`.

*Uses:*

* :doc:`seismic`



*Demonstrates:*

* Utilization of library seismic models within BurnMan
* Input of user-defined seismic models





example_composition
-------------------

This example shows how to create different minerals, how to compute seismic
velocities, and how to compare them to a seismic reference model.

There are many different ways in BurnMan to combine minerals into a
composition. Here we present a couple of examples:

1. Two minerals mixed in simple mole fractions. Can be chosen from the BurnMan
   libraries or from user defined minerals (see example_user_input_material)
2. Example with three minerals
3. Using preset solid solutions
4. Defining your own solid solution


To turn a method of mineral creation "on" the first if statement above the
method must be set to True, with all others set to False.

Note: These minerals can include a spin transition in (Mg,Fe)O, see
example_spintransition.py for explanation of how to implement this

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.mineral.Mineral`
* :class:`burnman.solidsolution.SolidSolution`

*Demonstrates:*

* Different ways to define a composite
* Using minerals and solid solutions
* Compare computations to seismic models




example_user_input_material
---------------------------

Shows user how to input a mineral of his/her choice without usint the library and which physical values
need to be input for BurnMan to calculate :math:`V_P, V_\Phi, V_S` and density at depth.

*Specifically uses:*


* :class:`burnman.mineral.Mineral`

*Demonstrates:*

* how to create your own minerals



example_averaging
-----------------

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:

1. Voight-Reuss-Hill
2. Voight averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See :cite:`Watt1976` Journal of Geophysics and Space Physics for explanations
of each averaging scheme.

*Specifically uses:*

* :class:`burnman.averaging_schemes.VoigtReussHill`
* :class:`burnman.averaging_schemes.Voigt`
* :class:`burnman.averaging_schemes.Reuss`
* :class:`burnman.averaging_schemes.HashinShtrikmanUpper`
* :class:`burnman.averaging_schemes.HashinShtrikmanLower`

*Demonstrates:*

* implemented averaging schemes



example_woutput
---------------

This example explains how to perform the basic i/o of BurnMan. A method of
calculation is chosen, a composite mineral/material (see
example_composition.py for explanation of this process) is created in the
class "rock," finally a geotherm is created and seismic velocities calculated.

Post-calculation, the results are written to a simple text file to
plot/manipulate at the user's whim.

requires:
- creating minerals
- compute seismic velocities
- geotherms

teaches:
- output computed seismic data to file




example_compare_all_methods
---------------------------

This example demonstrates how to call each of the individual calculation
methodologies that exist within BurnMan. See below for current options. This
example calculates seismic velocity profiles for the same set of minerals and
a plot of :math:`V_s, V_\phi` and :math:`\rho` is produce for the user to compare each of the
different methods.

*Specifically uses:*

* :doc:`eos`


*Demonstrates:*

* Each method for calculating velocity profiles currently included within BurnMan




example_optimize_pv
-------------------

Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM. For more extensive comments on this setup, see tutorial/step_2.py

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.material.Material.evaluate`
* :func:`burnman.main.compare_l2`

*Demonstrates:*

* compare errors between models
* loops over models




example_fit_data
----------------

This example demonstrates BurnMan's functionality to fit thermoelastic data to
both 2nd and 3rd orders using the EoS of the user's choice at 300 K. User's
must create a file with :math:`P, T` and :math:`V_s`. See input_minphys/ for example input
files.

requires:
- compute seismic velocities

teaches:
- averaging



example_grid
------------

This example shows how to evaluate seismic quantities on a :math:`P,T` grid.



example_chemical_potentials
---------------------------

This example shows how to use the chemical potentials library of functions.

*Demonstrates:*

* How to calculate chemical potentials
* How to compute fugacities and relative fugacities




example_solid_solution
----------------------

This example shows how to create different solid solution models and output
thermodynamic and thermoelastic quantities.

There are four main types of solid solution currently implemented in
BurnMan:

1. Ideal solid solutions
2. Symmmetric solid solutions
3. Asymmetric solid solutions
4. Subregular solid solutions

These solid solutions can potentially deal with:

* Disordered endmembers (more than one element on a crystallographic site)
* Site vacancies
* More than one valence/spin state of the same element on a site

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.solidsolution.SolidSolution`
* :class:`burnman.solutionmodel.SolutionModel`


*Demonstrates:*

* Different ways to define a solid solution
* How to set composition and state
* How to output thermodynamic and thermoelastic properties



example_build_planet
--------------------

For Earth we have well-constrained one-dimensional density models.  This allows us to
calculate pressure as a funcion of depth.  Furthermore, petrologic data and assumptions
regarding the convective state of the planet allow us to estimate the temperature.

For planets other than Earth we have much less information, and in particular we
know almost nothing about the pressure and temperature in the interior.  Instead, we tend
to have measurements of things like mass, radius, and moment-of-inertia.  We would like
to be able to make a model of the planet's interior that is consistent with those
measurements.

However, there is a difficulty with this.  In order to know the density of the planetary
material, we need to know the pressure and temperature.  In order to know the pressure,
we need to know the gravity profile.  And in order to the the gravity profile, we need
to know the density.  This is a nonlinear problem which requires us to iterate to find
a self-consistent solution.

Here we show an example that does this, using the planet Mercury as motivation.


*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :func:`burnman.material.Material.evaluate`

