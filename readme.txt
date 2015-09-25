    This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences.
    Copyright (C) 2012 - 2015 by the BurnMan team
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Contact: Ian Rose 
ian.rose@berkeley.edu

*** About

BurnMan is a lower mantle shear velocity generator constrained by mineral physics.

Work started by a work-group of the CIDER 2012 workshop in Santa Barbara.

BurnMan is released under the GNU GPL v2 or newer

Homepage: http://burnman.org
Documentation: http://burnman.org/current-doc
Source code: https://github.com/geodynamics/burnman

Authors (as of 2015, listed alphabetically by first name):
	Bob Myhill
	Cayman Unterborn
	Ian Rose
	Sanne Cottaar
	Timo Heister

*** Requirements

- Python 2.7 (not Python 3.x)
- Python modules:
  NumPy, SciPy, matplotlib


*** Install under Ubuntu

1. sudo apt-get install python python-scipy python-numpy python-matplotlib 
2. run with "python main.py" in a shell


*** Install on a MAC

0. get Xcode
1. If you don't have Python yet, download it (for free) from
   python.org/download . Make sure to use the latest version 2.x version (I
   used 2.7). To check your version of python, type the following in a
   terminal: 
     python --version
2. Install the latest Numpy version: http://sourceforge.net/projects/numpy/files/NumPy/
3. Install the latest Scipy at http://sourceforge.net/projects/scipy/files/
4. Install the latest Matplotlib from http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.1.1/
5. Go to the main BurnMan directory and type:
	python example_composition.py
    Figures should show up.

Problems you might run into:

  - Installing numpy/scipy/matplotlib for a different python version than the
    one on your computer

  - Having matplotlib for 32-bit instead of 64-bit (for me this got fixed by
    installing the very latest version). This will give you the error 'no
    matching architecture in universal wrapper'. You can check if your python
    is 32 or 64 bit with the following lines:
	python 
	>>> import platform
	>>> print platform.architecture()


*** Install under Windows

make Python 2.7.3 (for example) running under windows (do not use Python 3.x, but 2.7.x):

1. Download Python from http://www.python.org/ and install the version at C:\Python27\; the 32-bit version is recommended
2. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy, download "numpy-MKL-1.6.2.win32-py2.7.exe" and install
3. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy, download "scipy-0.10.1.win32-py2.7.exe" and install
4. Go to http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib, download "matplotlib-1.1.1.win32-py2.7.exe" and install
5. Open Python Shell (IDLE Python GUI)
6. File -- Open -- find one of the example files
7. Run the module (or press F5)

*** Start Here

To begin, the user may want to look at these examples to begin to understand
what tools are available in BurnMan and how values are calculated. Below is a
suggested order of examples that begin by introducing each of the user inputs
possible as well as each of the helpers involved with each example.

1. example_beginner.py (Creating a composite and computing seismic properties)
1. example_geotherms.py (Explains each built in geothermal and how to create
   your own) run this example by typing 'python example_geotherms.py' this
   should result in a plot showing various geotherms available
2. example_seismic.py (Explains the various seismic models included in
   BurnMan)
3. example_composition.py (Explains each of the possible input mineral
   creation schemes) by changing the True/False flags, one can plot the
   different compositions defined in this script
4. example_user_input_materials.py (Explains how to create user-defined
   minerals)
5. example_averaging.py (Explains how moduli and density are averaged to
   calculate seismic velocities)


*** About scripting in Python

Burnman has the advantage of being adaptable and extendable in easy scripts. The downside might be that we do not
provide a graphical user interface. For those of you who are not familiar  with python, we suspect it will still be 
relatively easy to adapt the scripts for computations and plotting. 
Here are some specific features and pitfalls on python:
    - Python uses specific identation. A script might fail if a code block is not idented correctly. We use four spaces and no tabs, 
      mixing these can give trouble.
    - Indices require square brackets and function or method calls parentheses (mainly different from Matlab).
    - The first index of an array is 0 (e.g. x[0])
    - Put dots after numbers to make them floats instead of integers (e.g. 5/3 will give 1 (rounded downward), while 5./3. will give 1.66666666667)


*** Examples 
    (Note: the following text is automatically generated from misc/gen_doc.py,
    do not edit)

* example_beginner.py

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
* :func:`burnman.main.velocities_from_rock`


*Demonstrates:*

* creating basic composites
* calculating thermoelastic properties
* seismic comparison


example_geotherms
-----------------

This example shows each of the geotherms currently possible with BurnMan.
These are:

1. Brown and Shankland, 1981
2. Anderson, 1982
3. Watson and Baxter, 2007
4. linear extrapolation
5. Read in from file from user
6. Adiabatic from potential temperature and choice of mineral (pyrolite in this example)

*Uses:*

* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.geotherm.anderson`
* input geotherm file *input_geotherm/example_geotherm.txt* (optional)
* :class:`burnman.composite.Composite` for adiabat

*Demonstrates:*

* the available geotherms



example_seismic
---------------

Shows the various ways to input seismic models (Vs, Vp, Vphi, Density) as a
function of depth (or P) as well as different velocity model libraries
available within Burnman:

1. PREM (Dziewonski & Anderson, 1981)
2. Reference model for fast regions (outside the LLSVP's) in the lower mantle
   (Lekic et al. 2012)
3. Reference model for slow regions (LLSVP's) in the lower mantle (Lekic et la. 2012)

This example will first calculate or read in a seismic model and plot the
model along the defined pressure range.

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
2. Two minerals mixed in simple mole fractions with user-defined Fe
   partitioning
3. The user can input wt% of each cation (Mg, Fe and Si) and BurnMan will
   calculate Fe partitioning along a P, T profile (see
   example_partition_coef.py)
4. A mixture of three minerals.

In compositions 2, 3, and 4 of the above inputs, BurnMan will mix the mineral
physical paremeters of end member minerals (pure Mg and Fe) of the user's
choice using either volumetric (moduli) or molar averaging (all others) at
room pressure and temperature (see example_user_input_material.py for
information on these parameters).

To turn a method of mineral creation "on" the first if statement above the
method must be set to True, with all others set to False.

Note: These minerals can include a spin transition in (Mg,Fe)O, see
example_spintransition.py for explanation of how to implement this

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`


*Demonstrates:*

* Different ways to define a composite
* Compare computations to seismic models



    
example_user_input_material
---------------------------

The main focus of this example is to show the mineral physical input constants
necessary for BurnMan to calculate seismic velocity profiles. Furht

Shows user how to input a mineral of his/her choice and which physical values
need to be input for BurnMan to calculate Vs, Vp, Vphi and density at depth.

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

See Watt et al., 1976 Journal of Geophysics and Space Physics for explanations
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
a plot of Vs, Vphi and density is produce for the user to compare each of the
different methods.

*Specifically uses:*

* :doc:`eos`


*Demonstrates:*

* Each method for calculating velocity profiles currently included within BurnMan



    
example_spintransition
----------------------

This example shows the different minerals that are implemented with a spin
transition.  Minerals with spin transition are implemented by defining two
separate minerals (one for the low and one for the high spin state).  Then a
third dynamic mineral is created that switches between the two previously
defined minerals by comparing the current pressure to the transition pressure.

*Specifically uses:*

* :func:`burnman.mineral_helpers.HelperSpinTransition`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase_HS`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase_LS`


*Demonstrates:*

* implementation of spin transition in (Mg,Fe)O at user defined pressure



example_parition_coef
---------------------

This example shows how to vary the distribution coefficient of the
perovskite/ferropericlase system. The user sets Kd_0 and BurnMan scales Kd as
a function of P and T adopting the formalism of Nakajima et al.,
2012. Specifically we adopt equation 5 of Nakajima et al., 2012 with DeltaV_0
= 0.2 cc/mol, and calculating the partition coefficient of Fe in each phase
from stoichiometry.

This example will calculate mineral input parameters from Mg and Fe endmembers
from Stixrude and Lithgow-bertelloni, 2005 with weighting determined by the
calculated partition coefficients. Finally, output plots of X_Fe pv and X_Fe
fp our output as well as the user's choice of geotherm

requires:
- geotherms
-input distribution coefficient Kd_0

teaches:
- creating varying proportions of Fe and its effect on seismic velocities




example_compare_enstpyro
------------------------

This example shows you how to create two materials from wt% determines the
optimum mixing between the two to match the seismic model of your choice.
Currently it compares two end member meteorite groups among the chondrites:
carbonaceous and enstatite. Velocities are calculated for each set of minerals
and plotted for comparison.

requires:
- geotherms
- seismic models
- compute seismic velocities
- creating minerals

teaches:
- weight percent materials



    
example_optimize_pv
-------------------

Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM.

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.composite.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.main.velocities_from_rock`
* :func:`burnman.main.compare_l2`

*Demonstrates:*

* compare errors between models
* loops over models


