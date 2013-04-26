    BurnMan - a lower mantle toolkit
    Copyright (C) 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.

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

see: http://code.google.com/p/BurnMan/

BurnMan is released under the GNU GPL v2 or newer

see burnman.org for more information


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

*** Examples

* example_geotherms.py

Shows the various ways to input geotherms: Built-in geotherms (geotherm1 and 2), basic linear (geotherm3),
loaded in from a data file (geotherm4) of your choice. Geotherm 1 is from Brown & Shankland (1981) and 
geotherm2 from Watson & Baxter (2007).

requires:

teaches:
- geotherms


* example_geotherms.py

Shows the various ways to input geotherms: Built-in geotherms (geotherm1 and 2), basic linear (geotherm3),
loaded in from a data file (geotherm4) of your choice. Geotherm 1 is from Brown & Shankland (1981) and 
geotherm2 from Watson & Baxter (2007).

requires:

teaches:
- geotherms


* example_seismic.py

Shows the various ways to input seismic models (Vs, Vp, Vphi, Density) as a
function of depth (or P) as well as different velocity models available:
PREM (Dziewonski & Anderson, 1981)
reference model for fast regionsi (outside the LLSVP's) in the lower mantle (Lekic et al. 2012)
reference model for slow regions (LLSVP's) in the lower mantle (Lekic et la. 2012)

requires:

teaches:
- seismic models


* example_compare_two_models.py

Calculates and plots two models for different minerals or methods and plots
the results. Calculates basic percent difference statistics as well.

requires:
- geotherms
- creating minerals

teaches:
- compute seismic velocities and compare

* example_composition.py

This example shows how to create different minerals, how to compute seismic
velocities, and how to compare them to a seismic reference model.

requires:
- geotherms
- seismic models
- compute seismic velocities

teaches:
- creating minerals
- seismic comparison


* example_woutput.py

Compute properties of minerals and creates a table of outputs in a text format
that could be used with other programs.

requires:
- creating minerals
- compute seismic velocities
- geotherms

teaches:
- output computed seismic data to file 


* example_user_input_material.py

Shows user how to input a mineral of his/her choice and which physical values
need to be input for BurnMan to calculate Vs, Vp, Vphi and density at depth.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- how to create your own minerals


* example_compare_enstpyro.py

This example shows you how to create two materials from wt% determines the
optimum mixing between the two to match the seismic model of your choice.

requires:
- geotherms
- seismic models
- compute seismic velocities
- creating minerals

teaches:
- weight percent materials


* example_optimize_pv.py

Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- compare errors between models
- loops over models


* example_spintransition.py

This example shows the different minerals that are implemented with a spin
transition.  Minerals with spin transition can be included in
burnman/minerals.py by defining parameters for the low spin state. Regular
parameters are by definition high spin and the second set of paramaters must
be named 'self.params_LS'. This set of parameters should include a transition
pressure called 'P_LS' in GPa. This example shows the minerals for which spin
transitions are implemented.

requires:
- geotherms
- seismic models
- compute seismic velocities

teaches:
- spin transitions


