# BurnMan - a thermoelastic and thermodynamic toolkit for Earth and planetary sciences

## About

BurnMan is a Python library for generating thermodynamic and thermoelastic models of planetary interiors.

It began as a working-group at the 2012 CIDER workshop in Santa Barbara.

BurnMan is released under the GNU GPL v2 or newer

Homepage: http://burnman.org

Documentation: http://burnman.org/current-doc

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

* Python 2.7 (not Python 3.x)
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
   python.org/download . Make sure to use the latest version 2.x version (I
   used 2.7). To check your version of python, type the following in a
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

make Python 2.7.3 (for example) running under windows (do not use Python 3.x, but 2.7.x):

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
Here are some specific features and pitfalls on python:

* Python uses specific identation. A script might fail if a code block is not idented correctly. We use four spaces and no tabs, 
  mixing these can give trouble.
* Indices require square brackets and function or method calls parentheses (mainly different from Matlab).
* The first index of an array is 0 (e.g. x[0])
* Put dots after numbers to make them floats instead of integers (e.g. 5/3 will give 1 (rounded downward), while 5./3. will give 1.66666666667)


## Examples 
