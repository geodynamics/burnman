*** About

burnman is a lower mantle shear velocity generator constrained by mineral physics.

Work started by a work-group of the CIDER 2012 workshop in Santa Barbara.

see: https://github.com/tjhei/burnman

*** Requirements

- Python 2.7 (not Python 3.x)
- Python modules:
  NumPy, SciPy, matplotlib


*** Install under Ubuntu

1. sudo apt-get install python python-scipy python-numpy python-matplotlib 
2. run with "python main.py" in a shell


*** Install on a MAC

0. get Xcode
1. If you don't have Python yet, download it (for free) from python.org/download . Make sure you have the latest version (I used 2.7). To check your version of python, type the following in a terminal:
	python --version
2. Install the latest Numpy version: http://sourceforge.net/projects/numpy/files/NumPy/
3. Install the latest Scipy at http://sourceforge.net/projects/scipy/files/
4. Install the latest Matplotlib from http://sourceforge.net/projects/matplotlib/files/matplotlib/matplotlib-1.1.1/
5. Go to the main BurnMan directory and type:
	python main.py
    Figures should show up.


Problems you might run into:
	- Installing numpy/scipy/matplotlib for a different python version than the one on your 	computer
	- Having matplotlib for 32-bit instead of 64-bit (for me this got fixed by installing the very latest 	version). This will give you the error 'no matching architecture in universal wrapper'. You can check if your python is 32 or 64 bit with the following lines:
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
6. File -- Open -- find 'main.py'
7. Run the module (or press F5)
