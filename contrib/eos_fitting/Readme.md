# A contributed file for equation of state fitting

The python scripts contained in this directory illustrate the equation of
state fitting routines in BurnMan.

The scripts contained in this directory are as follows:

eos_fitting.py
--------------
The equation of state fitting script itself. It has been written in
such a way that users can supply a different input file to
fit their own experimental data. A few lines at the beginning of the file
declare the mineral to use as a starting guess, and the parameters to optimize
in the inversion.

read_data.py
------------
Reads the formatted experimental data in test.dat.

test.dat
--------
A file containing experimental data on the enthalpy and volume of
periclase, taken from Victor and Douglas (1963), Hazen (1976) and
Dewaele et al. (2000). Full references given at the end of that file.
