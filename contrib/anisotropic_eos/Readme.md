## Companion code to
# An anisotropic equation of state for high pressure, high temperature applications
## R. Myhill

The python scripts contained in this directory accompany the paper
published by Myhill (2022; https://doi.org/10.31223/X5561D).
There are comments in each of the scripts, so it is hoped that by reading
these in tandem with the original paper, the reader will get a feel
for how to use the code for their own purposes.

The scripts contained in this directory are as follows:

cubic_fitting.py
----------------
Creates an AnisotropicMineral object corresponding to periclase
using a Taylor expansion in f and Pth to create Psi.

cubic_fitting_second_form.py
----------------------------
Creates an AnisotropicMineral object corresponding to periclase
using an improved formulation for Psi. The parameter starting guesses
come from the following two scripts.

cubic_fitting_second_form_atherm.py
-----------------------------------
Creates an AnisotropicMineral object corresponding to periclase
using an improved formulation for Psi.

Only the athermal parameters are refined.

cubic_fitting_second_form_only_therm.py
---------------------------------------
Creates an AnisotropicMineral object corresponding to periclase
using an improved formulation for Psi.

Only the thermal terms are refined, using the athermal parameters
from the previous script.

orthorhombic_fitting.py
------------------------
Creates an AnisotropicMineral object corresponding to San Carlos olivine
using a Taylor expansion in f and Pth to create Psi.

orthorhombic_fitting_second_form.py
-----------------------------------

Creates an AnisotropicMineral object corresponding to San Carlos olivine
using an improved formulation for Psi.

There is also a data directory that contains the data used to fit the
model parameters.
