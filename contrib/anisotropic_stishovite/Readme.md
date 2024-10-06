## Companion code to
# A model of elastic softening and second order phase transitions in anisotropic phases, with application to stishovite and post-stishovite
## R. Myhill

The python scripts contained in this directory accompany the paper
submitted by Myhill (2024).
There are comments in each of the scripts, so it is hoped that by reading
these in tandem with the original paper, the reader will get a feel
for how to use the code for their own purposes.

The scripts contained in this directory are as follows:

stishovite_data.py
------------------
This script packages all of the data in the data directory into convenient
dictionaries.

stishovite_parameters.py
------------------------
This script provides the final fitted parameters for the model,
along with the bounds used during fitting.

stishovite_model.py
-------------------
This file makes the burnman solution objects from the parameters.
The function get_models() also modifies the Zhang elasticity data
based on the model parameters, as described in the paper.

stishovite_fit_eos.py
---------------------
This script allows the user to refine the model parameters
based on the experimental data.

stishovite_model_plots.py
-------------------------
This script creates all of the plots presented in the paper.

stishovite_model_covariance_matrix.py
-------------------------------------
This script outputs the model uncertainties and correlation matrix
from the covariance matrix.

stishovite_model_Carpenter_2000.py
----------------------------------
Plots the Gibbs energy, volumes, Reuss bulk modulus and elastic moduli
from the model of Carpenter et al. (2000)

stishovite_model_Zhang_2021.py
------------------------------
Plots the Gibbs energy, volumes, Reuss bulk modulus and elastic moduli
from the model of Zhang et al. (2021)
