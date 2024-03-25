## Companion code to
# An anisotropic equation of state for solid solutions
## R. Myhill

The python scripts contained in this directory accompany the paper
submitted by Myhill (2024).
There are comments in each of the scripts, so it is hoped that by reading
these in tandem with the original paper, the reader will get a feel
for how to use the code for their own purposes.

The scripts contained in this directory are as follows:

plagioclase_data.py
-------------------
This script packages all of the data in the data directory into convenient
dictionaries.

plagioclase_parameters.py
-------------------------
This script provides the final fitted parameters for the model,
along with the bounds used during fitting.

plagioclase_model.py
--------------------
This file makes the burnman solution objects from the parameters.

plagioclase_fit_eos.py
----------------------
This script allows the user to refine the model parameters
based on the experimental data.

plagioclase_model_plots.py
--------------------------
This script creates all of the plots presented in the paper.
