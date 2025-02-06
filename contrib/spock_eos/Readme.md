## Companion code to
# The SPOCK equation of state for condensed phases under arbitrary compression
## R. Myhill

The python scripts contained in this directory accompany the paper
submitted by Myhill (2025).
There are comments in each of the scripts, and it is hoped that by reading
these in tandem with the original paper, the reader will get a feel
for how to use the code for their own purposes.

The scripts contained in this directory are as follows:

fit_spock.py
------------
This script fits various equations of state to the Au and Pt data of
Fratanduono et al., 2021. Plots Figures 1 and 2 from the paper, as well
as corner plots demonstrating the parameter uncertainties.

Au_comparison.py
----------------
This script compares various equations of state given the standard state
and infinite pressure parameters optimised for the Au SPOCK equation of state.
Plots Figure 3.

SPOCK_EoS.xlsx
--------------
This Excel spreadsheet contains an implementation of the SPOCK equation
of state.
