# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.
"""
This example script produces an input files for AXISEM and Mineos while
replacing the lower mantle with a chosen composition

Mineos computes normal mode synthetics
Mineos Masters et al. 2007
https://geodynamics.org/cig/software/mineos/
Thanks for input on this script from Jessica Irving.

Axisem computes 2.5D synthetic seismograms
Nissen-Meyer et al. (2014), AxiSEM: broadband 3-D seismic wavefields in axisymmetric media, Solid Earth, 5, 425-445. doi:10.5194/se-5-425-2014
www.axisem.info

"""
# Import supporting libraries
# Imports to be compatible with Python2 and Python3
from __future__ import absolute_import
from __future__ import print_function
import os
import sys  # Library used to interact with your operating system
import numpy as np  # Library used for general array
import matplotlib.pyplot as plt  # Library used for plotting
# Import BurnMan
sys.path.insert(1, os.path.abspath('../'))  # add path to burnman
import burnman
from burnman import minerals  # import mineral library seperately



if __name__ == "__main__":

    #-Defining the rocks-
    #Here we test the models for a Pyrolitic and Chondritic lower mantle.

    # Perovksite solide solution
    frac_mg = 0.94
    frac_fe = 0.06
    frac_al = 0.00
    mg_fe_perovskite = minerals.SLB_2011.mg_fe_perovskite()
    mg_fe_perovskite.set_composition([frac_mg, frac_fe, frac_al])

    # ferropericlase solid solution
    frac_mg = 0.8
    frac_fe = 0.2
    mg_fe_periclase = minerals.SLB_2011.ferropericlase()
    mg_fe_periclase.set_composition([frac_mg, frac_fe])

    # Ca Perovskite
    ca_perovskite = minerals.SLB_2011.ca_perovskite()

    # Pyrolitic composition
    pyr_pv = 0.75
    pyr_fp = 0.18
    pyr_capv = 0.07
    pyrolitic_mantle = burnman.Composite(
        [mg_fe_perovskite, mg_fe_periclase, ca_perovskite], [pyr_pv, pyr_fp, pyr_capv])

    # Chondritic composition
    chon_pv = 0.88
    chon_fp = 0.05
    chon_capv = 0.07
    chondritic_mantle = burnman.Composite(
        [mg_fe_perovskite, mg_fe_periclase, ca_perovskite], [chon_pv, chon_fp, chon_capv])



    # Writing axisem input file
    burnman.output_seismo.write_axisem_input(pyrolitic_mantle, filename='axisem_pyrolite.txt', plotting=True)
    burnman.output_seismo.write_axisem_input(chondritic_mantle, filename='axisem_chondrite.txt', plotting=True)
    # Write mineous input file
    burnman.output_seismo.write_mineos_input(pyrolitic_mantle, filename='mineos_pyrolite.txt', plotting=True)
    burnman.output_seismo.write_mineos_input(chondritic_mantle, filename='mineos_chondrite.txt', plotting=True)
