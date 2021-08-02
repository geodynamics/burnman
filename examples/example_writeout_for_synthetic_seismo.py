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

import numpy as np  # Library used for general array
import matplotlib.pyplot as plt  # Library used for plotting

# Import BurnMan
import burnman_path  # adds the local burnman directory to the path
import burnman
from burnman import minerals  # import mineral library seperately

assert burnman_path  # silence pyflakes warning

if __name__ == "__main__":

    # -Defining the rocks-
    # Here we test the models for a Pyrolitic and Chondritic lower mantle.

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

    # We create an array of 20 depths at which we want to evaluate PREM, and then
    # query the seismic model for the pressure, density, P wave speed, S wave
    # speed, and bulk sound velocity at those depths
    depths = np.linspace(2890e3, 671e3, 20)

    # Here we define the lower mantle as a Layer(). .
    lower_mantle = burnman.Layer(name='Lower Mantle', radii=6371.e3 - depths)
    lower_mantle.set_temperature_mode(
        temperature_mode='adiabatic',
        temperature_top=1900.)
    lower_mantle.set_pressure_mode(pressure_mode='self-consistent',
                                   pressure_top=2.4e10, gravity_bottom=10.0)

    lower_mantle.set_material(pyrolitic_mantle)
    lower_mantle.make()

    # Writing axisem input file
    burnman.output_seismo.write_axisem_input(
        [lower_mantle], modelname='lowermantle_pyrolite', plotting=True)
    # Write mineos input file
    burnman.output_seismo.write_mineos_input(
        [lower_mantle], modelname='lowermantle_pyrolite', plotting=True)

    # change composition
    lower_mantle.set_material(chondritic_mantle)
    lower_mantle.make()
    # Writing axisem input file
    burnman.output_seismo.write_axisem_input(
        [lower_mantle], modelname='lowermantle_chondritic', plotting=True)
    # Write mineos input file
    burnman.output_seismo.write_mineos_input(
        [lower_mantle], modelname='lowermantle_chondritic', plotting=True)
