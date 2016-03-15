# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_woutput
---------------

This example explains how to perform the basic i/o of BurnMan. A method of
calculation is chosen, a composite mineral/material (see
example_composition.py for explanation of this process) is created in the
class "rock," finally a geotherm is created and seismic velocities calculated.

Post-calculation, the results are written to a simple text file to
plot/manipulate at the user's whim.

requires:
- creating minerals
- compute seismic velocities
- geotherms

teaches:
- output computed seismic data to file

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":

# See bottom for output tool

    # input variables ###
    #
    # specify material
    amount_perovskite = 0.95
    fe_pv = 0.05
    fe_pc = 0.2
    pv = minerals.SLB_2011.mg_fe_perovskite()
    pc = minerals.SLB_2011.ferropericlase()
    pv.set_composition([1. - fe_pv, fe_pv, 0.])
    pc.set_composition([1. - fe_pc, fe_pc])
    rock = burnman.Composite(
        [pv, pc], [amount_perovskite, 1.0 - amount_perovskite])

    # define some pressure range
    pressures = np.arange(25e9, 130e9, 5e9)

    temperature = burnman.geotherm.brown_shankland(pressures)

    # Begin calculating velocities and density as depth
    print("Calculations are done for:")
    rock.debug_print()

    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
        rock.evaluate(
            ['density', 'v_p', 'v_s', 'v_phi', 'K_S', 'G'], pressures, temperature)
    # write to file:
    output_filename = "example_woutput.txt"
    f = open(output_filename, 'wb')
    header = "#Pressure\tTemperature\tmat_rho\tmat_vs\tmat_vp\tmat_vphi\tmat_K\tmat_G\n"
    f.write(header.encode('utf-8'))

    data = list(
        zip(pressures, temperature, mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_G))
    np.savetxt(f, data, fmt='%.10e', delimiter='\t')

    print("\nYour data has been saved as: ", output_filename)
