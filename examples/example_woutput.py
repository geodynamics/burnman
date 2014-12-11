# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

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

import os, sys, numpy as np, matplotlib.pyplot as plt
#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":

#See bottom for output tool


    ### input variables ###
    #######################


    #specify material
    amount_perovskite = 0.95
    rock = burnman.Composite([amount_perovskite, 1.0-amount_perovskite],
                             [minerals.SLB_2005.mg_fe_perovskite(0.7),
                              minerals.SLB_2005.ferropericlase(0.5)])

    #define some pressure range
    pressures = np.arange(25e9,130e9,5e9)

    temperature = burnman.geotherm.brown_shankland(pressures)


    #Begin calculating velocities and density as depth
    print "Calculations are done for:"
    rock.debug_print()

    mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = \
        burnman.velocities_from_rock(rock, pressures, temperature, \
                                     burnman.averaging_schemes.VoigtReussHill())

    #write to file:
    output_filename = "example_woutput.txt"
    f = open(output_filename, 'wb')
    f.write("#Pressure\tTemperature\tmat_rho\tmat_vs\tmat_vp\tmat_vphi\tmat_K\tmat_G\n")

    data = zip(pressures,temperature,mat_rho, mat_vs, mat_vp, mat_vphi, mat_K, mat_G)
    np.savetxt(f, data, fmt='%.10e', delimiter='\t')


    print "\nYour data has been saved as: ",output_filename

