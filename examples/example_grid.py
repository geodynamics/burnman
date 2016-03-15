# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_grid
------------

This example shows how to evaluate seismic quantities on a :math:`P,T` grid.
"""
from __future__ import absolute_import
from __future__ import print_function

# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":

    rock = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                              minerals.SLB_2011.periclase()], [0.8, 0.2])

    seismic_model = burnman.seismic.PREM()

    depths = np.linspace(750e3, 2800e3, 10)
    [p] = seismic_model.evaluate(['pressure'], depths)

    # Now we get an array of temperatures at which will be used for computing
    # the seismic properties of the rock.
    T = np.linspace(1900, 2400, 15)

    print("pressures:\n", p)
    print("temperatures:\n", T)

    # turn grid into array:
    tarray = np.tile(T, len(p))
    parray = np.repeat(p, len(T))

    [vs] = rock.evaluate(['v_s'], parray, tarray)

    mat_vs = np.reshape(vs, [len(p), len(T)])

    # print mat_vs

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X, Y = np.meshgrid(p / 1e9, T)

    surf = ax.plot_surface(
        X, Y, mat_vs.transpose(), rstride=1, cstride=1, linewidth=1, cmap=cm.coolwarm)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Temperature")
    ax.set_zlabel("Vs")
    ax.view_init(22, 119)

    plt.show()
