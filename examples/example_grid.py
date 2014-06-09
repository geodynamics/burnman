# BurnMan - a lower mantle toolkit
# Copyright (C) 2012, 2013, Heister, T., Unterborn, C., Rose, I. and Cottaar, S.
# Released under GPL v2 or later.

"""
This example shows how to evaluate seismic quantities on a p,T grid.
"""

# Here we import standard python modules that are required for
# usage of BurnMan.  In particular, numpy is used for handling
# numerical arrays and mathematical operations on them, and
# matplotlib is used for generating plots of results of calculations
import os, sys, numpy as np, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm


#hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1,os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":

    rock = burnman.Composite([0.8, 0.2],
                             [minerals.SLB_2011.mg_perovskite(),
                              minerals.SLB_2011.periclase()])

    seismic_model = burnman.seismic.PREM()

    depths = np.linspace(750e3, 2800e3, 10)
    p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate_all_at(depths)

    # Now we get an array of temperatures at which will be used for computing
    # the seismic properties of the rock.
    T = np.linspace(1900,2400,15)

    print "pressures:\n", p
    print "temperatures:\n", T

    # turn grid into array:
    tarray=np.tile(T,len(p))
    parray=np.repeat(p,len(T))

    rock.set_method('slb3')

    density, vp, vs, vphi, K, G = burnman.velocities_from_rock(rock, parray, tarray)

    mat_vs = np.reshape(vs,[len(p),len(T)]);

    print mat_vs

    fig = plt.figure()
    ax = fig.gca(projection='3d')


    X,Y = np.meshgrid(p/1e9, T)
    print X.shape, Y.shape, mat_vs.shape

    surf = ax.plot_surface(X,Y, mat_vs.transpose(), rstride=1, cstride=1, linewidth=1, cmap=cm.coolwarm)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Temperature")
    ax.set_zlabel("Vs")
    ax.view_init(22, 119)

    plt.show()
