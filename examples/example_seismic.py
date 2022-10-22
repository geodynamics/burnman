# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_seismic
---------------

Shows the various ways to input seismic models
(:math:`V_s, V_p, V_{\\phi}, \\rho`) as a
function of depth (or pressure) as well as different velocity model libraries
available within Burnman:

1. PREM :cite:`dziewonski1981`
2. STW105 :cite:`kustowski2008`
3. AK135 :cite:`kennett1995`
4. IASP91 :cite:`kennett1991`

This example will first calculate or read in a seismic model and plot the
model along the defined pressure range. The example also illustrates how to
import a seismic model of your choice,
here shown by importing AK135 :cite:`kennett1995`.

*Uses:*

* :doc:`seismic`



*Demonstrates:*

* Utilization of library seismic models within BurnMan
* Input of user-defined seismic models


"""
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt

import burnman
import warnings


if __name__ == "__main__":

    # List of seismic 1D models
    models = [
        burnman.seismic.PREM(),
        burnman.seismic.STW105(),
        burnman.seismic.AK135(),
        burnman.seismic.IASP91(),
    ]
    colors = ["r", "b", "m", "k"]
    # Variables to plot
    variables = ["pressure", "gravity", "v_p", "v_s", "v_phi", "density"]
    units = ["Pa", "m/s^2", "m/s", "m/s", "m/s", "kg/m^3", "Pa", "m/s^2"]

    fig = plt.figure(figsize=(10, 9))
    ax = [fig.add_subplot(3, 2, i) for i in range(1, 7)]

    # Run through models and variables
    for variable_index in range(len(variables)):
        for model_index in range(len(models)):

            # specify where we want to evaluate, here we map from pressure
            # to depth
            # 1. format p = np.arange (starting pressure, ending pressure,
            # pressure step) (in Pa)
            # p = np.arange(1.0e9,360.0e9,1.e9)
            # depths = np.array([models[model_index].depth(pr) for pr in p])
            # 2. we could also just specify some depth levels directly like
            # this
            # depths = np.arange(700e3,2800e3,100e3)
            # 3. we could also use the data points where the seismic model is
            # specified over a depth range,
            # this will bring out any discontinuities
            # this is the preferred way to plot seismic discontinuities
            # correctly
            depths = models[model_index].internal_depth_list(
                mindepth=0, maxdepth=6371e3
            )
            # now evaluate everything at the given depths levels
            # (using linear interpolation)
            # try to get and plot values for given model, if this fails the
            # variable is likely not defined for that model
            try:
                with warnings.catch_warnings(record=True) as wrn:
                    values = getattr(models[model_index], variables[variable_index])(
                        depths
                    )
                    if (len(wrn) == 1) or (len(wrn) == 2):
                        for w in wrn:
                            print(w.message)
                    elif len(wrn) > 2:
                        raise Exception("Unexpected number of warnings")

                ax[variable_index].plot(
                    depths / 1.0e3,
                    values,
                    color=colors[model_index],
                    linestyle="-",
                    label=models[model_index].__class__.__name__,
                )
            except ValueError:
                # write out warning that the variable failed for given
                # model
                print(
                    variables[variable_index]
                    + " is not defined for "
                    + models[model_index].__class__.__name__
                )

        ax[variable_index].set_title(variables[variable_index])
        if variable_index == 3:
            ax[variable_index].legend()
        if variable_index > 3:
            ax[variable_index].set_xlabel("depth in km")
        ax[variable_index].set_ylabel(units[variable_index])
        ax[variable_index].set_xticks([660, 2891, 5150])

    fig.set_tight_layout(True)
    plt.show()

    # Alternatively one is able to evaluate all the variables for a model in a
    # single line
    eval = models[0].evaluate(
        ["pressure", "gravity", "v_p", "v_s", "v_phi", "density"],
        models[0].internal_depth_list(mindepth=-1.0e3, maxdepth=6372.1e3),
    )
    pressure, gravity, v_p, v_s, v_phi, density = eval

    # The following shows how to read in your own model from a file
    # Model needs to be defined with increasing depth and decreasing radius.
    # In this case the table is switched.
    class ak135_table(burnman.seismic.SeismicTable):
        def __init__(self):
            burnman.seismic.SeismicTable.__init__(self)
            # In format: radius, pressure, density, v_p, v_s
            table = burnman.utils.misc.read_table("input_seismic/ak135_lowermantle.txt")
            table = np.array(table)
            self.table_radius = table[:, 0][::-1]
            self.table_pressure = table[:, 1][::-1]
            self.table_density = table[:, 2][::-1]
            self.table_vp = table[:, 3][::-1]
            self.table_vs = table[:, 4][::-1]

            # self.table_depth needs to be defined and needs to be increasing
            self.table_depth = self.earth_radius - self.table_radius

    ak = ak135_table()
    # specify where we want to evaluate, here we map from pressure to depth
    depths = np.linspace(700e3, 2800e3, 40)
    # now evaluate everything at the given depths levels (using interpolation)
    pressures, density, v_p, v_s, v_phi = ak.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi"], depths
    )
    # plot vs and vp and v_phi (note that v_phi is computed!)
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 2, 1)
    plt.title("ak135")
    plt.plot(depths / 1.0e3, v_p / 1.0e3, "+-r", label="v_p")
    plt.plot(depths / 1.0e3, v_s / 1.0e3, "+-b", label="v_s")
    plt.plot(depths / 1.0e3, v_phi / 1.0e3, "--g", label="v_phi")
    plt.legend(loc="lower left")
    plt.xlabel("depth in km")
    plt.ylabel("km/s")

    # plot pressure,density vs depth from prem:
    plt.subplot(1, 2, 2)
    plt.title("ak135")
    plt.plot(depths / 1.0e3, pressures / 1.0e9, "-r", label="pressure")
    plt.ylabel("GPa")
    plt.xlabel("depth in km")
    plt.legend(loc="upper left")
    plt.twinx()
    plt.ylabel("g/cc")
    plt.plot(depths / 1.0e3, density / 1.0e3, "-b", label="density")
    plt.legend(loc="lower right")
    plt.show()
