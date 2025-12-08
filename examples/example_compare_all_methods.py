# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_compare_all_methods
---------------------------

This example demonstrates how to call each of the individual calculation
methodologies that exist within BurnMan. See below for current options. This
example calculates seismic velocity profiles for the same set of minerals and
a plot of :math:`V_s, V_\\phi` and :math:`\\rho` is produce for the user to compare each of the
different methods.

*Specifically uses:*

* :ref:`ref-eos`


*Demonstrates:*

* Each method for calculating velocity profiles currently included within BurnMan

"""

import warnings
import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman import minerals


if __name__ == "__main__":
    # Input composition.

    amount_perovskite = 0.95
    rock = burnman.Composite(
        [
            minerals.Murakami_etal_2012.fe_perovskite(),
            minerals.Murakami_etal_2012.fe_periclase_LS(),
        ],
        [amount_perovskite, 1.0 - amount_perovskite],
    )

    # (min pressure, max pressure, pressure step)
    seis_p = np.arange(25e9, 125e9, 5e9)

    # Input adiabat potential temperature
    T0 = 1500.0

    # Now we'll calculate the models by forcing the rock to use a method. The
    # preset equation of state for the Murakami_etal_2012 minerals is 'slb2'

    """ 'slb2' (finite-strain 2nd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debye 3rd order shear modulus,
        matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debye 2nd order shear modulus,
        matas et al. 2007)
    or 'bm3shear2' (birch-murnaghan 3rd order, 2nd order in shear,
        if you choose to ignore temperature
        (your choice in geotherm will not matter in this case))
    or 'bm3' (birch-murnaghan 3rd order, if you choose to ignore temperature
        (your choice in geotherm will not matter in this case))"""

    methods = ["bm3", "bm3shear2", "mgd3", "mgd2", "slb3", "slb2"]
    colors = ["r", "k", "g", "b", "y", "m"]
    markers = ["+", "x", ">", "^", "<", "v"]

    fig = plt.figure(figsize=(12, 10))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

    for m in range(len(methods)):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            rock.set_method(methods[m])
        temperature = burnman.geotherm.adiabatic_profile(seis_p, rock, T0)

        print("Calculations are done for:")
        rock.debug_print()

        mat_rho_1, mat_vs_1, mat_vphi_1 = rock.evaluate(
            ["density", "v_s", "v_phi"], seis_p, temperature
        )

        # Now let's plot the comparison. You can conversely just output to a data file
        # (see example_woutput.py)

        # plot Vs
        ax[0].plot(
            seis_p / 1.0e9,
            mat_vs_1 / 1.0e3,
            color=colors[m],
            linestyle="-",
            marker=markers[m],
            markerfacecolor=colors[m],
            markersize=4,
        )

        # plot Vphi
        ax[1].plot(
            seis_p / 1.0e9,
            mat_vphi_1 / 1.0e3,
            color=colors[m],
            linestyle="-",
            marker=markers[m],
            markerfacecolor=colors[m],
            markersize=4,
        )

        # plot density
        ax[2].plot(
            seis_p / 1.0e9,
            mat_rho_1 / 1.0e3,
            color=colors[m],
            linestyle="-",
            marker=markers[m],
            markerfacecolor=colors[m],
            markersize=4,
        )

        # plot temperature
        ax[3].plot(
            seis_p / 1.0e9,
            temperature,
            color=colors[m],
            linestyle="-",
            marker=markers[m],
            markerfacecolor=colors[m],
            markersize=4,
            label=methods[m],
        )
        ax[3].legend(loc="upper left")

    for i in range(4):
        ax[0].set_xlabel("Pressure (GPa)")

    ax[0].set_ylabel("Vs (km/s)")
    ax[1].set_ylabel("Vphi (km/s)")
    ax[2].set_ylabel("Density ($\\cdot 10^3$ kg/m^3)")
    ax[3].set_ylabel("Temperature (K)")
    fig.savefig("output_figures/example_compare_all_methods.png")
    plt.show()
