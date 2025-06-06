# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_composite_seismic_velocities
------------------------------------

This example shows how to create different minerals, how to compute seismic
velocities, and how to compare them to a seismic reference model.

There are many different ways in BurnMan to combine minerals into a
composition. Here we present a couple of examples:

1. Two minerals mixed in simple mole fractions. Can be chosen from the BurnMan
   libraries or from user defined minerals (see example_user_input_material)
2. Example with three minerals
3. Using preset solutions
4. Defining your own solution


To turn a method of mineral creation "on" the first if statement above the
method must be set to True, with all others set to False.

Note: These minerals can include a spin transition in (Mg,Fe)O, see
example_spintransition.py for explanation of how to implement this

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Composite`
* :class:`burnman.Mineral`
* :class:`burnman.Solution`

*Demonstrates:*

* Different ways to define a composite
* Using minerals and solutions
* Compare computations to seismic models

"""

import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman import minerals
from burnman.classes.solution import Solution
from burnman.classes.solutionmodel import IdealSolution


if __name__ == "__main__":
    # To compute seismic velocities and other properties, we need to supply
    # burnman with a list of minerals (phases) and their molar abundances.
    # Minerals are classes found in burnman.minerals and are derived from
    # burnman.minerals.material.
    # Here are a few ways to define phases and molar_abundances:
    # Example 1: two simple fixed minerals
    if True:
        amount_perovskite = 0.95
        rock = burnman.Composite(
            [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.periclase()],
            [amount_perovskite, 1 - amount_perovskite],
        )

    # Example 2: three materials
    if False:
        rock = burnman.Composite(
            [
                minerals.SLB_2011.fe_perovskite(),
                minerals.SLB_2011.periclase(),
                minerals.SLB_2011.stishovite(),
            ],
            [0.7, 0.2, 0.1],
        )

    # Example 3: Mixing solutions
    if False:
        # Defining a rock using a predefined solution from the mineral
        # library database.
        preset_solution = minerals.SLB_2011.mg_fe_perovskite()
        # The line below is optional to see which endmembers
        # (and in which order) are in the solution
        # print preset_solution.endmembers
        # Set molar_fraction of mg_perovskite, fe_perovskite and al_perovskite
        preset_solution.set_composition([0.9, 0.1, 0.0])
        rock = burnman.Composite(
            [preset_solution, minerals.SLB_2011.periclase()], [0.8, 0.2]
        )

    # Example 4: Defining your own solution
    if False:
        # Define a new Solution with mg and fe perovskite endmembers
        mpv = minerals.SLB_2011.mg_perovskite()
        fpv = minerals.SLB_2011.fe_perovskite()
        new_solution = Solution(
            name="New Mg-Fe bridgmanite",
            solution_model=IdealSolution(
                endmembers=[[mpv, "[Mg]SiO3"], [fpv, "[Fe]SiO3"]]
            ),
        )

        # Set molar fraction of endmembers
        new_solution.set_composition([0.9, 0.1])
        rock = burnman.Composite(
            [new_solution, minerals.SLB_2011.periclase()], [0.8, 0.2]
        )

    # seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()
    # set on how many depth slices the computations should be done
    number_of_points = 20
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    # alternatively, we could use the values where prem is defined:
    # depths = seismic_model.internal_depth_list(mindepth=700.e3,
    # maxdepth=2800.e3)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi"], depths
    )

    temperature = burnman.geotherm.brown_shankland(depths)

    print("Calculations are done for:")
    rock.debug_print()

    mat_rho, mat_vp, mat_vphi, mat_vs, mat_K, mat_G = rock.evaluate(
        ["density", "v_p", "v_phi", "v_s", "K_S", "G"], seis_p, temperature
    )

    [vs_err, vphi_err, rho_err] = burnman.utils.math.chisqr_profiles(
        [mat_vs, mat_vphi, mat_rho], [seis_vs, seis_vphi, seis_rho]
    )

    # PLOTTING
    # plot vs
    plt.subplot(2, 2, 1)
    plt.plot(
        seis_p / 1.0e9,
        mat_vs / 1.0e3,
        color="b",
        linestyle="-",
        marker="o",
        markerfacecolor="b",
        markersize=4,
        label="computation",
    )
    plt.plot(
        seis_p / 1.0e9,
        seis_vs / 1.0e3,
        color="k",
        linestyle="-",
        marker="o",
        markerfacecolor="k",
        markersize=4,
        label="reference",
    )
    plt.title("Vs (km/s)")
    plt.xlim(min(seis_p) / 1.0e9, max(seis_p) / 1.0e9)
    plt.ylim(5.1, 7.6)
    plt.legend(loc="lower right")
    plt.text(40, 7.3, "misfit= %3.3f" % vs_err)

    # plot Vphi
    plt.subplot(2, 2, 2)
    plt.plot(
        seis_p / 1.0e9,
        mat_vphi / 1.0e3,
        color="b",
        linestyle="-",
        marker="o",
        markerfacecolor="b",
        markersize=4,
    )
    plt.plot(
        seis_p / 1.0e9,
        seis_vphi / 1.0e3,
        color="k",
        linestyle="-",
        marker="o",
        markerfacecolor="k",
        markersize=4,
    )
    plt.title("Vphi (km/s)")
    plt.xlim(min(seis_p) / 1.0e9, max(seis_p) / 1.0e9)
    plt.ylim(7, 12)
    plt.text(40, 11.5, "misfit= %3.3f" % vphi_err)

    # plot density
    plt.subplot(2, 2, 3)
    plt.plot(
        seis_p / 1.0e9,
        mat_rho / 1.0e3,
        color="b",
        linestyle="-",
        marker="o",
        markerfacecolor="b",
        markersize=4,
    )
    plt.plot(
        seis_p / 1.0e9,
        seis_rho / 1.0e3,
        color="k",
        linestyle="-",
        marker="o",
        markerfacecolor="k",
        markersize=4,
    )
    plt.title("density ($\\cdot 10^3$ kg/m$^3$)")
    plt.xlim(min(seis_p) / 1.0e9, max(seis_p) / 1.0e9)
    plt.text(40, 4.3, "misfit= %3.3f" % rho_err)
    plt.xlabel("Pressure (GPa)")

    # plot geotherm
    plt.subplot(2, 2, 4)
    plt.plot(
        seis_p / 1e9,
        temperature,
        color="r",
        linestyle="-",
        marker="o",
        markerfacecolor="r",
        markersize=4,
    )
    plt.title("Geotherm (K)")
    plt.xlim(min(seis_p) / 1.0e9, max(seis_p) / 1.0e9)
    plt.xlabel("Pressure (GPa)")

    plt.savefig("output_figures/example_composition.png")
    plt.show()
