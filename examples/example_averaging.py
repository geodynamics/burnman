# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_averaging
-----------------

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:

1. Voight-Reuss-Hill
2. Voight averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See :cite:`Watt1976` Journal of Geophysics and Space Physics for explanations
of each averaging scheme.

*Specifically uses:*

* :class:`burnman.averaging_schemes.VoigtReussHill`
* :class:`burnman.averaging_schemes.Voigt`
* :class:`burnman.averaging_schemes.Reuss`
* :class:`burnman.averaging_schemes.HashinShtrikmanUpper`
* :class:`burnman.averaging_schemes.HashinShtrikmanLower`

*Demonstrates:*

* implemented averaging schemes

"""
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import burnman
from burnman import minerals


if __name__ == "__main__":

    # Create a rock out of MgSiO3 perovskite and MgO periclase
    amount_perovskite = 0.6

    rock = burnman.Composite(
        [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.periclase()],
        [amount_perovskite, 1.0 - amount_perovskite],
    )

    perovskitite = minerals.SLB_2011.mg_perovskite()
    perovskitite.name = "Mg Perovskite"

    periclasite = minerals.SLB_2011.periclase()
    periclasite.name = "Periclase"

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
    pressures, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi"], depths
    )

    temperatures = burnman.geotherm.brown_shankland(depths)

    print("Calculations are done for:")
    rock.debug_print()

    # calculate the seismic velocities of the rock using a whole battery of
    # averaging schemes:

    # Voigt Reuss Hill / Voigt / Reuss averaging
    rock_v = rock.copy()
    rock_v.set_averaging_scheme(burnman.averaging_schemes.Voigt())
    rock_v.name = "Voigt"
    rock_r = rock.copy()
    rock_r.set_averaging_scheme(burnman.averaging_schemes.Reuss())
    rock_r.name = "Reuss"
    rock_vrh = rock.copy()
    rock_vrh.set_averaging_scheme(burnman.averaging_schemes.VoigtReussHill())
    rock_vrh.name = "Voigt-Reuss-Hill"

    # Upper/lower bound for Hashin-Shtrikman averaging
    rock_hsu = rock.copy()
    rock_hsu.set_averaging_scheme(burnman.averaging_schemes.HashinShtrikmanUpper())
    rock_hsu.name = "Hashin-Shtrikman upper"
    rock_hsl = rock.copy()
    rock_hsl.set_averaging_scheme(burnman.averaging_schemes.HashinShtrikmanLower())
    rock_hsl.name = "Hashin-Shtrikman lower"

    rocks = [rock_v, rock_r, rock_vrh, rock_hsu, rock_hsl, perovskitite, periclasite]
    markers = ["^", "v", "x", "x", "x", "x", "x"]
    colors = ["c", "k", "b", "r", "r", "y", "g"]

    v_ss = [
        rock.evaluate(["shear_wave_velocity"], pressures, temperatures)[0]
        for rock in rocks
    ]

    # PLOTTING
    # plot vs
    fig = plt.figure(figsize=(10, 7))

    for i, rock in enumerate(rocks):
        plt.plot(
            pressures / 1.0e9,
            v_ss[i] / 1.0e3,
            color=colors[i],
            linestyle="-",
            marker=markers[i],
            markersize=4,
            label=rock.name,
        )

    plt.xlim(min(pressures) / 1.0e9, max(pressures) / 1.0e9)
    plt.legend(loc="upper left", prop={"size": 11}, frameon=False)
    plt.xlabel("pressure (GPa)")
    plt.ylabel("Vs (km/s)")

    v_s_norms = [(v_s - v_ss[-1]) / (v_ss[-2] - v_ss[-1]) for v_s in v_ss]

    ax = fig.add_axes([0.58, 0.18, 0.3, 0.3])

    for i, rock in enumerate(rocks):
        plt.plot(
            pressures / 1.0e9,
            v_s_norms[i],
            color=colors[i],
            linestyle="-",
            marker=markers[i],
            markersize=4,
            label=rock.name,
        )

    ax.tick_params(labelsize=10)
    plt.title("normalized by mixture endmembers", fontsize=10)
    plt.xlim(min(pressures) / 1.0e9, max(pressures) / 1.0e9)
    plt.ylim(-0.005, 1.005)
    plt.xlabel("pressure (GPa)", fontsize=10)
    plt.ylabel("normalized Vs", fontsize=10)

    plt.savefig("output_figures/example_averaging_normalized.png")
    plt.show()
