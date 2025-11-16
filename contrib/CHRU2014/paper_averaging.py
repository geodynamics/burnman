# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_averaging
---------------

This script reproduces :cite:`Cottaar2014`, Figure 2.

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:
1. Voigt-Reuss-Hill
2. Voigt averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See :cite:`Watt1976` for explanations
of each averaging scheme.

requires:
- geotherms
- compute seismic velocities

teaches:
- averaging

"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists("burnman") and os.path.exists("../../burnman"):
    sys.path.insert(1, os.path.abspath("../.."))

import burnman
from burnman import minerals
import contrib.CHRU2014.colors as colors

if __name__ == "__main__":
    figsize = (8, 6)
    prop = {"size": 12}
    # plt.rc('text', usetex=True)
    plt.rc("font", family="sans-serif")
    figure = plt.figure(dpi=100, figsize=figsize)

    """ choose 'slb2' (finite-strain 2nd order shear modulus,
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

    amount_perovskite = 0.6

    rock = burnman.Composite(
        [minerals.SLB_2011.mg_perovskite(), minerals.SLB_2011.wuestite()],
        [amount_perovskite, 1.0 - amount_perovskite],
    )

    perovskitite = burnman.Composite([minerals.SLB_2011.mg_perovskite()], [1.0])

    periclasite = burnman.Composite([minerals.SLB_2011.wuestite()], [1.0])

    # seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()
    # set on how many depth slices the computations should be done
    number_of_points = 20
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    # alternatively, we could use the values where prem is defined:
    # depths = seismic_model.internal_depth_list()
    pressures, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi"], depths
    )

    temperatures = burnman.geotherm.brown_shankland(depths)

    print("Calculations are done for:")
    rock.debug_print()

    # calculate the seismic velocities of the rock using a whole battery of
    # averaging schemes:

    # evaluate the end members
    rho_pv, vp_pv, vs_pv, vphi_pv, K_pv, G_pv = perovskitite.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    rho_fp, vp_fp, vs_fp, vphi_fp, K_fp, G_fp = periclasite.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # Voigt Reuss Hill averaging
    rock.set_averaging_scheme(burnman.averaging_schemes.VoigtReussHill())
    rho_vrh, vp_vrh, vs_vrh, vphi_vrh, K_vrh, G_vrh = rock.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # Voigt averaging
    rock.set_averaging_scheme(burnman.averaging_schemes.Voigt())
    rho_v, vp_v, vs_v, vphi_v, K_v, G_v = rock.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # Reuss averaging
    rock.set_averaging_scheme(burnman.averaging_schemes.Reuss())
    rho_r, vp_r, vs_r, vphi_r, K_r, G_r = rock.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # Upper bound for Hashin-Shtrikman averaging
    rock.set_averaging_scheme(burnman.averaging_schemes.HashinShtrikmanUpper())
    rho_hsu, vp_hsu, vs_hsu, vphi_hsu, K_hsu, G_hsu = rock.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # Lower bound for Hashin-Shtrikman averaging
    rock.set_averaging_scheme(burnman.averaging_schemes.HashinShtrikmanLower())
    rho_hsl, vp_hsl, vs_hsl, vphi_hsl, K_hsl, G_hsl = rock.evaluate(
        ["rho", "v_p", "v_s", "v_phi", "K_eff", "G_eff"], pressures, temperatures
    )

    # linear fit
    vs_lin = vs_pv * amount_perovskite + vs_fp * (1.0 - amount_perovskite)

    # PLOTTING

    # plot vs
    ax = figure.add_subplot(1, 1, 1)
    ax.plot(
        pressures / 1.0e9,
        vs_v / 1.0e3,
        color=colors.color(0),
        linewidth=2,
        linestyle="-",
        marker="^",
        markersize=4,
        label="Voigt",
    )
    ax.plot(
        pressures / 1.0e9,
        vs_r / 1.0e3,
        color=colors.color(5),
        linewidth=2,
        linestyle="-",
        marker="v",
        markersize=4,
        label="Reuss",
    )
    ax.plot(
        pressures / 1.0e9,
        vs_vrh / 1.0e3,
        color=colors.color(1),
        linestyle="-",
        marker="*",
        markersize=6,
        label="Voigt-Reuss-Hill",
    )
    ax.fill_between(
        pressures / 1.0e9,
        vs_hsu / 1.0e3,
        vs_hsl / 1.0e3,
        facecolor="red",
        alpha=0.8,
        label="Hashin-Shtrikman bounds",
        interpolate=False,
    )

    # plt.plot(pressures/1.e9,vs_hsu/1.e3,color='r',linestyle='-',\
    #    markersize=4,label='Hashin-Shtrikman')
    # plt.plot(pressures/1.e9,vs_hsl/1.e3,color='r',linestyle='-',marker='x',\
    #    markersize=4)
    ax.plot(
        pressures / 1.0e9,
        vs_lin / 1.0e3,
        color="k",
        linewidth=2,
        linestyle="--",
        markersize=4,
        label="linear",
    )
    ax.plot(
        pressures / 1.0e9,
        vs_pv / 1.0e3,
        color=colors.color(2),
        linewidth=2,
        linestyle="-",
        marker="d",
        markersize=4,
        label="Mg Perovskite",
    )
    ax.plot(
        pressures / 1.0e9,
        vs_fp / 1.0e3,
        color=colors.color(4),
        linewidth=2,
        linestyle="-",
        marker="x",
        markersize=6,
        label="Wuestite",
    )

    ax.set_ylim(3.0, 7.5)
    ax.set_xlim(min(pressures) / 1.0e9, max(pressures) / 1.0e9)

    ax.legend(loc="lower right", ncol=2, prop=prop)

    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Shear velocity $V_s$ (km/s)")
    if "RUNNING_TESTS" not in globals():
        figure.savefig("example_averaging.pdf", bbox_inches="tight")
#    plt.show()
