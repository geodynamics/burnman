# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2024 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
fit_spock
---------
"""

from copy import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import burnman
from burnman.utils.unitcell import molar_volume_from_unit_cell_volume
from burnman.utils.misc import pretty_string_values
from tabulate import tabulate

# Colorblind friendly colors
# https://www.nature.com/articles/nmeth.1618.pdf
colors = {
    "SPOCK": "#000000",
    "BM3": "#E69F00",
    "none": "#56B4E9",
    "MT": "#009E73",
    "Vinet": "#F0E442",
    "MACAW": "#0072B2",
    "BM4": "#D55E00",
    "RK": "#CC79A7",
}

# Loop over the gold and platinum datasets from Fratanduono et al. (2021)
for basename in ["Au", "Pt"]:
    print(f"Fitting {basename} equation of state parameters")

    # Load the data
    d = np.loadtxt(f"data/{basename}.dat")

    # Convert the data into SI units
    V = molar_volume_from_unit_cell_volume(d[:, 1], 1.0)
    P = d[:, -2] * 1.0e9
    T = P * 0.0 + 298.15
    P_unc = np.max(np.array([d[:, -1] * 1.0e9, 1.0e6 + d[:, -1] * 0.0]), axis=0)

    # Approximate the bulk modulus from the data
    KT_estimates = -np.gradient(P, np.log(V), edge_order=2)

    # Create the PTV array and covariance matrices for the data
    PTV = np.array([P, T, V]).T
    nul = 0.0 * PTV.T[0]
    PTV_covariances = np.array(
        [
            [P_unc * P_unc, nul, nul],
            [nul, nul, nul],
            [nul, nul, 1.0e-6 * V * V],
        ]
    ).T

    # Start plotting the data
    fig = plt.figure(figsize=(12, 6))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

    ax[0].errorbar(
        PTV[:, 0] / 1.0e9,
        PTV[:, 2] * 1.0e6,
        xerr=np.sqrt(PTV_covariances.T[0][0]) / 1.0e9,
        yerr=np.sqrt(PTV_covariances.T[2][2]) * 1.0e6,
        linestyle="None",
        marker="o",
        markersize=2,
        c="red",
        label="Fratanduono et al. (2021)",
        zorder=-200,
    )

    ax[1].scatter(
        PTV[:, 0] / 1.0e9,
        KT_estimates / 1.0e9,
        linestyle="None",
        marker="o",
        s=2,
        c="red",
        label="Fratanduono et al. (2021)",
        zorder=-200,
    )

    # Define the step size for the fitting routines
    deltas = {
        "V_0": 1.0e-7,
        "K_0": 1.0e8,
        "Kprime_0": 0.1,
        "Kprime_inf": 0.1,
        "Kdprime_0": 1.0e-12,
        "Kprime_prime_0": 1.0e-14,
    }

    # Define the labels for plotting
    labels = {
        "V_0": "$V_0$",
        "K_0": "$K_0$",
        "Kprime_0": "$K_0'$",
        "Kprime_inf": "$K_{{\\infty}}'$",
        "Kdprime_0": "$K_{{0}}''$",
        "Kprime_prime_0": "$K_{{0}}''$",
    }

    # Initialise the values by first fitting with the
    # Vinet equation of state, which has the benefit of being robust
    vinet_params = {
        "V_0": V[0],
        "K_0": 170.9e9,
        "Kprime_0": 5.88,
        "equation_of_state": "vinet",
    }
    fit_params = ["V_0", "K_0", "Kprime_0"]
    m = burnman.Mineral(vinet_params)
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        m,
        fit_params,
        PTV,
        PTV_covariances,
        delta_params=[deltas[p] for p in fit_params],
        verbose=False,
        param_tolerance=1.0e-8,
    )

    # Instantiate the dictionary to be filled with optimized values
    fitted_value_dict = {}

    # Define the pressure range over which to plot the models
    P_plot = np.linspace(PTV[0, 0], PTV[-1, 0], 1001)
    T_plot = 298.15 + 0.0 * P_plot

    # Loop over all the equations of state
    for sname, ename, eos, fit_params in [
        ("Vinet", "Vinet", "vinet", ["V_0", "K_0", "Kprime_0"]),
        ("BM3", "BM3", "bm3", ["V_0", "K_0", "Kprime_0"]),
        ("BM4", "BM4", "bm4", ["V_0", "K_0", "Kprime_0", "Kprime_prime_0"]),
        ("MT", "Modified Tait", "mt", ["V_0", "K_0", "Kprime_0", "Kdprime_0"]),
        ("RK", "RK", "rkprime", ["V_0", "K_0", "Kprime_0", "Kprime_inf"]),
        ("MACAW", "MACAW", "macaw", ["V_0", "K_0", "Kprime_0", "Kprime_inf"]),
        (
            "SPOCK",
            "SPOCK",
            "spock",
            ["V_0", "K_0", "Kprime_0", "Kprime_inf", "Kdprime_0"],
        ),
    ]:

        # Define some plotting values
        if sname == "SPOCK":
            linestyle = "-"
            linewidth = 2.0
            alpha = 1.0
        else:
            linestyle = "-"
            linewidth = 1.5
            alpha = 1.0

        params = copy(vinet_params)
        params["Kprime_inf"] = 2.87
        params["equation_of_state"] = eos
        params["Kdprime_0"] = -2.0 * params["Kprime_0"] / params["K_0"]
        params["Kprime_prime_0"] = -2.0 * params["Kprime_0"] / params["K_0"]

        print(f"Optimizing {ename} equation of state...")
        m = burnman.Mineral(params)
        fitted_eos = burnman.eos_fitting.fit_PTV_data(
            m,
            fit_params,
            PTV,
            PTV_covariances,
            delta_params=[deltas[p] for p in fit_params],
            verbose=False,
            param_tolerance=1.0e-8,
        )

        # Store the optimized parameters
        vals, sigs, scales = pretty_string_values(
            fitted_eos.popt,
            fitted_eos.pcov,
            extra_decimal_places=1,
            combine_value_and_sigma=True,
        )
        vals1, sigs1, scales1 = pretty_string_values(
            fitted_eos.popt,
            fitted_eos.pcov,
            extra_decimal_places=0,
            combine_value_and_sigma=True,
        )

        fitted_value_dict[ename] = {}
        for i, p in enumerate(fitted_eos.fit_params):
            if sigs[i][-2] == "1" or sigs[i][-2] == "0":
                fitted_value_dict[ename][labels[p]] = vals[i] + scales[i][1:]
            else:
                fitted_value_dict[ename][labels[p]] = vals1[i] + scales1[i][1:]
        fitted_value_dict[ename]["WSS"] = fitted_eos.WSS

        if ename == "BM4":
            fitted_value_dict[ename]["$-K_0'' K_0$"] = float(
                -m.params["Kprime_prime_0"] * m.params["K_0"]
            )

        if ename == "SPOCK" or ename == "Modified Tait":
            fitted_value_dict[ename]["$-K_0'' K_0$"] = float(
                -m.params["Kdprime_0"] * m.params["K_0"]
            )

        V, KT = m.evaluate(["V", "K_T"], P_plot, T_plot)
        ax[0].plot(
            P_plot / 1.0e9,
            V * 1.0e6,
            c=colors[sname],
            alpha=alpha,
            linestyle=linestyle,
            linewidth=linewidth,
            label=f"{ename} (WSS={fitted_eos.WSS:.2f})",
        )
        ax[1].plot(
            P_plot / 1.0e9,
            KT / 1.0e9,
            c=colors[sname],
            alpha=alpha,
            linestyle=linestyle,
            linewidth=linewidth,
            label=f"{ename}",
        )

        diff_KT = m.evaluate(["K_T"], PTV[:, 0], PTV[:, 1])[0] - KT_estimates
        diff_V = m.evaluate(["V"], PTV[:, 0], PTV[:, 1])[0] - PTV[:, 2]
        ax[2].plot(
            PTV[:, 0] / 1.0e9,
            diff_V / PTV[:, 2],
            c=colors[sname],
            alpha=alpha,
            linestyle=linestyle,
            linewidth=linewidth,
            label=f"{ename}/obs - 1",
        )
        ax[3].plot(
            PTV[:, 0] / 1.0e9,
            diff_KT / KT_estimates,
            c=colors[sname],
            alpha=alpha,
            linestyle=linestyle,
            linewidth=linewidth,
            label=f"{ename}/obs - 1",
        )

        # Create corner plots of the covariances
        fig2 = burnman.nonlinear_fitting.corner_plot(
            fitted_eos.popt, fitted_eos.pcov, fit_params
        )
        fig2[0].savefig(f"figures/{basename}_{eos}_corner_plot.pdf")

    # Set y limits
    ax[2].set_ylim(-0.0025, 0.0025)
    ax[3].set_ylim(-0.1, 0.1)

    # Add axis labels, subplot labels and legend
    for i in range(4):
        ax[i].set_xlabel("P (GPa)")

    ax[0].set_ylabel("Volume (cm$^3$/mol)")
    ax[1].set_ylabel("Bulk modulus (GPa)")
    ax[2].set_ylabel("$V_{{model}}/V_{{obs}} - 1$")
    ax[3].set_ylabel("$K_{{model}}/K_{{obs}} - 1$")

    bbox = dict(facecolor="white", edgecolor="none", alpha=1.0)
    for i, label in enumerate(["a", "b", "c", "d"]):
        ax[i].text(
            0.04,
            0.93,
            label,
            fontsize=12,
            ha="center",
            va="center",
            transform=ax[i].transAxes,
            bbox=bbox,
        )

    ax[0].legend()
    fig.set_layout_engine("tight")
    fig.savefig(f"figures/{basename}_eos_fit.pdf")

    # Create a pandas dataframe from the data and print it
    df = pd.DataFrame.from_dict(fitted_value_dict, orient="index").fillna("-")
    df = df[
        [
            "$V_0$",
            "$K_0$",
            "$K_0'$",
            "$K_{{\\infty}}'$",
            "$K_{{0}}''$",
            "$-K_0'' K_0$",
            "WSS",
        ]
    ]
    print(basename)
    print(tabulate(df, headers="keys", tablefmt="latex_raw", floatfmt=".2f"))
    plt.show()
