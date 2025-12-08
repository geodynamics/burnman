import numpy as np
import itertools
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from copy import deepcopy
from scipy.optimize import minimize
from scipy.integrate import cumulative_simpson
from stishovite_data import common_data, get_data, other_data
from stishovite_model import get_models
from stishovite_parameters import all_args
from burnman.tools.eos import check_anisotropic_eos_consistency
from stishovite_data import (
    CN_GPa,
    P_for_CN_new,
    KNR_GPa,
    CN_err_GPa,
    SN_invGPa,
    V_for_CN,
    P_for_CN,
    J2009_pressures,
    J2009_temperatures,
    J2009_dPsidfs,
    J2009_stiffness_matrices,
    C2024_pressures,
    C2024_temperatures,
    C2024_Vp,
    C2024_Vp_unc,
    C2024_Vs,
    C2024_Vs_unc,
    C2024_Vphi,
    C2024_Vphi_unc,
)
from stishovite_fit_eos import transition_pressure
from burnman.tools.plot import plot_projected_elastic_properties
from burnman.utils.math import is_positive_definite
import burnman

CN_GPa_orig = deepcopy(CN_GPa)

models = get_models(all_args)
_, stishovite_scalar, stishovite_unrelaxed, stishovite_relaxed = models

stishovite_relaxed.set_composition([1.0])
stishovite_unrelaxed.set_composition([0.5, 0.5])

plot_isotropic_velocities = True
plot_seismic_properties = False
low_res_seismic_properties = False
plot_cell_properties = False
plot_relaxed = False
plot_relaxed_SN = False
plot_Q_V_abc = False
plot_G = False
plot_lnabc = False
plot_Fischer_Andrault_splitting = False

show_plots_one_by_one = False

if plot_Fischer_Andrault_splitting:

    data = common_data()

    PTV = np.concatenate(
        (
            data["PTV"]["stv"]["Andrault_2003"],
            data["PTV"]["poststv"]["Andrault_2003"],
            data["PTV"]["stv"]["Fischer_2018"],
            data["PTV"]["poststv"]["Fischer_2018"],
        )
    )
    PTV_err = np.concatenate(
        (
            data["PTV_err"]["stv"]["Andrault_2003"],
            data["PTV_err"]["poststv"]["Andrault_2003"],
            data["PTV_err"]["stv"]["Fischer_2018"],
            data["PTV_err"]["poststv"]["Fischer_2018"],
        )
    )
    abc = np.concatenate(
        (
            data["abc"]["stv"]["Andrault_2003"],
            data["abc"]["poststv"]["Andrault_2003"],
            data["abc"]["stv"]["Fischer_2018"],
            data["abc"]["poststv"]["Fischer_2018"],
        )
    )
    abc_err = np.concatenate(
        (
            data["abc_err"]["stv"]["Andrault_2003"],
            data["abc_err"]["poststv"]["Andrault_2003"],
            data["abc_err"]["stv"]["Fischer_2018"],
            data["abc_err"]["poststv"]["Fischer_2018"],
        )
    )

    ln_abc = np.log(abc)
    P = PTV[:, 0]
    T = PTV[:, 1]
    d_ba = ln_abc[:, 1] - ln_abc[:, 0]

    mask = np.abs(d_ba) < 1.0e-10

    fig = plt.figure(figsize=(6, 4))
    ax = [fig.add_subplot(1, 1, 1)]

    scatter = plt.scatter(P / 1.0e9, T, c=d_ba, cmap="rainbow", s=50, edgecolor="k")
    ax[0].scatter(P[mask] / 1.0e9, T[mask], c="white", s=50, edgecolor="k")
    cbar = fig.colorbar(scatter, ax=ax[0], pad=0.1, label="$\\ln (b/a)$")

    args = deepcopy(all_args)
    for ddebye_0 in [args[6] + 0.6, args[6] + 5.6, args[6] + 10.6]:
        args[6] = ddebye_0
        models = get_models(args)
        _, stishovite_scalar, stishovite_unrelaxed, stishovite_relaxed = models
        debye_0 = stishovite_scalar.endmembers[0][0].params["Debye_0"]
        temperatures = np.linspace(5.0, 4000.0, 41)
        pressures = np.array(
            [transition_pressure(stishovite_scalar, T) for T in temperatures]
        )
        ax[0].plot(
            pressures / 1.0e9, temperatures, label=f"$\\Theta_D$: {debye_0:.1f} K"
        )

    ax[0].legend()
    ax[0].set_xlabel("Pressure (GPa)")
    ax[0].set_ylabel("Temperature (K)")
    ax[0].set_xlim(0.0, 125.0)
    ax[0].set_ylim(0.0, 4000.0)
    fig.set_layout_engine("tight")
    fig.savefig("figures/Fischer_Andrault_splitting.pdf")
    plt.show()


if plot_cell_properties:
    data = get_data()
    models = get_models(all_args)
    _, stishovite_scalar, stishovite_unrelaxed, stishovite_relaxed = models

    marker_cycle = itertools.cycle(Line2D.filled_markers)
    colour_cycle = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    factor = 1.00
    KTR_GPa = factor * KNR_GPa
    beta = np.sum(SN_invGPa[:, :3, :3], axis=1)
    lnV_for_CN = np.log(V_for_CN)

    fig = plt.figure(figsize=(10, 8))
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

    stishovite_relaxed.set_composition([1.0])
    J2009_V = stishovite_relaxed.evaluate(["V"], J2009_pressures, J2009_temperatures)[0]

    i = 0
    a_ln = ax[1].scatter(
        -np.log(J2009_V / J2009_V[0]),
        (J2009_dPsidfs[:, i, 0] + J2009_dPsidfs[:, i, 1] + J2009_dPsidfs[:, i, 2])
        - 1.0 / 3.0,
        label="a, b (J09)",
        marker="^",
    )
    i = 2
    c_ln = ax[1].scatter(
        -np.log(J2009_V / J2009_V[0]),
        (J2009_dPsidfs[:, i, 0] + J2009_dPsidfs[:, i, 1] + J2009_dPsidfs[:, i, 2])
        - 1.0 / 3.0,
        label="c (J09)",
        marker="^",
    )

    beta_ab = (beta[:, 0] + beta[:, 1]) / 2.0
    ax[1].scatter(
        lnV_for_CN[0] - lnV_for_CN,
        beta[:, 0] * KNR_GPa - 1.0 / 3.0,
        label="a, b (Z21)",
        c=a_ln.get_facecolor(),
    )
    ax[1].scatter(
        lnV_for_CN[0] - lnV_for_CN,
        beta[:, 1] * KNR_GPa - 1.0 / 3.0,
        c=a_ln.get_facecolor(),
    )
    ax[1].scatter(
        lnV_for_CN[0] - lnV_for_CN,
        beta[:, 2] * KNR_GPa - 1.0 / 3.0,
        label="c (Z21)",
        c=c_ln.get_facecolor(),
    )

    def dlncdlnV_linear(lnV, grad=1.2):
        return beta[0, 2] * KNR_GPa[0] + grad * (lnV_for_CN[0] - lnV)

    stishovite_relaxed.set_composition([1.0])
    stishovite_relaxed.set_state(1.0e5, 298.15)
    lnaxis0 = np.log(stishovite_relaxed.cell_parameters[:3])

    pressures_model = np.linspace(1.0e5, 150.0e9, 301)
    temperatures_model = 300.0 + 0.0 * pressures_model
    stishovite_relaxed.set_composition([1.0])

    prps = stishovite_relaxed.evaluate(
        [
            "V",
            "cell_parameters",
            "isothermal_compressibility_tensor",
            "isothermal_bulk_modulus_reuss",
        ],
        pressures_model,
        temperatures_model,
    )
    V_model = prps[0]
    cell_parameters_model = prps[1]
    isothermal_compressibility_model = prps[2]
    bulk_moduli_model = prps[3]

    lnV_model = np.log(V_model)
    for grad, linestyle in [(0.0, "--"), (1.2, ":")]:
        label = f"$d^2\\Psi_3/df^2$ = {grad}"
        ax[1].plot(
            lnV_for_CN[0] - lnV_model,
            dlncdlnV_linear(lnV_model, grad) - 1.0 / 3.0,
            c=c_ln.get_facecolor(),
            linestyle=linestyle,
            label=label,
        )

        cs = np.exp(
            lnaxis0[2]
            - cumulative_simpson(
                dlncdlnV_linear(lnV_model, grad), x=-lnV_model, initial=0
            )
        )
        ax[3].plot(
            V_model * 1.0e6,
            cs * 1.0e2,
            c=c_ln.get_facecolor(),
            linestyle=linestyle,
            label=label,
        )

    ax[1].set_ylabel("$(\\beta_{{axis}}/\\beta_{{RT}}) - 1/3$")
    ax[1].set_xlabel("$-\\ln (V/V_{{0}})$")
    ax[1].set_xlim(-0.01, np.max(lnV_for_CN[0] - lnV_for_CN) + 0.1)

    marker = next(marker_cycle)
    custom_legend = []
    labels = []
    for study in [
        "Andrault_2003",
        # "Murakami_2003",
        "Nishihara_2005",
        "Wang_2012",
        # "Grocholski_2013",
        # "Fischer_2018", # No RT data
        # "Sun_2019",
        # "Zhang_2021",
        "Zhang_2023",
    ]:
        marker = next(marker_cycle)
        colour = next(colour_cycle)

        custom_legend.append(
            (
                Line2D(
                    [0],
                    [0],
                    marker=marker,
                    color=colour,
                    markeredgecolor=colour,
                    markerfacecolor="none",
                    linestyle="",
                ),
                Line2D([0], [0], marker=marker, color=colour, linestyle=""),
            )
        )
        labels.append(f"{study.split('_')[0]} et al. ({study.split('_')[1]})")

        for phase in ["stv", "poststv"]:
            try:
                V = data[study][phase]["V"]
                a = data[study][phase]["a"]
                b = data[study][phase]["b"]
                c = data[study][phase]["c"]

                P = data[study][phase]["P"]
                T = data[study][phase]["T"]

                mask = T < 310.0
                if phase == "poststv":
                    ax[0].scatter(
                        P[mask] / 1.0e9, V[mask] * 1.0e6, marker=marker, c=colour
                    )
                    ax[2].scatter(
                        V[mask] * 1.0e6, a[mask] * 1.0e2, marker=marker, c=colour
                    )
                    ax[2].scatter(
                        V[mask] * 1.0e6, b[mask] * 1.0e2, marker=marker, c=colour
                    )
                    ax[2].scatter(
                        V[mask] * 1.0e6,
                        np.sqrt(a[mask] * b[mask]) * 1.0e2,
                        marker=marker,
                        c=colour,
                        s=3,
                    )
                    ax[3].scatter(
                        V[mask] * 1.0e6, c[mask] * 1.0e2, marker=marker, c=colour
                    )
                else:
                    ax[0].scatter(
                        P[mask] / 1.0e9,
                        V[mask] * 1.0e6,
                        marker=marker,
                        edgecolors=colour,
                        facecolors="none",
                    )
                    ax[2].scatter(
                        V[mask] * 1.0e6,
                        a[mask] * 1.0e2,
                        marker=marker,
                        edgecolors=colour,
                        facecolors="none",
                    )
                    ax[2].scatter(
                        V[mask] * 1.0e6,
                        b[mask] * 1.0e2,
                        marker=marker,
                        edgecolors=colour,
                        facecolors="none",
                    )
                    ax[3].scatter(
                        V[mask] * 1.0e6,
                        c[mask] * 1.0e2,
                        marker=marker,
                        edgecolors=colour,
                        facecolors="none",
                    )
            except KeyError:
                pass

    ax[1].plot(
        lnV_model[0] - lnV_model,
        bulk_moduli_model * isothermal_compressibility_model[:, 0, 0] - 1.0 / 3.0,
        linestyle="-",
        c=a_ln.get_facecolor(),
    )

    ax[1].plot(
        lnV_model[0] - lnV_model,
        bulk_moduli_model * isothermal_compressibility_model[:, 1, 1] - 1.0 / 3.0,
        linestyle="-",
        c=a_ln.get_facecolor(),
    )

    ax[1].plot(
        lnV_model[0] - lnV_model,
        bulk_moduli_model * isothermal_compressibility_model[:, 2, 2] - 1.0 / 3.0,
        linestyle="-",
        c=c_ln.get_facecolor(),
        label="model",
    )

    ax[2].plot(
        V_model * 1.0e6,
        cell_parameters_model[:, 0] * 1.0e2,
        linestyle="-",
        label="model",
        c=a_ln.get_facecolor(),
    )
    ax[2].plot(
        V_model * 1.0e6,
        cell_parameters_model[:, 1] * 1.0e2,
        linestyle="-",
        label="model",
        c=a_ln.get_facecolor(),
    )
    ax[3].plot(
        V_model * 1.0e6,
        cell_parameters_model[:, 2] * 1.0e2,
        linestyle="-",
        label="model",
        c=c_ln.get_facecolor(),
    )

    ln = ax[0].plot(P_for_CN / 1.0e9, V_for_CN * 1.0e6, linestyle="--", c="purple")
    custom_legend.append(ln[0])
    labels.append("Zhang et al. (2021), BLS")

    ln = ax[0].plot(P_for_CN_new / 1.0e9, V_for_CN * 1.0e6, c="purple")
    custom_legend.append(ln[0])
    labels.append("Zhang et al. (2021), BLS cal.-free")

    ln = ax[0].plot(pressures_model / 1.0e9, V_model * 1.0e6, label="model", c="black")
    custom_legend.append(ln[0])
    labels.append("model")

    ax[0].set_xlabel("P (GPa)")
    ax[0].set_ylabel("V (cm$^3$/mol)")
    ax[0].legend(
        custom_legend,
        labels,
        handler_map={tuple: HandlerTuple(ndivide=None)},
        fontsize="small",
        markerscale=0.8,
        labelspacing=0.5,
    )

    for i in range(2, 4):
        ax[i].set_xlim(14.1, 10.9)
        ax[i].set_xlabel("V (cm$^3$/mol)")
    ax[2].set_ylabel("a,b (cm/mol$^{{1/3}}$)")
    ax[3].set_ylabel("c (cm/mol$^{{1/3}}$)")
    ax[2].set_ylim(2.8 * 0.88, 2.8 * 1.01)
    ax[3].set_ylim(1.786 * 0.88, 1.786 * 1.01)

    i = 1
    ax_twin = ax[i].twiny()
    Ps = np.linspace(0.0, 120.0e9, 7)
    Vs = stishovite_relaxed.evaluate(["V"], Ps, Ps * 0.0 + 298.15)[0]
    ax_twin.set_xticks(-np.log(Vs / Vs[0]))
    ax_twin.set_xticklabels(["%d" % P for P in Ps / 1.0e9])
    ax_twin.set_xlim(ax[i].get_xlim())
    ax_twin.set_xlabel("P (model, GPa)")

    for i in range(2, 4):
        ax_twin = ax[i].twiny()
        ax_twin.set_xticks(Vs * 1.0e6)
        ax_twin.set_xticklabels(["%d" % P for P in Ps / 1.0e9])
        ax_twin.set_xlim(ax[i].get_xlim())
        ax_twin.set_xlabel("P (model, GPa)")

    ax[1].set_ylim(-1.0, 1.0)
    ax[1].legend(fontsize=8)
    fig.set_layout_engine("tight")
    fig.savefig("figures/stv_cell_properties.pdf")
    plt.show()

if plot_isotropic_velocities:
    pressures = np.linspace(9.0e9, 128e9, 501)
    seismic_model = burnman.seismic.PREM()
    depths = seismic_model.depth(pressures)
    geotherm = burnman.geotherm.BrownShankland()
    temperatures = geotherm.temperatures(depths)

    fig = plt.figure()
    ax = [fig.add_subplot(1, 1, 1)]

    for T, colour in [(300, "grey"), (1073, "black")]:
        mask = np.abs(C2024_temperatures - T) < 1.0
        ax[0].errorbar(
            seismic_model.depth(C2024_pressures[mask]) / 1.0e3,
            C2024_Vp[mask] / 1.0e3,
            C2024_Vp_unc[mask] / 1.0e3,
            c=colour,
            label=f"$V_{{VRH}}$ ({T:.0f} K; Chen et al., 2024)",
        )
        ax[0].errorbar(
            seismic_model.depth(C2024_pressures[mask]) / 1.0e3,
            C2024_Vs[mask] / 1.0e3,
            C2024_Vs_unc[mask] / 1.0e3,
            c=colour,
        )
        ax[0].errorbar(
            seismic_model.depth(C2024_pressures[mask]) / 1.0e3,
            C2024_Vphi[mask] / 1.0e3,
            C2024_Vphi_unc[mask] / 1.0e3,
            c=colour,
        )

    stishovite_relaxed.set_composition([1.0])
    v = stishovite_relaxed.evaluate(
        [
            "isentropic_shear_modulus_voigt",
            "isentropic_shear_modulus_reuss",
            "isentropic_shear_modulus_vrh",
            "isentropic_bulk_modulus_voigt",
            "isentropic_bulk_modulus_reuss",
            "isentropic_bulk_modulus_vrh",
            "density",
        ],
        pressures,
        temperatures,
    )

    G_V, G_R, G_VRH, K_V, K_R, K_VRH, rho = v

    VS_V = np.sqrt(G_V / rho)
    VS_R = np.sqrt(G_R / rho)
    VS_VRH = np.sqrt(G_VRH / rho)

    Vphi_V = np.sqrt(K_V / rho)
    Vphi_R = np.sqrt(K_R / rho)
    Vphi_VRH = np.sqrt(K_VRH / rho)

    VP_V = np.sqrt((K_V + 4.0 / 3.0 * G_V) / rho)
    VP_R = np.sqrt((K_R + 4.0 / 3.0 * G_R) / rho)
    VP_VRH = np.sqrt((K_VRH + 4.0 / 3.0 * G_VRH) / rho)

    ax[0].fill_between(
        depths / 1.0e3, VP_R / 1.0e3, VP_V / 1.0e3, alpha=0.3, color="blue"
    )
    ax[0].plot(depths / 1.0e3, VP_VRH / 1.0e3, color="blue", label="$V_P$")

    ax[0].fill_between(
        depths / 1.0e3, Vphi_R / 1.0e3, Vphi_V / 1.0e3, alpha=0.3, color="purple"
    )
    ax[0].plot(depths / 1.0e3, Vphi_VRH / 1.0e3, color="purple", label="$V_{\\Phi}$")

    ax[0].fill_between(
        depths / 1.0e3, VS_R / 1.0e3, VS_V / 1.0e3, alpha=0.3, color="red"
    )
    ax[0].plot(depths / 1.0e3, VS_VRH / 1.0e3, color="red", label="$V_S$")

    ax[0].set_ylim(
        0.0,
    )
    ax[0].set_xlim(np.min(depths) / 1.0e3, np.max(depths) / 1.0e3)

    T_labels = np.linspace(100.0, 3000.0, 30)
    spl = interpolate.splrep(temperatures, depths, s=0, k=4)
    depth_spline = interpolate.BSpline(*spl, extrapolate=False)

    ax_twin = ax[0].twiny()
    ax_twin.set_xticks(depth_spline(T_labels) / 1.0e3)
    xlabels = [f"{int(T)}" if np.abs(T % 200.0) < 1.0 else "" for T in T_labels]
    ax_twin.set_xticklabels(xlabels)
    ax_twin.set_xlim(ax[0].get_xlim())
    ax_twin.set_xlabel("Temperatures (K)", labelpad=10.0)

    P_labels = np.linspace(10, 140, 14)
    spl = interpolate.splrep(pressures, depths, s=0, k=4)
    depth_spline = interpolate.BSpline(*spl, extrapolate=False)

    ax_twin = ax[0].twiny()
    ax_twin.tick_params(axis="x", direction="in", pad=-15)
    ax_twin.set_xticks(depth_spline(P_labels * 1.0e9) / 1.0e3)
    xlabels = [f"{int(P)}" if np.abs(P % 20.0) < 1.0 else "" for P in P_labels]
    ax_twin.set_xticklabels(xlabels)
    ax_twin.set_xlim(ax[0].get_xlim())
    ax_twin.set_xlabel("Pressures (GPa)", labelpad=-30.0)

    ax[0].set_ylim(0, 18)
    ax[0].set_xlabel("Depth (km)")
    ax[0].set_ylabel("Velocity (km/s)")
    ax[0].legend(fontsize=8.0)
    fig.set_layout_engine("tight")
    fig.savefig("figures/stv_isotropic_velocities.pdf")
    if show_plots_one_by_one:
        plt.show()

if plot_seismic_properties:
    T = 2200.0
    # for delta_P in [0.1e9]:
    # for delta_P in [-40.0e9, -1.0e9, -0.1e9, 0.1e9, 1.0e9, 40.0e9]:
    for delta_P in [
        -40.0e9,
        -30.0e9,
        -20.0e9,
        -10.0e9,
        -1.0e9,
        -0.1e9,
        0.1e9,
        1.0e9,
        10.0e9,
        20.0e9,
        30.0e9,
        40.0e9,
    ]:
        print(
            f"Plotting seismic properties at {delta_P/1.e9:.2f} GPa above the transition"
        )
        P = transition_pressure(stishovite_scalar, T) + delta_P

        stishovite_unrelaxed.set_composition([0.5, 0.5])
        stishovite_relaxed.set_composition([1.0])

        stishovite_relaxed.set_state(P, T)
        stishovite_unrelaxed.set_state(P, T)

        print(
            f"Is elastic tensor positive definite?: {is_positive_definite(stishovite_relaxed.isentropic_compliance_tensor)}"
        )
        print(
            f"Isotropic Poisson ratio: {stishovite_relaxed.isentropic_isotropic_poisson_ratio:.2f}"
        )

        plot_types = [
            "vp",
            "vs1",
            "vs2",
            "linear compressibility",
            "minimum poisson ratio",
            "maximum poisson ratio",
        ]

        fig = plt.figure(figsize=(12, 7))
        ax = [fig.add_subplot(2, 3, i, projection="polar") for i in range(1, 7)]

        if low_res_seismic_properties:
            contour_sets, ticks, lines = plot_projected_elastic_properties(
                stishovite_relaxed, plot_types, ax
            )
        else:
            contour_sets, ticks, lines = plot_projected_elastic_properties(
                stishovite_relaxed, plot_types, ax, 181, 721, 100, 361
            )

        for i in range(len(contour_sets)):
            cbar = fig.colorbar(contour_sets[i], ax=ax[i], ticks=ticks[i], pad=0.1)
            cbar.add_lines(lines[i])

        fig.set_layout_engine("tight")
        fig.savefig(
            f"figures/stishovite_seismic_properties_{P/1.e9:.2f}_GPa_{int(T)}_K.pdf"
        )
        fig.savefig(
            f"figures/stishovite_seismic_properties_{P/1.e9:.2f}_GPa_{int(T)}_K.jpg"
        )
        if show_plots_one_by_one:
            plt.show()

data = common_data()
data2 = other_data()

PTV = np.concatenate(
    (
        data["PTV"]["stv"]["Ito_1974"],
        data["PTV"]["stv"]["Wang_2012"],
        data["PTV"]["stv"]["Nishihara_2005"],
        data["PTV"]["stv"]["Zhang_2021"],
        data["PTV"]["poststv"]["Zhang_2021"],
        data["PTV"]["stv"]["Andrault_2003"],
        data["PTV"]["poststv"]["Andrault_2003"],
        data["PTV"]["stv"]["Fischer_2018"],
        data["PTV"]["poststv"]["Fischer_2018"],
        data2["PTV"]["poststv"]["Murakami_2003"],
        data2["PTV"]["poststv"]["Grocholski_2013"],
        data2["PTV"]["stv"]["Zhang_2023"],
        data2["PTV"]["poststv"]["Zhang_2023"],
    )
)
PTV_err = np.concatenate(
    (
        data["PTV_err"]["stv"]["Ito_1974"],
        data["PTV_err"]["stv"]["Wang_2012"],
        data["PTV_err"]["stv"]["Nishihara_2005"],
        data["PTV_err"]["stv"]["Zhang_2021"],
        data["PTV_err"]["poststv"]["Zhang_2021"],
        data["PTV_err"]["stv"]["Andrault_2003"],
        data["PTV_err"]["poststv"]["Andrault_2003"],
        data["PTV_err"]["stv"]["Fischer_2018"],
        data["PTV_err"]["poststv"]["Fischer_2018"],
        data2["PTV_err"]["poststv"]["Murakami_2003"],
        data2["PTV_err"]["poststv"]["Grocholski_2013"],
        data2["PTV_err"]["stv"]["Zhang_2023"],
        data2["PTV_err"]["poststv"]["Zhang_2023"],
    )
)
abc = np.concatenate(
    (
        data["abc"]["stv"]["Ito_1974"],
        data["abc"]["stv"]["Wang_2012"],
        data["abc"]["stv"]["Nishihara_2005"],
        data["abc"]["stv"]["Zhang_2021"],
        data["abc"]["poststv"]["Zhang_2021"],
        data["abc"]["stv"]["Andrault_2003"],
        data["abc"]["poststv"]["Andrault_2003"],
        data["abc"]["stv"]["Fischer_2018"],
        data["abc"]["poststv"]["Fischer_2018"],
        data2["abc"]["poststv"]["Murakami_2003"],
        data2["abc"]["poststv"]["Grocholski_2013"],
        data2["abc"]["stv"]["Zhang_2023"],
        data2["abc"]["poststv"]["Zhang_2023"],
    )
)
abc_err = np.concatenate(
    (
        data["abc_err"]["stv"]["Ito_1974"],
        data["abc_err"]["stv"]["Wang_2012"],
        data["abc_err"]["stv"]["Nishihara_2005"],
        data["abc_err"]["stv"]["Zhang_2021"],
        data["abc_err"]["poststv"]["Zhang_2021"],
        data["abc_err"]["stv"]["Andrault_2003"],
        data["abc_err"]["poststv"]["Andrault_2003"],
        data["abc_err"]["stv"]["Fischer_2018"],
        data["abc_err"]["poststv"]["Fischer_2018"],
        data2["abc_err"]["poststv"]["Murakami_2003"],
        data2["abc_err"]["poststv"]["Grocholski_2013"],
        data2["abc_err"]["stv"]["Zhang_2023"],
        data2["abc_err"]["poststv"]["Zhang_2023"],
    )
)

is_stv = np.concatenate(
    (
        data["abc"]["stv"]["Ito_1974"] > 0.0,
        data["abc"]["stv"]["Wang_2012"] > 0.0,
        data["abc"]["stv"]["Nishihara_2005"] > 0.0,
        data["abc"]["stv"]["Zhang_2021"] > 0.0,
        data["abc"]["poststv"]["Zhang_2021"] < 0.0,
        data["abc"]["stv"]["Andrault_2003"] > 0.0,
        data["abc"]["poststv"]["Andrault_2003"] < 0.0,
        data["abc"]["stv"]["Fischer_2018"] > 0.0,
        data["abc"]["poststv"]["Fischer_2018"] < 0.0,
        data2["abc"]["poststv"]["Murakami_2003"] < 0.0,
        data2["abc"]["poststv"]["Grocholski_2013"] < 0.0,
        data2["abc"]["stv"]["Zhang_2023"] > 0.0,
        data2["abc"]["poststv"]["Zhang_2023"] < 0.0,
    )
)
is_stv = is_stv[:, 0]

id = np.concatenate(
    (
        np.ones(len(data["abc_err"]["stv"]["Ito_1974"])) * 1,
        np.ones(len(data["abc_err"]["stv"]["Wang_2012"])) * 4,
        np.ones(len(data["abc_err"]["stv"]["Nishihara_2005"])) * 3,
        np.ones(len(data["abc_err"]["stv"]["Zhang_2021"])) * 6,
        np.ones(len(data["abc_err"]["poststv"]["Zhang_2021"])) * 6,
        np.ones(len(data["abc_err"]["stv"]["Andrault_2003"])) * 2,
        np.ones(len(data["abc_err"]["poststv"]["Andrault_2003"])) * 2,
        np.ones(len(data["abc_err"]["stv"]["Fischer_2018"])) * 5,
        np.ones(len(data["abc_err"]["poststv"]["Fischer_2018"])) * 5,
        np.ones(len(data2["abc_err"]["poststv"]["Murakami_2003"])) * 7,
        np.ones(len(data2["abc_err"]["poststv"]["Grocholski_2013"])) * 8,
        np.ones(len(data2["abc_err"]["stv"]["Zhang_2023"])) * 9,
        np.ones(len(data2["abc_err"]["poststv"]["Zhang_2023"])) * 9,
    )
)

Ito_mask = id == 1
Andrault_mask = id == 2
Nishihara_mask = id == 3
Wang_mask = id == 4
Fischer_mask = id == 5
Zhang_mask = id == 6
Murakami_mask = id == 7
Grocholski_mask = id == 8
Zhang_mask2 = id == 9

id_labels = ["I74", "A03", "N05", "W12", "F18", "Z21"]
id_masks = [id == i for i in range(1, 6)]

P_for_CN_orig = deepcopy(P_for_CN_new)

for T in [298.15, 1000.0, 2000.0, 3000.0]:
    print(
        f"Transition pressure at {T} K: {transition_pressure(stishovite_scalar, T) / 1.0e9} GPa"
    )

stishovite_relaxed.set_composition([1])
print(f"Consistent?: {check_anisotropic_eos_consistency(stishovite_relaxed)}")

if plot_lnabc:
    fig_lnabc = plt.figure(figsize=(10, 4))
    ax_lnabc = [fig_lnabc.add_subplot(1, 2, i) for i in range(1, 3)]

    labels = [
        "Zhang et al. (2023)",
        "Zhang et al. (2021)",
        "Fischer et al. (2018)",
        "Grocholski et al. (2013)",
        "Wang et al. (2012)",
        "Nishihara et al. (2005)",
        "Murakami et al. (2003)",
        "Andrault et al. (2003)",
        "Ito et al. (1974)",
    ]
    marker_cycle = itertools.cycle(Line2D.filled_markers)
    colour_cycle = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    for i, mask in enumerate(
        [
            Zhang_mask2,
            Zhang_mask,
            Fischer_mask,
            Grocholski_mask,
            Wang_mask,
            Nishihara_mask,
            Murakami_mask,
            Andrault_mask,
            Ito_mask,
        ]
    ):
        c = next(colour_cycle)
        marker = next(marker_cycle)

        lnV = deepcopy(np.log(PTV[mask, 2]))
        if i == 0:
            stishovite_scalar.set_composition([0.5, 0.5])
            stishovite_scalar.set_state(1.0e5, 298.15)
            lnV += np.log(stishovite_scalar.V / PTV[mask, 2][0])
        s = 15
        ax_lnabc[0].scatter(lnV, np.log(abc[mask, 0]), color=c, marker=marker, s=s)
        ax_lnabc[0].errorbar(
            lnV,
            np.log(abc[mask, 0]),
            xerr=PTV_err[mask, 2] / PTV[mask, 2],
            yerr=abc_err[mask, 0] / abc[mask, 0],
            ls="none",
            color=c,
        )
        ax_lnabc[0].scatter(lnV, np.log(abc[mask, 1]), color=c, marker=marker, s=s)
        ax_lnabc[0].errorbar(
            lnV,
            np.log(abc[mask, 1]),
            xerr=PTV_err[mask, 2] / PTV[mask, 2],
            yerr=abc_err[mask, 1] / abc[mask, 1],
            ls="none",
            color=c,
        )

        ax_lnabc[0].scatter(
            lnV, 0.5 * np.log(abc[mask, 0] * abc[mask, 1]), s=2, color=c, marker=marker
        )
        ax_lnabc[1].scatter(
            lnV, np.log(abc[mask, 2]), color=c, label=labels[i], marker=marker, s=s
        )
        ax_lnabc[1].errorbar(
            lnV,
            np.log(abc[mask, 2]),
            xerr=PTV_err[mask, 2] / PTV[mask, 2],
            yerr=abc_err[mask, 2] / abc[mask, 2],
            ls="none",
            color=c,
        )

    for i in range(2):
        ax_lnabc[i].set_xlabel("$\\ln(V)$")

    ax_lnabc[0].set_ylabel("$\\ln(a)$, $\\ln(b)$ (cm/mol$^{\\frac{1}{3}}$)")
    ax_lnabc[1].set_ylabel("$\\ln(c)$ (cm/mol$^{\\frac{1}{3}}$)")
    handles, labels = ax_lnabc[1].get_legend_handles_labels()
    ax_lnabc[1].legend(handles[::-1], labels[::-1], fontsize=8)
    fig_lnabc.set_layout_engine("tight")
    fig_lnabc.savefig("figures/stv_lnabc.pdf")
    if show_plots_one_by_one:
        plt.show()


if plot_G:
    fig_G = plt.figure(figsize=(5, 4))
    ax_G = [fig_G.add_subplot(1, 1, 1)]
    pressures = np.linspace(40.0e9, 100.0e9, 13)
    ps = np.linspace(-0.25, 1.25, 201)
    pressure_grid, p_grid = np.meshgrid(pressures, ps)
    molar_fractions = np.moveaxis(np.array([p_grid, 1.0 - p_grid]), 0, 2)
    T_grid = pressure_grid * 0.0 + 298.15
    G_xs = stishovite_scalar.evaluate(
        ["excess_gibbs"], pressure_grid, T_grid, molar_fractions
    )[0]

    n_lines = len(pressures)
    cmap = matplotlib.colormaps["rainbow"]

    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))

    for i in range(len(pressures)):
        if i % 2 == 0:
            ax_G[0].plot(
                p_grid[:, i],
                G_xs[:, i] / 1000.0,
                label=f"{pressures[i]/1.e9} GPa",
                color=colors[i],
            )
        else:
            ax_G[0].plot(
                p_grid[:, i], G_xs[:, i] / 1000.0, linestyle=":", color=colors[i]
            )
        imin = np.argmin(G_xs[:, i])
        G_min = G_xs[imin, i] / 1.0e3
        p_min = p_grid[imin, i]
        ax_G[0].scatter([p_min, 1.0 - p_min], [G_min, G_min], color=colors[i])
    ax_G[0].set_xlim(-0.25, 1.25)
    ax_G[0].set_xlabel("$p_{Q1}$")
    ax_G[0].set_ylabel("$\\mathcal{G}_{xs}$ (kJ/mol)")
    ax_G[0].legend()

    # instantiate a second axes that shares the same y-axis
    ax_G2 = ax_G[0].twiny()
    ax_G2.set_xlim(-1.5, 1.5)
    ax_G2.tick_params(axis="x")
    ax_G2.set_xlabel("$Q$")
    fig_G.set_layout_engine("tight")
    fig_G.savefig("figures/stv_G.pdf")
    if show_plots_one_by_one:
        plt.show()


if plot_Q_V_abc:
    marker_cycle = itertools.cycle(Line2D.filled_markers)
    colour_cycle = itertools.cycle(plt.rcParams["axes.prop_cycle"].by_key()["color"])

    fig_Q = plt.figure(figsize=(5, 4))
    # gs = GridSpec(1, 2, width_ratios=[20, 1])  # 1 row, 2 columns
    # ax_Q = [fig_Q.add_subplot(gs[i]) for i in range(2)]
    ax_Q = [fig_Q.add_subplot(1, 1, 1)]

    fig_V = plt.figure(figsize=(8, 4))
    ax_V = [fig_V.add_subplot(1, 2, i) for i in range(1, 3)]

    fig_abc = plt.figure(figsize=(8, 4))
    ax_abc = [fig_abc.add_subplot(1, 2, i) for i in range(1, 3)]

    Tmin = np.min(PTV[:, 1])
    Tmax = 4000.0

    cmap = "rainbow"
    norm = matplotlib.colors.Normalize(vmin=Tmin, vmax=Tmax, clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cmap)

    temperatures = np.array([4000.0, 3000.0, 2000.0, 1000.0, 298.15])
    stishovite_relaxed.set_composition([1.0])
    stishovite_relaxed.set_state(1.0e5, 300)

    for i, T in enumerate(temperatures):
        print(f"Plotting Qs and volumes at {T} K")
        color = mapper.to_rgba(T)
        Pmin = np.max([1.0e5, (T - 300.0) * 1.0e7 - 10.0e9])
        pressures = np.linspace(Pmin, 150.0e9, 151)
        Ts = pressures * 0.0 + T
        (volumes, cell_params, molar_fractions) = stishovite_relaxed.evaluate(
            ["V", "cell_parameters", "molar_fractions"], pressures, Ts
        )
        stishovite_unrelaxed.set_composition([0.5, 0.5])
        (volumes_Q0, cell_params_Q0) = stishovite_unrelaxed.evaluate(
            ["V", "cell_parameters"], pressures, Ts
        )
        ax_Q[0].plot(
            pressures / 1.0e9,
            molar_fractions[:, 0] - molar_fractions[:, 1],
            c=color,
            label=f"{T} K",
        )
        a = cell_params[:, 0]
        b = cell_params[:, 1]
        c = cell_params[:, 2]
        ax_V[0].plot(pressures / 1.0e9, volumes * 1.0e6, c=color, label=f"{T} K")
        ax_V[1].plot(
            pressures / 1.0e9, (volumes - volumes_Q0) * 1.0e6, c=color, label=f"{T} K"
        )
        ax_abc[0].plot(
            pressures / 1.0e9, cell_params[:, 0] * 1.0e2, c=color, label=f"{T} K"
        )
        ax_abc[0].plot(pressures / 1.0e9, cell_params[:, 1] * 1.0e2, c=color)
        ax_abc[1].plot(
            pressures / 1.0e9, cell_params[:, 2] * 1.0e2, c=color, label=f"{T} K"
        )

    temperatures = np.linspace(298.15, 4000.0, 41)
    pressures = np.array(
        [transition_pressure(stishovite_scalar, T) for T in temperatures]
    )
    volumes = stishovite_scalar.evaluate(["V"], pressures, temperatures)[0]
    ax_V[0].plot(pressures / 1.0e9, volumes * 1.0e6, c="black", label="transition")

    for i in range(4):
        marker = next(marker_cycle)

    custom_legend = []
    labels = []

    for i in range(5):  # don't plot Zhang et al (2021)
        label = id_labels[i]
        mask = id_masks[i]
        marker = next(marker_cycle)
        s = 15

        # reversing the order makes for a better choice of label colours
        scatter_V = ax_V[0].scatter(
            PTV[mask, 0][::-1] / 1.0e9,
            PTV[mask, 2][::-1] * 1.0e6,
            c=PTV[mask, 1][::-1],
            cmap=cmap,
            norm=norm,
            s=s,
            marker=marker,
            label=label,
        )
        ax_abc[0].scatter(
            PTV[mask, 0][::-1] / 1.0e9,
            abc[mask, 0][::-1] * 1.0e2,
            c=PTV[mask, 1][::-1],
            cmap=cmap,
            norm=norm,
            s=s,
            marker=marker,
            label=label,
        )
        ax_abc[0].scatter(
            PTV[mask, 0][::-1] / 1.0e9,
            abc[mask, 1][::-1] * 1.0e2,
            c=PTV[mask, 1][::-1],
            cmap=cmap,
            norm=norm,
            s=s,
            marker=marker,
        )
        scatter_c = ax_abc[1].scatter(
            PTV[mask, 0][::-1] / 1.0e9,
            abc[mask, 2][::-1] * 1.0e2,
            c=PTV[mask, 1][::-1],
            cmap=cmap,
            norm=norm,
            s=s,
            marker=marker,
            label=label,
        )

    ax_Q[0].legend(fontsize=8.0)
    ax_abc[0].legend(fontsize=8.0)
    ax_abc[1].legend(fontsize=8.0)
    ax_V[0].legend(fontsize=8.0)
    ax_V[1].legend(fontsize=8.0)

    ax_Q[0].set_ylim(0.0, 1.5)

    # plot rotation angles on Q plot
    d = np.loadtxt("data/Zhang_2021_stishovite_angles.dat", unpack=True)
    pressures_GPa = d[0]  # *(fP_Zhang + d[0]*1.e9*fP2_Zhang)
    Ts = pressures_GPa * 0.0 + 298.15
    angles, angles_err = [d[-2], d[-1]]
    molar_fractions = stishovite_relaxed.evaluate(
        ["molar_fractions"], pressures_GPa * 1.0e9, Ts
    )[0]
    Q_model = molar_fractions[:, 0] - molar_fractions[:, 1]

    def Q_misfit(args):
        return np.sum(np.power((angles * args[0] - Q_model) / angles_err, 2.0))

    sol = minimize(Q_misfit, [0.0])
    print(f"scaling parameter for Q to theta={1./sol.x[0]}")
    # instantiate a second axes that shares the same x-axis
    color = "tab:blue"
    ax_Q2 = ax_Q[0].twinx()
    ax_Q2.scatter(pressures_GPa, angles, color=color)
    ax_Q2.errorbar(pressures_GPa, angles, yerr=angles_err, ls="none", color=color)
    ax_Q2.set_ylim(0.0, 1.5 / sol.x)
    ax_Q2.tick_params(axis="y", labelcolor=color)
    ax_Q2.set_ylabel("SiO$_6$ rotation angle ($^{\\circ}$)", color=color)

    # fig_Q.colorbar(
    #     scatter, cax=ax_Q[1], orientation="vertical", label="Temperature (K)"
    # )

    fig_V.colorbar(scatter_V, ax=ax_V[1], label="Temperature (K)")
    fig_abc.colorbar(scatter_c, ax=ax_abc[1], label="Temperature (K)")

    for i in range(2):
        ax_V[i].set_xlabel("Pressure (GPa)")
        ax_abc[i].set_xlabel("Pressure (GPa)")

    ax_Q[0].set_xlabel("Pressure (GPa)")
    ax_Q[0].set_ylabel("Q")
    ax_V[0].set_ylabel("Volume (cm$^{3}$/mol)")
    ax_V[1].set_ylabel("Excess volume (cm$^{3}$/mol)")
    ax_abc[0].set_ylabel("$a$, $b$ length (cm/mol$^{\\frac{1}{3}}$)")
    ax_abc[1].set_ylabel("$c$ length (cm/mol$^{\\frac{1}{3}}$)")

    ax_V[0].set_xlim(0.0, 150.0)
    ax_V[1].set_xlim(0.0, 150.0)
    ax_abc[0].set_xlim(0.0, 150.0)
    ax_abc[1].set_xlim(0.0, 150.0)
    ax_Q[0].set_xlim(0.0, 150.0)

    fig_Q.set_layout_engine("tight")
    fig_V.set_layout_engine("tight")
    fig_abc.set_layout_engine("tight")

    fig_Q.savefig("figures/stv_Q.pdf")
    fig_V.savefig("figures/stv_V.pdf")
    fig_abc.savefig("figures/stv_abc.pdf")
    if show_plots_one_by_one:
        plt.show()

if plot_relaxed:
    fig_relaxed = plt.figure(figsize=(12, 4))
    fig_Q0 = plt.figure(figsize=(12, 4))

    ax_relaxed = [fig_relaxed.add_subplot(1, 3, i) for i in range(1, 4)]
    ax_Q0 = [fig_Q0.add_subplot(1, 3, i) for i in range(1, 4)]

    T = 298.15
    pressures = np.linspace(1.0e5, 70.0e9, 101)
    Ts = pressures * 0.0 + T

    (molar_fractions, C_N, K_NR) = stishovite_relaxed.evaluate(
        [
            "molar_fractions",
            "isentropic_stiffness_tensor",
            "isentropic_bulk_modulus_reuss",
        ],
        pressures,
        Ts,
    )
    (C_N_u, K_NR_u) = stishovite_relaxed.unrelaxed.evaluate(
        ["isentropic_stiffness_tensor", "isentropic_bulk_modulus_reuss"],
        pressures,
        Ts,
        molar_fractions,
    )

    stishovite_unrelaxed.set_composition([0.5, 0.5])
    (C_N_Q0, K_NR_Q0) = stishovite_unrelaxed.evaluate(
        ["isentropic_stiffness_tensor", "isentropic_bulk_modulus_reuss"], pressures, Ts
    )
    stishovite_unrelaxed.set_composition([1.0, 0.0])
    (C_N_Q1, K_NR_Q1) = stishovite_unrelaxed.evaluate(
        ["isentropic_stiffness_tensor", "isentropic_bulk_modulus_reuss"], pressures, Ts
    )

    c = ax_relaxed[1].plot(pressures / 1.0e9, K_NR / 1.0e9, label="$K_{SR}$")
    ax_relaxed[1].plot(
        pressures / 1.0e9, K_NR_u / 1.0e9, linestyle=":", c=c[0].get_color()
    )
    ax_relaxed[1].scatter(P_for_CN_new / 1.0e9, KNR_GPa, label="$K_{SR}$")

    ax_Q0[1].plot(pressures / 1.0e9, K_NR_Q0 / 1.0e9, label="$K_{SR} (Q0)$")
    ax_Q0[1].plot(pressures / 1.0e9, K_NR_Q1 / 1.0e9, label="$K_{SR} (Q1)$")

    for axi, i, j in (
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, 1),
        (1, 0, 2),
        (1, 1, 2),
        (1, 2, 2),
        (2, 3, 3),
        (2, 4, 4),
        (2, 5, 5),
    ):
        c = ax_relaxed[axi].plot(
            pressures / 1.0e9,
            C_N[:, i, j] / 1.0e9,
            label=f"$\\mathbb{{C}}_{{S{i+1}{j+1}}}$",
        )
        color = c[0].get_color()
        ax_relaxed[axi].plot(
            pressures / 1.0e9, C_N_Q0[:, i, j] / 1.0e9, linestyle=":", c=color
        )
        ax_relaxed[axi].scatter(
            J2009_pressures / 1.0e9,
            J2009_stiffness_matrices[:, i, j] / 1.0e9,
            label=f"$\\mathbb{{C}}_{{S{i+1}{j+1}}}$ (J09)",
            color=color,
            marker="^",
        )
        ax_relaxed[axi].scatter(
            P_for_CN_new / 1.0e9,
            CN_GPa[:, i, j],
            label=f"$\\mathbb{{C}}_{{S{i+1}{j+1}}}$ (Z21)",
            color=color,
        )

        # Plot original data
        ax_relaxed[axi].errorbar(
            P_for_CN_orig / 1.0e9,
            CN_GPa_orig[:, i, j],
            yerr=CN_err_GPa[:, i, j],
            alpha=0.5,
            color=color,
            ls="none",
        )
        ax_relaxed[axi].scatter(
            P_for_CN_orig / 1.0e9, CN_GPa_orig[:, i, j], s=5, alpha=0.5, color=color
        )

        ax_Q0[axi].plot(
            pressures / 1.0e9, C_N_Q0[:, i, j] / 1.0e9, label=f"{i+1}{j+1} (Q0)"
        )
        ax_Q0[axi].plot(
            pressures / 1.0e9, C_N_Q1[:, i, j] / 1.0e9, label=f"{i+1}{j+1} (Q1)"
        )

    for i in range(3):
        ax_relaxed[i].legend(fontsize=8)
        ax_Q0[i].legend()

        ax_relaxed[i].set_xlabel("Pressure (GPa)")
        ax_Q0[i].set_xlabel("Pressure (GPa)")

        ax_relaxed[i].set_ylabel("Modulus (GPa)")
        ax_Q0[i].set_ylabel("Modulus (GPa)")

    fig_relaxed.set_layout_engine("tight")
    fig_Q0.set_layout_engine("tight")

    fig_relaxed.savefig("figures/stv_relaxed.pdf")

    if show_plots_one_by_one:
        plt.show()

if plot_relaxed_SN:
    fig_relaxed = plt.figure(figsize=(12, 4))
    fig_Q0 = plt.figure(figsize=(12, 4))

    ax_relaxed = [fig_relaxed.add_subplot(1, 3, i) for i in range(1, 4)]
    ax_Q0 = [fig_Q0.add_subplot(1, 3, i) for i in range(1, 4)]

    T = 298.15
    pressures = np.linspace(1.0e5, 70.0e9, 101)
    Ts = pressures * 0.0 + T

    (molar_fractions, S_N, K_NR) = stishovite_relaxed.evaluate(
        [
            "molar_fractions",
            "isentropic_compliance_tensor",
            "isentropic_bulk_modulus_reuss",
        ],
        pressures,
        Ts,
    )
    (S_N_u, K_NR_u) = stishovite_relaxed.unrelaxed.evaluate(
        ["isentropic_compliance_tensor", "isentropic_bulk_modulus_reuss"],
        pressures,
        Ts,
        molar_fractions,
    )

    stishovite_unrelaxed.set_composition([0.5, 0.5])
    (S_N_Q0, K_NR_Q0) = stishovite_unrelaxed.evaluate(
        ["isentropic_compliance_tensor", "isentropic_bulk_modulus_reuss"], pressures, Ts
    )
    stishovite_unrelaxed.set_composition([1.0, 0.0])
    (S_N_Q1, K_NR_Q1) = stishovite_unrelaxed.evaluate(
        ["isentropic_compliance_tensor", "isentropic_bulk_modulus_reuss"], pressures, Ts
    )

    c = ax_relaxed[1].plot(pressures / 1.0e9, 1.0e9 / K_NR, label="$1/K_{SR}$")
    ax_relaxed[1].plot(
        pressures / 1.0e9, 1.0e9 / K_NR_u, linestyle=":", c=c[0].get_color()
    )
    ax_relaxed[1].scatter(P_for_CN_new / 1.0e9, 1.0 / KNR_GPa, label="$1/K_{SR}$")

    ax_Q0[1].plot(pressures / 1.0e9, 1.0e9 / K_NR_Q0, label="$1/K_{SR} (Q0)$")
    ax_Q0[1].plot(pressures / 1.0e9, 1.0e9 / K_NR_Q1, label="$1/K_{SR} (Q1)$")

    for axi, i, j in (
        (0, 0, 0),
        (0, 0, 1),
        (0, 1, 1),
        (1, 0, 2),
        (1, 1, 2),
        (1, 2, 2),
        (2, 3, 3),
        (2, 4, 4),
        (2, 5, 5),
    ):
        c = ax_relaxed[axi].plot(
            pressures / 1.0e9, K_NR * S_N[:, i, j], label=f"$S_{{S{i+1}{j+1}}}$"
        )
        color = c[0].get_color()
        ax_relaxed[axi].plot(
            pressures / 1.0e9, K_NR_Q0 * S_N_Q0[:, i, j], linestyle=":", c=color
        )
        ax_relaxed[axi].scatter(
            P_for_CN_new / 1.0e9,
            KNR_GPa * SN_invGPa[:, i, j],
            label=f"S_{{{i+1}{j+1}}} (Z21)",
            color=color,
        )
        ax_relaxed[axi].scatter(
            J2009_pressures / 1.0e9,
            J2009_dPsidfs,
            label=f"S_{{{i+1}{j+1}}} (J09)",
            color=color,
        )

        # Plot original data
        ax_relaxed[axi].scatter(
            P_for_CN_orig / 1.0e9,
            KNR_GPa * SN_invGPa[:, i, j],
            s=5,
            alpha=0.5,
            color=color,
        )

        ax_Q0[axi].plot(
            pressures / 1.0e9, K_NR_Q0 * S_N_Q0[:, i, j], label=f"{i+1}{j+1} (Q0)"
        )
        ax_Q0[axi].plot(
            pressures / 1.0e9, K_NR_Q1 * S_N_Q1[:, i, j], label=f"{i+1}{j+1} (Q1)"
        )

    for i in range(3):
        ax_relaxed[i].legend(fontsize=8)
        ax_Q0[i].legend()

        ax_relaxed[i].set_xlabel("Pressure (GPa)")
        ax_Q0[i].set_xlabel("Pressure (GPa)")

        ax_relaxed[i].set_ylabel("Compliance (1/GPa)")
        ax_Q0[i].set_ylabel("Compliance (1/GPa)")

    fig_relaxed.set_layout_engine("tight")
    fig_Q0.set_layout_engine("tight")

    fig_relaxed.savefig("figures/stv_relaxed_ST.pdf")

    if show_plots_one_by_one:
        plt.show()


if not show_plots_one_by_one:
    plt.show()
