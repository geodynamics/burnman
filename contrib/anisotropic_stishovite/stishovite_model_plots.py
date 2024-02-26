import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm
from matplotlib.gridspec import GridSpec
from copy import deepcopy
from scipy.optimize import minimize
from stishovite_data import common_data
from stishovite_model import get_models
from burnman.tools.eos import check_anisotropic_eos_consistency
from stishovite_data import CN_GPa, P_for_CN, KNR_GPa, CN_err_GPa
from stishovite_fit_eos import transition_pressure
from burnman.tools.plot import plot_projected_elastic_properties
from burnman.utils.math import is_positive_definite
import burnman

CN_GPa_orig = deepcopy(CN_GPa)

models = get_models()
_, stishovite_scalar, stishovite_unrelaxed, stishovite_relaxed = models[0]
fP_Zhang, fP_Andrault, fP2_Zhang, fP2_Andrault = models[1:]

stishovite_relaxed.set_composition([1.0])
stishovite_unrelaxed.set_composition([0.5, 0.5])

plot_isotropic_velocities = True
plot_seismic_properties = True
low_res_seismic_properties = True
plot_relaxed = True
plot_Q_V_abc = True
plot_G = True
plot_lnabc = True

if plot_isotropic_velocities:
    pressures = np.linspace(9.0e9, 128e9, 501)
    seismic_model = burnman.seismic.PREM()
    depths = seismic_model.depth(pressures)
    temperatures = burnman.geotherm.brown_shankland(depths)

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

    fig = plt.figure()
    ax = [fig.add_subplot(1, 1, 1)]

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

    ax[0].set_xlabel("Depths (km)")
    ax[0].set_ylabel("Velocities (km/s)")
    ax[0].legend()
    fig.set_tight_layout(True)
    fig.savefig("figures/stv_isotropic_velocities.pdf")
    plt.show()

if plot_seismic_properties:
    T = 2200.0
    # for delta_P in [0.1e9]:
    for delta_P in [-40.0e9, -1.0e9, -0.1e9, 0.1e9, 1.0e9, 40.0e9]:
        print(
            f"Plotting seismic properties at {delta_P/1.e9:.2f} GPa above the transition"
        )
        P = transition_pressure(stishovite_scalar, T) + delta_P

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

        fig.set_tight_layout(True)
        fig.savefig(
            f"figures/stishovite_seismic_properties_{P/1.e9:.2f}_GPa_{int(T)}_K.pdf"
        )
        plt.show()

data = common_data()

PTV = np.concatenate(
    (
        data["PTV"]["stv"]["Ito_1974"],
        data["PTV"]["stv"]["Zhang_2021"],
        data["PTV"]["poststv"]["Zhang_2021"],
        data["PTV"]["stv"]["Andrault_2003"],
        data["PTV"]["poststv"]["Andrault_2003"],
        data["PTV"]["stv"]["Fischer_2018"],
        data["PTV"]["poststv"]["Fischer_2018"],
    )
)
PTV_err = np.concatenate(
    (
        data["PTV_err"]["stv"]["Ito_1974"],
        data["PTV_err"]["stv"]["Zhang_2021"],
        data["PTV_err"]["poststv"]["Zhang_2021"],
        data["PTV_err"]["stv"]["Andrault_2003"],
        data["PTV_err"]["poststv"]["Andrault_2003"],
        data["PTV_err"]["stv"]["Fischer_2018"],
        data["PTV_err"]["poststv"]["Fischer_2018"],
    )
)
abc = np.concatenate(
    (
        data["abc"]["stv"]["Ito_1974"],
        data["abc"]["stv"]["Zhang_2021"],
        data["abc"]["poststv"]["Zhang_2021"],
        data["abc"]["stv"]["Andrault_2003"],
        data["abc"]["poststv"]["Andrault_2003"],
        data["abc"]["stv"]["Fischer_2018"],
        data["abc"]["poststv"]["Fischer_2018"],
    )
)
abc_err = np.concatenate(
    (
        data["abc_err"]["stv"]["Ito_1974"],
        data["abc_err"]["stv"]["Zhang_2021"],
        data["abc_err"]["poststv"]["Zhang_2021"],
        data["abc_err"]["stv"]["Andrault_2003"],
        data["abc_err"]["poststv"]["Andrault_2003"],
        data["abc_err"]["stv"]["Fischer_2018"],
        data["abc_err"]["poststv"]["Fischer_2018"],
    )
)

id = np.concatenate(
    (
        np.ones(len(data["abc_err"]["stv"]["Ito_1974"])) * 1,
        np.ones(len(data["abc_err"]["stv"]["Zhang_2021"])) * 2,
        np.ones(len(data["abc_err"]["poststv"]["Zhang_2021"])) * 2,
        np.ones(len(data["abc_err"]["stv"]["Andrault_2003"])) * 3,
        np.ones(len(data["abc_err"]["poststv"]["Andrault_2003"])) * 3,
        np.ones(len(data["abc_err"]["stv"]["Fischer_2018"])) * 4,
        np.ones(len(data["abc_err"]["poststv"]["Fischer_2018"])) * 4,
    )
)

Ito_mask = id == 1
Zhang_mask = id == 2
Andrault_mask = id == 3
Fischer_mask = id == 4

PTV[Zhang_mask, 0] *= fP_Zhang + PTV[Zhang_mask, 0] * fP2_Zhang
PTV[Andrault_mask, 0] *= fP_Andrault + PTV[Andrault_mask, 0] * fP2_Andrault
P_for_CN *= fP_Zhang + P_for_CN * fP2_Zhang

print(
    f"Transition pressure at 3000 K: {transition_pressure(stishovite_scalar, 3000.0) / 1.0e9} GPa"
)

stishovite_relaxed.set_composition([1])
print(f"Consistent?: {check_anisotropic_eos_consistency(stishovite_relaxed)}")

if plot_lnabc:
    fig_lnabc = plt.figure(figsize=(8, 4))
    ax_lnabc = [fig_lnabc.add_subplot(1, 2, i) for i in range(1, 3)]

    colors = ["tab:blue", "tab:red", "tab:purple", "tab:orange"]
    labels = [
        "Ito et al. (1974)",
        "Zhang et al. (2021)",
        "Andrault et al. (2003)",
        "Fischer et al. (2018)",
    ]
    for i, mask in enumerate([Ito_mask, Zhang_mask, Andrault_mask, Fischer_mask]):
        c = colors[i]

        lnV = deepcopy(np.log(PTV[mask, 2]))
        if i == 0:
            stishovite_scalar.set_composition([0.5, 0.5])
            stishovite_scalar.set_state(1.0e5, 298.15)
            lnV += np.log(stishovite_scalar.V / PTV[mask, 2][0])

        ax_lnabc[0].scatter(lnV, np.log(abc[mask, 0]), color=c)
        ax_lnabc[0].errorbar(
            lnV,
            np.log(abc[mask, 0]),
            xerr=PTV_err[mask, 2] / PTV[mask, 2],
            yerr=abc_err[mask, 0] / abc[mask, 0],
            ls="none",
            color=c,
        )
        ax_lnabc[0].scatter(lnV, np.log(abc[mask, 1]), color=c)
        ax_lnabc[0].errorbar(
            lnV,
            np.log(abc[mask, 1]),
            xerr=PTV_err[mask, 2] / PTV[mask, 2],
            yerr=abc_err[mask, 1] / abc[mask, 1],
            ls="none",
            color=c,
        )

        ax_lnabc[0].scatter(
            lnV, 0.5 * np.log(abc[mask, 0] * abc[mask, 1]), s=5, color=c
        )
        ax_lnabc[1].scatter(lnV, np.log(abc[mask, 2]), color=c, label=labels[i])
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
    ax_lnabc[1].legend()
    fig_lnabc.set_tight_layout(True)
    fig_lnabc.savefig("figures/stv_lnabc.pdf")
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
    for i in range(len(pressures)):
        if i % 2 == 0:
            ax_G[0].plot(
                p_grid[:, i], G_xs[:, i] / 1000.0, label=f"{pressures[i]/1.e9} GPa"
            )
        else:
            ax_G[0].plot(p_grid[:, i], G_xs[:, i] / 1000.0, linestyle=":")
        imin = np.argmin(G_xs[:, i])
        G_min = G_xs[imin, i] / 1.0e3
        p_min = p_grid[imin, i]
        ax_G[0].scatter([p_min, 1.0 - p_min], [G_min, G_min])
    ax_G[0].set_xlim(-0.25, 1.25)
    ax_G[0].set_xlabel("$p_{Q1}$")
    ax_G[0].set_ylabel("$\\mathcal{G}_{xs}$ (kJ/mol)")
    ax_G[0].legend()

    # instantiate a second axes that shares the same y-axis
    ax_G2 = ax_G[0].twiny()
    ax_G2.set_xlim(-1.5, 1.5)
    ax_G2.tick_params(axis="x")
    ax_G2.set_xlabel("$Q$")
    fig_G.set_tight_layout(True)
    fig_G.savefig("figures/stv_G.pdf")
    plt.show()


if plot_Q_V_abc:
    fig_Q = plt.figure(figsize=(5, 4))
    gs = GridSpec(1, 2, width_ratios=[20, 1])  # 1 row, 2 columns
    ax_Q = [fig_Q.add_subplot(gs[i]) for i in range(2)]

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

    for i, T in enumerate(temperatures):
        print(f"Plotting Qs and volumes at {T} K")
        color = mapper.to_rgba(T)
        Pmin = np.max([1.0e5, (T - 300.0) * 1.0e7 - 10.0e9])
        pressures = np.linspace(Pmin, 120.0e9, 151)
        Ts = pressures * 0.0 + T
        (volumes, cell_params, molar_fractions) = stishovite_relaxed.evaluate(
            ["V", "cell_parameters", "molar_fractions"], pressures, Ts
        )
        stishovite_unrelaxed.set_composition([0.5, 0.5])
        (volumes_Q0, cell_params_Q0) = stishovite_unrelaxed.evaluate(
            ["V", "cell_parameters"], pressures, Ts
        )
        ax_Q[0].plot(
            pressures / 1.0e9, molar_fractions[:, 0] - molar_fractions[:, 1], c=color
        )
        a = cell_params[:, 0]
        b = cell_params[:, 1]
        c = cell_params[:, 2]
        ax_V[0].plot(pressures / 1.0e9, volumes * 1.0e6, c=color)
        ax_V[1].plot(pressures / 1.0e9, (volumes - volumes_Q0) * 1.0e6, c=color)
        ax_abc[0].plot(pressures / 1.0e9, cell_params[:, 0] * 1.0e2, c=color)
        ax_abc[0].plot(pressures / 1.0e9, cell_params[:, 1] * 1.0e2, c=color)
        ax_abc[1].plot(pressures / 1.0e9, cell_params[:, 2] * 1.0e2, c=color)

    scatter = ax_V[0].scatter(
        PTV[:, 0] / 1.0e9, PTV[:, 2] * 1.0e6, c=PTV[:, 1], cmap=cmap, norm=norm
    )
    ax_abc[0].scatter(
        PTV[:, 0] / 1.0e9, abc[:, 0] * 1.0e2, c=PTV[:, 1], cmap=cmap, norm=norm
    )
    ax_abc[0].scatter(
        PTV[:, 0] / 1.0e9, abc[:, 1] * 1.0e2, c=PTV[:, 1], cmap=cmap, norm=norm
    )
    ax_abc[1].scatter(
        PTV[:, 0] / 1.0e9, abc[:, 2] * 1.0e2, c=PTV[:, 1], cmap=cmap, norm=norm
    )

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

    # instantiate a second axes that shares the same x-axis
    color = "tab:blue"
    ax_Q2 = ax_Q[0].twinx()
    ax_Q2.scatter(pressures_GPa, angles, color=color)
    ax_Q2.errorbar(pressures_GPa, angles, yerr=angles_err, ls="none", color=color)
    ax_Q2.set_ylim(0.0, 1.5 / sol.x)
    ax_Q2.tick_params(axis="y", labelcolor=color)
    ax_Q2.set_ylabel("SiO$_6$ rotation angle ($^{\\circ}$)", color=color)

    fig_Q.colorbar(
        scatter, cax=ax_Q[1], orientation="vertical", label="Temperature (K)"
    )

    fig_V.colorbar(scatter, ax=ax_V[1], label="Temperature (K)")
    fig_abc.colorbar(scatter, ax=ax_abc[1], label="Temperature (K)")

    for i in range(2):
        ax_V[i].set_xlabel("Pressure (GPa)")
        ax_abc[i].set_xlabel("Pressure (GPa)")

    ax_Q[0].set_xlabel("Pressure (GPa)")
    ax_Q[0].set_ylabel("Q")
    ax_V[0].set_ylabel("Volume (cm$^{3}$/mol)")
    ax_V[1].set_ylabel("Excess volume (cm$^{3}$/mol)")
    ax_abc[0].set_ylabel("$a$, $b$ length (cm/mol$^{\\frac{1}{3}}$)")
    ax_abc[1].set_ylabel("$c$ length (cm/mol$^{\\frac{1}{3}}$)")

    fig_Q.set_tight_layout(True)
    fig_V.set_tight_layout(True)
    fig_abc.set_tight_layout(True)

    fig_Q.savefig("figures/stv_Q.pdf")
    fig_V.savefig("figures/stv_V.pdf")
    fig_abc.savefig("figures/stv_abc.pdf")

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

    c = ax_relaxed[1].plot(pressures / 1.0e9, K_NR / 1.0e9, label="$K_{NR}$")
    ax_relaxed[1].plot(
        pressures / 1.0e9, K_NR_u / 1.0e9, linestyle=":", c=c[0].get_color()
    )
    ax_relaxed[1].scatter(P_for_CN / 1.0e9, KNR_GPa, label="$K_{NR}$")

    ax_Q0[1].plot(pressures / 1.0e9, K_NR_Q0 / 1.0e9, label="$K_{NR} (Q0)$")
    ax_Q0[1].plot(pressures / 1.0e9, K_NR_Q1 / 1.0e9, label="$K_{NR} (Q1)$")

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
            pressures / 1.0e9, C_N[:, i, j] / 1.0e9, label=f"$C_{{N{i+1}{j+1}}}$"
        )
        color = c[0].get_color()
        ax_relaxed[axi].plot(
            pressures / 1.0e9, C_N_Q0[:, i, j] / 1.0e9, linestyle=":", c=color
        )
        ax_relaxed[axi].scatter(
            P_for_CN / 1.0e9, CN_GPa[:, i, j], label=f"{i+1}{j+1}", color=color
        )
        ax_relaxed[axi].errorbar(
            P_for_CN / 1.0e9,
            CN_GPa_orig[:, i, j],
            yerr=CN_err_GPa[:, i, j],
            alpha=0.5,
            color=color,
            ls="none",
        )
        ax_relaxed[axi].scatter(
            P_for_CN / 1.0e9, CN_GPa_orig[:, i, j], s=5, alpha=0.5, color=color
        )

        ax_Q0[axi].plot(
            pressures / 1.0e9, C_N_Q0[:, i, j] / 1.0e9, label=f"{i+1}{j+1} (Q0)"
        )
        ax_Q0[axi].plot(
            pressures / 1.0e9, C_N_Q1[:, i, j] / 1.0e9, label=f"{i+1}{j+1} (Q1)"
        )

    for i in range(3):
        ax_relaxed[i].legend()
        ax_Q0[i].legend()

        ax_relaxed[i].set_xlabel("Pressure (GPa)")
        ax_Q0[i].set_xlabel("Pressure (GPa)")

        ax_relaxed[i].set_ylabel("Modulus (GPa)")
        ax_Q0[i].set_ylabel("Modulus (GPa)")

    fig_relaxed.set_tight_layout(True)
    fig_Q0.set_tight_layout(True)

    fig_relaxed.savefig("figures/stv_relaxed.pdf")

    plt.show()
