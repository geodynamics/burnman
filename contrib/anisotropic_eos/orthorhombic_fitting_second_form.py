# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
orthorhombic_fitting
--------------------

This script creates an AnisotropicMineral object corresponding to
San Carlos olivine (an orthorhombic mineral). If run_fitting is set to True,
the script uses experimental data to find the optimal anisotropic parameters.
If set to False, it uses pre-optimized parameters.
The data is used to optimize both the isotropic (volumetric) and
anisotropic parameters.

The script ends by making three plots; one with the linear and volumetric
thermal expansivities at 1 bar, one with components of the
isentropic elastic stiffness tensor at high pressure, and one with
selected seismic properties at a fixed pressure and temperature.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

import burnman

from burnman import AnisotropicMineral

from tools import print_table_for_mineral_constants_2
from burnman.tools.plot import plot_projected_elastic_properties


run_fitting = False
plot_SLB = False

formula = "Mg1.8Fe0.2SiO4"
formula = burnman.tools.chemistry.dictionarize_formula(formula)
formula_mass = burnman.tools.chemistry.formula_mass(formula)

# Define the unit cell lengths and unit cell volume.
# These are taken from Abramson et al., 1997
Z = 4.0
cell_lengths_angstrom = np.array([4.7646, 10.2296, 5.9942])
cell_lengths_0_guess = cell_lengths_angstrom * np.cbrt(
    burnman.constants.Avogadro / Z / 1.0e30
)
V_0_guess = np.prod(cell_lengths_0_guess)

ol_data = np.loadtxt("data/Mao_et_al_2015_ol.dat")
ol_1bar_lattice_data_Suzuki = np.loadtxt("data/Suzuki_1975_ol_Kenya_expansion.dat")

fo = burnman.minerals.SLB_2011.forsterite()
fa = burnman.minerals.SLB_2011.fayalite()

CN_data = np.zeros((6, 6, len(ol_data)))
CNerr_data = np.ones((6, 6, len(ol_data))) * 1.0e20
CN_data[0, 0] = ol_data.T[4]
CN_data[1, 1] = ol_data.T[6]
CN_data[2, 2] = ol_data.T[8]
CN_data[3, 3] = ol_data.T[10]
CN_data[4, 4] = ol_data.T[12]
CN_data[5, 5] = ol_data.T[14]
CN_data[0, 1] = ol_data.T[16]
CN_data[1, 0] = ol_data.T[16]
CN_data[0, 2] = ol_data.T[18]
CN_data[2, 0] = ol_data.T[18]
CN_data[1, 2] = ol_data.T[20]
CN_data[2, 1] = ol_data.T[20]
CNerr_data[0, 0] = ol_data.T[5]
CNerr_data[1, 1] = ol_data.T[7]
CNerr_data[2, 2] = ol_data.T[9]
CNerr_data[3, 3] = ol_data.T[11]
CNerr_data[4, 4] = ol_data.T[13]
CNerr_data[5, 5] = ol_data.T[15]
CNerr_data[0, 1] = ol_data.T[17]
CNerr_data[1, 0] = ol_data.T[17]
CNerr_data[0, 2] = ol_data.T[19]
CNerr_data[2, 0] = ol_data.T[19]
CNerr_data[1, 2] = ol_data.T[21]
CNerr_data[2, 1] = ol_data.T[21]

CN_data = CN_data.T
CNerr_data = CNerr_data.T
SN_data = np.linalg.inv(CN_data)
betaN_data = np.sum(SN_data[:, :3, :3], axis=(1, 2))
SNoverbetaN_data = (SN_data.T / betaN_data).T


def psi_func(f, Pth, params):
    dPsidf = (
        params["a"]
        + params["b_1"] * params["c_1"] * np.exp(params["c_1"] * f)
        + params["b_2"] * params["c_2"] * np.exp(params["c_2"] * f)
    )
    Psi = (
        0.0
        + params["a"] * f
        + params["b_1"] * (np.exp(params["c_1"] * f) - 1.0)
        + params["b_2"] * (np.exp(params["c_2"] * f) - 1.0)
        + params["d"] * Pth / 1.0e9
    )
    dPsidPth = params["d"] / 1.0e9
    return (Psi, dPsidf, dPsidPth)


def make_orthorhombic_mineral_from_parameters(x):
    # First, make the scalar model
    san_carlos_params = {
        "name": "San Carlos olivine",
        "formula": formula,
        "equation_of_state": "slb3",
        "F_0": 0.0,
        "V_0": V_0_guess,  # we overwrite this in a second
        "K_0": 1.263e11,  # Abramson et al. 1997
        "Kprime_0": 4.28,  # Abramson et al. 1997
        "Debye_0": 760.0,  # Robie, forsterite
        "grueneisen_0": 0.99282,  # Fo in SLB2011
        "q_0": 2.10672,  # Fo in SLB2011
        "G_0": 81.6e9,
        "Gprime_0": 1.46257,
        "eta_s_0": 2.29972,
        "n": 7.0,
        "molar_mass": formula_mass,
    }

    R = burnman.constants.gas_constant
    san_carlos_property_modifiers = [
        [
            "linear",
            {
                "delta_E": 0.0,
                "delta_S": 26.76 * 0.1
                - 2.0 * R * (0.1 * np.log(0.1) + 0.9 * np.log(0.9)),
                "delta_V": 0.0,
            },
        ]
    ]

    ol = burnman.Mineral(
        params=san_carlos_params, property_modifiers=san_carlos_property_modifiers
    )

    # Overwrite some properties
    y = x

    ol.params["V_0"] = y[0] * V_0_guess  # Abramson et al. 1997
    ol.params["K_0"] = y[1] * 1.263e11  # Abramson et al. 1997
    ol.params["Kprime_0"] = y[2] * 4.28  # Abramson et al. 1997
    ol.params["grueneisen_0"] = y[3] * 0.99282  # Fo in SLB2011
    ol.params["q_0"] = y[4] * 2.10672  # Fo in SLB2011
    # ol.params['Debye_0'] = x[5]*809.1703 # Fo in SLB2011 strong tendency to 0

    # Next, each of the eight independent elastic tensor component get
    # their turn.
    # We arbitrarily choose S[2,3] as the ninth component,
    # which is determined by the others.
    i = 5
    anisotropic_parameters = {
        "a": np.zeros((6, 6)),
        "b_1": np.zeros((6, 6)),
        "c_1": np.ones((6, 6)),
        "d": np.zeros((6, 6)),
        "b_2": np.zeros((6, 6)),
        "c_2": np.ones((6, 6)),
    }

    for p, q in ((1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3)):
        anisotropic_parameters["a"][p - 1, q - 1] = x[i]
        anisotropic_parameters["a"][q - 1, p - 1] = x[i]
        i = i + 1
        anisotropic_parameters["b_1"][p - 1, q - 1] = x[i]
        anisotropic_parameters["b_1"][q - 1, p - 1] = x[i]
        i = i + 1
        anisotropic_parameters["d"][p - 1, q - 1] = x[i]
        anisotropic_parameters["d"][q - 1, p - 1] = x[i]
        i = i + 1

    anisotropic_parameters["c_1"][:3, :3] = x[i]
    i = i + 1
    for j in range(3):
        anisotropic_parameters["c_1"][3 + j, 3 + j] = x[i]
        i = i + 1

    anisotropic_parameters["b_2"][3, 3] = x[i]
    i = i + 1
    anisotropic_parameters["b_2"][4, 4] = x[i]
    i = i + 1
    anisotropic_parameters["b_2"][5, 5] = x[i]
    i = i + 1
    anisotropic_parameters["c_2"][3, 3] = x[i]
    i = i + 1
    anisotropic_parameters["c_2"][4, 4] = x[i]
    i = i + 1
    anisotropic_parameters["c_2"][5, 5] = x[i]
    i = i + 1

    assert len(x) == i

    # Fill the values for the dependent element c[2,3]
    anisotropic_parameters["a"][1, 2] = (
        1.0 - np.sum(anisotropic_parameters["a"][:3, :3])
    ) / 2.0
    anisotropic_parameters["b_1"][1, 2] = (
        0.0 - np.sum(anisotropic_parameters["b_1"][:3, :3])
    ) / 2.0
    anisotropic_parameters["d"][1, 2] = (
        0.0 - np.sum(anisotropic_parameters["d"][:3, :3])
    ) / 2.0

    anisotropic_parameters["a"][2, 1] = anisotropic_parameters["a"][1, 2]
    anisotropic_parameters["b_1"][2, 1] = anisotropic_parameters["b_1"][1, 2]
    anisotropic_parameters["d"][2, 1] = anisotropic_parameters["d"][1, 2]

    cell_lengths = cell_lengths_0_guess * np.cbrt(ol.params["V_0"] / V_0_guess)
    ol_cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], 90, 90, 90]
    )

    m = AnisotropicMineral(
        ol, ol_cell_parameters, anisotropic_parameters, psi_func, orthotropic=True
    )
    return m


sol = []
if run_fitting:

    def orthorhombic_misfit(x, imin):
        m = make_orthorhombic_mineral_from_parameters(x)

        chisqr = 0.0
        try:
            for i, d in enumerate(ol_data):
                TK, PGPa, rho, rhoerr = d[:4]

                PPa = PGPa * 1.0e9

                m.set_state(PPa, TK)

                CN = m.isentropic_stiffness_tensor / 1.0e9

                chisqr += np.power((m.density / 1000.0 - rho) / rhoerr, 2.0)
                chisqr += np.sum(np.power((CN - CN_data[i]) / CNerr_data[i], 2.0))

            # Not San Carlos, fo92.3, not fo90.4
            for d in ol_1bar_lattice_data_Suzuki:
                m.set_state(1.0e5, d[0] + 273.15)  # T in C

                Y = (
                    (np.diag(m.cell_vectors) / np.diag(m.cell_vectors_0)) - 1.0
                ) * 1.0e4
                Y_expt = d[1:4]
                Y_err = 0.01 * Y_expt + 1.0
                for i in range(3):
                    chisqr += np.power((Y_expt[i] - Y[i]) / Y_err[i], 2.0)

            # if chisqr < 1500.:
            #    print(chisqr)
            # m.set_state(1.e5, 300)
            # print(np.diag(m.thermal_expansivity_tensor))

            if np.isnan(chisqr):
                print(d, "Noooo, there was a nan")
                chisqr = 1.0e7

        except Exception as e:
            print(e)
            print("There was an exception")
            chisqr = 1.0e7
        imin[0][0] += 1
        if chisqr < imin[0][1]:
            imin[0][1] = chisqr
            print(imin[0])
            print(repr(x))
        return chisqr

    guesses = np.array(
        [
            1.00263718e00,
            9.92582060e-01,
            9.94959405e-01,
            1.11986468e00,
            3.26968155e-01,
            1.36482754e00,
            -9.12627751e-01,
            1.42227141e-02,
            2.38325937e00,
            -1.60675570e00,
            1.97961677e-02,
            2.40863946e00,
            -1.76297463e00,
            1.58628282e-02,
            2.77600544e01,
            -8.63488144e01,
            -2.01661072e-01,
            1.16595606e01,
            -9.09004579e00,
            -1.89410920e-01,
            2.35013807e01,
            -5.64225909e01,
            -4.62157096e-01,
            -5.14209400e-01,
            3.78113506e-01,
            -9.03837805e-03,
            -6.29669950e-01,
            5.34427058e-01,
            -4.77601300e-03,
            1.00090172e00,
            3.13842027e-01,
            1.48687484e00,
            4.16041813e-01,
            2.19213171e-01,
            6.40447373e-01,
            2.53813715e-01,
            6.08105337e00,
            5.47261046e00,
            6.12582927e00,
        ]
    )

    i = 0
    min = 1.0e10
    # sol = minimize(orthorhombic_misfit, guesses, method='COBYLA', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 10000})
    # sol = minimize(orthorhombic_misfit, guesses, method='Powell', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 10000})
    sol = minimize(
        orthorhombic_misfit,
        guesses,
        method="Nelder-Mead",
        args=[[i, min]],
        options={"rhobeg": 0.2, "maxiter": 100000},
    )
    print(sol)

do_plotting = True
if do_plotting:
    if run_fitting:
        m = make_orthorhombic_mineral_from_parameters(sol.x)
    else:
        x = np.array(
            [
                1.00263718e00,
                9.92582060e-01,
                9.94959405e-01,
                1.11986468e00,
                3.26968155e-01,
                1.36482754e00,
                -9.12627751e-01,
                1.42227141e-02,
                2.38325937e00,
                -1.60675570e00,
                1.97961677e-02,
                2.40863946e00,
                -1.76297463e00,
                1.58628282e-02,
                2.77600544e01,
                -8.63488144e01,
                -2.01661072e-01,
                1.16595606e01,
                -9.09004579e00,
                -1.89410920e-01,
                2.35013807e01,
                -5.64225909e01,
                -4.62157096e-01,
                -5.14209400e-01,
                3.78113506e-01,
                -9.03837805e-03,
                -6.29669950e-01,
                5.34427058e-01,
                -4.77601300e-03,
                1.00090172e00,
                3.13842027e-01,
                1.48687484e00,
                4.16041813e-01,
                2.19213171e-01,
                6.40447373e-01,
                2.53813715e-01,
                6.08105337e00,
                5.47261046e00,
                6.12582927e00,
            ]
        )
        m = make_orthorhombic_mineral_from_parameters(x)

    assert burnman.tools.eos.check_anisotropic_eos_consistency(m)

    print(
        "The following parameters were used for the volumetric part of "
        f'the isotropic model: $V_0$: {m.params["V_0"]*1.e6:.5f} cm$^3$/mol, '
        f'$K_0$: {m.params["K_0"]/1.e9:.5f} GPa, '
        f'$K\'_0$: {m.params["Kprime_0"]:.5f}, '
        f'$\\Theta_0$: {m.params["Debye_0"]:.5f} K, '
        f'$\\gamma_0$: {m.params["grueneisen_0"]:.5f}, '
        f'and $q_0$: {m.params["q_0"]:.5f}.'
    )

    print_table_for_mineral_constants_2(
        m,
        ["a", "b_1", "c_1", "b_2", "c_2", "d"],
        [(1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (4, 4), (5, 5), (6, 6)],
    )

    # Plot thermal expansion figure
    fig = plt.figure(figsize=(4, 8))
    ax = [fig.add_subplot(2, 1, i) for i in range(1, 3)]

    temperatures = np.linspace(10.0, 1600.0, 101)
    alphas = np.empty((101, 4))
    extensions = np.empty((101, 3))
    vectors = np.empty((101, 4))

    labels = ["a", "b", "c", "V"]

    for i, T in enumerate(temperatures):
        m.set_state(1.0e5, T)
        alphas[i, :3] = np.diag(m.thermal_expansivity_tensor) * 1.0e5
        alphas[i, 3] = m.alpha * 1.0e5 / 3.0
        extensions[i] = (
            (np.diag(m.cell_vectors) / np.diag(m.cell_vectors_0)) - 1.0
        ) * 1.0e4
        vectors[i, :3] = np.diag(m.cell_vectors)

        vectors[i, 3] = m.V

    for i in range(4):
        label = f"$\\alpha_{{{labels[i]}}}$"
        if i == 3:
            ln = ax[0].plot(temperatures, alphas[:, i], label=label + "/3")
            ol_SLB = burnman.minerals.SLB_2011.mg_fe_olivine([0.903, 0.097])
            pressures = 1.0e5 + 0.0 * temperatures
            if plot_SLB:
                ax[0].plot(
                    temperatures,
                    ol_SLB.evaluate(["alpha"], pressures, temperatures)[0]
                    * 1.0e5
                    / 3.0,
                    label=label + "/3 (SLB2011)",
                    linestyle="--",
                    color=ln[0].get_color(),
                )

        else:
            ax[0].plot(temperatures, alphas[:, i], label=label)

    for i in range(3):
        l = ax[1].plot(temperatures, extensions[:, i], label=labels[i])
        ax[1].scatter(
            ol_1bar_lattice_data_Suzuki[:, 0] + 273.15,
            ol_1bar_lattice_data_Suzuki[:, 1 + i],
            color=l[0].get_color(),
        )

    Vthird_expansion = 1.0e4 * (
        np.power(np.prod(extensions * 1.0e-4 + 1, axis=1), 1.0 / 3.0) - 1.0
    )
    ln = ax[1].plot(temperatures, Vthird_expansion, label="$V^{1/3}$")
    ol_SLB = burnman.minerals.SLB_2011.mg_fe_olivine([0.903, 0.097])
    ol_SLB.set_state(1.0e5, 300)
    V_0 = ol_SLB.V
    pressures = 1.0e5 + 0.0 * temperatures
    if plot_SLB:
        ax[1].plot(
            temperatures,
            1.0e4
            * (
                np.power(
                    ol_SLB.evaluate(["V"], pressures, temperatures)[0] / V_0, 1.0 / 3.0
                )
                - 1.0
            ),
            label="$V^{1/3}$ (SLB2011)",
            linestyle="--",
            color=ln[0].get_color(),
        )

    Vthird_expansion = 1.0e4 * (
        np.power(
            np.prod(ol_1bar_lattice_data_Suzuki[:, 1:4] * 1.0e-4 + 1, axis=1), 1.0 / 3.0
        )
        - 1.0
    )
    ax[1].scatter(
        ol_1bar_lattice_data_Suzuki[:, 0] + 273.15,
        Vthird_expansion,
        color=ln[0].get_color(),
    )

    ax[0].set_ylim(
        0.0,
    )

    for i in range(2):
        ax[i].set_xlim(0.0, 1600.0)
        ax[i].set_xlabel("Temperature (K)")
        ax[i].legend()

    ax[0].set_ylabel("Thermal expansivity (10$^{-5}$/K)")
    ax[1].set_ylabel("Relative length change ($10^{4} (x/x_0 - 1)$)")

    fig.set_layout_engine("tight")
    fig.savefig("olivine_expansivities.pdf")
    plt.show()

    # Start plotting Cij figure
    fig = plt.figure(figsize=(12, 12))
    ax = [fig.add_subplot(3, 3, i) for i in range(1, 10)]

    pressures = np.linspace(1.0e7, 30.0e9, 101)
    G_iso = np.empty_like(pressures)
    G_aniso = np.empty_like(pressures)
    C = np.empty((len(pressures), 6, 6))
    SN = np.empty((len(pressures), 6, 6))
    betaN = np.empty_like(pressures)
    ST = np.empty((len(pressures), 6, 6))
    betaT = np.empty_like(pressures)

    f = np.empty_like(pressures)
    dXdf = np.empty_like(pressures)

    i_pq = ((1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3), (2, 3))

    m.set_state(1.0e5, 300.0)
    rho0 = m.rho
    temperatures = [300.0, 500.0, 750.0, 900.0]
    for T in temperatures:
        for i, P in enumerate(pressures):
            m.set_state(P, T)
            C[i] = m.isentropic_stiffness_tensor
            SN[i] = m.isentropic_compliance_tensor
            betaN[i] = 1.0 / m.isentropic_bulk_modulus_reuss
            ST[i] = m.isothermal_compliance_tensor
            betaT[i] = 1.0 / m.isothermal_bulk_modulus_reuss
            f[i] = np.log(rho0 / m.rho)
        # TK, PGPa, rho, rhoerr = d[:4]
        # C11, C11err = d[4:6]
        # C22, C22err = d[6:8]
        # C33, C33err = d[8:10]
        # C44, C44err = d[10:12]
        # C55, C55err = d[12:14]
        # C66, C66err = d[14:16]
        # C12, C12err = d[16:18]
        # C13, C13err = d[18:20]
        # C23, C23err = d[20:22]
        T_data = np.array([d for d in ol_data if np.abs(d[0] - T) < 1])

        data_mask = np.array([i for i, d in enumerate(ol_data) if np.abs(d[0] - T) < 1])

        if False:
            for i, (p, q) in enumerate(i_pq):
                ln = ax[i].plot(f, SN[:, p - 1, q - 1] / betaN, label=f"{T} K")
                # ln = ax[i].plot(f, ST[:, p-1, q-1]/betaT, linestyle=':',
                # color=ln[0].get_color())
                ax[i].scatter(
                    np.log(rho0 / 1000.0 / ol_data[data_mask, 2]),
                    SNoverbetaN_data[data_mask, p - 1, q - 1],
                    color=ln[0].get_color(),
                )
        else:
            for i, (p, q) in enumerate(i_pq):
                ln = ax[i].plot(
                    pressures / 1.0e9, C[:, p - 1, q - 1] / 1.0e9, label=f"{T} K"
                )
                j = 4 + 2 * i
                ax[i].scatter(T_data[:, 1], T_data[:, j], color=ln[0].get_color())
                ax[i].errorbar(
                    T_data[:, 1],
                    T_data[:, j],
                    yerr=T_data[:, j + 1],
                    linestyle="None",
                    color=ln[0].get_color(),
                )

    for i, (p, q) in enumerate(i_pq):
        ax[i].set_xlabel("Pressure (GPa)")
        ax[i].set_ylabel(f"$C_{{N {p}{q}}}$ (GPa)")
        ax[i].legend()

    fig.set_layout_engine("tight")
    fig.savefig("olivine_CNijs.pdf")
    plt.show()

    fig = plt.figure(figsize=(12, 7))
    ax = [fig.add_subplot(2, 3, i, projection="polar") for i in range(1, 7)]

    P = 3.0e9
    T = 1600.0
    m.set_state(P, T)
    plot_types = [
        "vp",
        "vs1",
        "vp/vs1",
        "s anisotropy",
        "linear compressibility",
        "youngs modulus",
    ]

    contour_sets, ticks, lines = plot_projected_elastic_properties(m, plot_types, ax)
    for i in range(len(contour_sets)):
        cbar = fig.colorbar(contour_sets[i], ax=ax[i], ticks=ticks[i], pad=0.1)
        cbar.add_lines(lines[i])

    fig.set_layout_engine("tight")
    fig.savefig(f"olivine_seismic_properties_{P/1.e9:.2f}_GPa_{int(T)}_K.pdf")
    plt.show()
