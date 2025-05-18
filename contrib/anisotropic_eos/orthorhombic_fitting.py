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

from tools import print_table_for_mineral_constants
from burnman.tools.plot import plot_projected_elastic_properties


run_fitting = False

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
CN_data = CN_data.T
SN_data = np.linalg.inv(CN_data)
betaN_data = np.sum(SN_data[:, :3, :3], axis=(1, 2))
SNoverbetaN_data = (SN_data.T / betaN_data).T


def make_orthorhombic_mineral_from_parameters(x):
    f_order = 3
    Pth_order = 2
    constants = np.zeros((6, 6, f_order + 1, Pth_order + 1))

    san_carlos_params = {
        "name": "San Carlos olivine",
        "formula": formula,
        "equation_of_state": "slb3",
        "F_0": 0.0,
        "V_0": V_0_guess,  # we overwrite this in a second
        "K_0": 1.263e11,  # Abramson et al. 1997
        "Kprime_0": 4.28,  # Abramson et al. 1997
        "Debye_0": fo.params["Debye_0"] * 0.9 + fa.params["Debye_0"] * 0.1,  #
        "grueneisen_0": 0.99282,  # Fo in SLB2011
        "q_0": 2.10672,  # Fo in SLB2011
        "G_0": 81.6e9,
        "Gprime_0": 1.46257,
        "eta_s_0": 2.29972,
        "n": 7.0,
        "molar_mass": formula_mass,
    }

    san_carlos_property_modifiers = [
        [
            "linear",
            {
                "delta_E": 0.0,
                "delta_S": 26.76 * 0.1
                - 2.0
                * burnman.constants.gas_constant
                * (0.1 * np.log(0.1) + 0.9 * np.log(0.9)),
                "delta_V": 0.0,
            },
        ]
    ]

    ol = burnman.Mineral(
        params=san_carlos_params, property_modifiers=san_carlos_property_modifiers
    )

    # Overwrite some properties
    ol.params["V_0"] = x[0] * V_0_guess  # Abramson et al. 1997
    ol.params["K_0"] = x[1] * 1.263e11  # Abramson et al. 1997
    ol.params["Kprime_0"] = x[2] * 4.28  # Abramson et al. 1997
    # ol.params['Debye_0'] = x[3]*809.1703 # Fo in SLB2011 strong tendency to 0
    ol.params["grueneisen_0"] = x[3] * 0.99282  # Fo in SLB2011
    ol.params["q_0"] = x[4] * 2.10672  # Fo in SLB2011

    # Next, each of the eight independent elastic tensor component get their turn.
    # We arbitrarily choose S[2,3] as the ninth component, which is determined by the others.
    i = 5
    for p, q in ((1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3)):
        for m, n in ((1, 0), (2, 0), (3, 0)):
            constants[p - 1, q - 1, m, n] = x[i]
            constants[q - 1, p - 1, m, n] = x[i]
            i += 1

        for m, n in ((0, 1), (1, 1), (2, 1), (3, 1)):
            constants[p - 1, q - 1, m, n] = x[i] * 1.0e-11
            constants[q - 1, p - 1, m, n] = x[i] * 1.0e-11
            i += 1

        for m, n in ((0, 2),):
            constants[p - 1, q - 1, m, n] = x[i] * 1.0e-22
            constants[q - 1, p - 1, m, n] = x[i] * 1.0e-22
            i += 1

    assert i == 69  # 40 parameters

    # Fill the values for the dependent element c[2,3]
    constants[1, 2, 1, 0] = (1.0 - np.sum(constants[:3, :3, 1, 0])) / 2.0
    constants[1, 2, 2:, 0] = -np.sum(constants[:3, :3, 2:, 0], axis=(0, 1)) / 2.0
    constants[1, 2, :, 1:] = -np.sum(constants[:3, :3, :, 1:], axis=(0, 1)) / 2.0

    # And for c[3,2]
    constants[2, 1, :, :] = constants[1, 2, :, :]

    cell_lengths = cell_lengths_0_guess * np.cbrt(ol.params["V_0"] / V_0_guess)
    ol_cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], 90, 90, 90]
    )

    m = AnisotropicMineral(ol, ol_cell_parameters, constants)
    return m


sol = []
if run_fitting:

    def orthorhombic_misfit(x, imin):
        m = make_orthorhombic_mineral_from_parameters(x)

        chisqr = 0.0
        try:
            for d in ol_data:
                TK, PGPa, rho, rhoerr = d[:4]
                C11, C11err = d[4:6]
                C22, C22err = d[6:8]
                C33, C33err = d[8:10]
                C44, C44err = d[10:12]
                C55, C55err = d[12:14]
                C66, C66err = d[14:16]
                C12, C12err = d[16:18]
                C13, C13err = d[18:20]
                C23, C23err = d[20:22]

                PPa = PGPa * 1.0e9

                m.set_state(PPa, TK)

                CN = m.isentropic_stiffness_tensor / 1.0e9

                chisqr += np.power((m.density / 1000.0 - rho) / rhoerr, 2.0)
                chisqr += np.power((CN[0, 0] - C11) / C11err, 2.0)
                chisqr += np.power((CN[1, 1] - C22) / C22err, 2.0)
                chisqr += np.power((CN[2, 2] - C33) / C33err, 2.0)
                chisqr += np.power((CN[3, 3] - C44) / C44err, 2.0)
                chisqr += np.power((CN[4, 4] - C55) / C55err, 2.0)
                chisqr += np.power((CN[5, 5] - C66) / C66err, 2.0)
                chisqr += np.power((CN[0, 1] - C12) / C12err, 2.0)
                chisqr += np.power((CN[0, 2] - C13) / C13err, 2.0)
                chisqr += np.power((CN[1, 2] - C23) / C23err, 2.0)

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

        except:
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
            1.00258391,
            0.99370699,
            0.99752955,
            1.12714466,
            0.37112938,
            0.44371596,
            -0.91751765,
            0.57212708,
            0.26245383,
            0.33363313,
            -5.73660271,
            13.61011798,
            0.50361811,
            0.77519018,
            -0.99654079,
            5.48884001,
            0.51512727,
            1.48952566,
            13.17545958,
            -8.79522221,
            0.31525521,
            0.65662634,
            -0.6299461,
            7.97883382,
            1.71310531,
            1.59569965,
            6.1824309,
            3.91279783,
            0.45100248,
            1.973719,
            -0.53564915,
            21.27700395,
            -3.15830421,
            2.62856492,
            -5.2336551,
            12.4248871,
            -6.05264045,
            1.61904353,
            -1.93685182,
            23.30968136,
            -2.23850942,
            2.7957687,
            1.79408332,
            3.61628653,
            -6.49969595,
            1.55704258,
            -1.38247656,
            12.99413311,
            -0.82259761,
            5.79739211,
            7.82201844,
            -4.90148767,
            -2.3374101,
            -0.12059504,
            0.47144612,
            -1.54276206,
            0.46036409,
            -0.19727193,
            8.82188458,
            8.52558995,
            -0.02754936,
            -0.10031719,
            0.25444599,
            -1.90038101,
            -0.69354727,
            -0.38390188,
            -4.39474523,
            12.20397732,
            0.06899027,
        ]
    )

    i = 0
    min = 1.0e10
    sol = minimize(
        orthorhombic_misfit,
        guesses,
        method="COBYLA",
        args=[[i, min]],
        options={"rhobeg": 0.2, "maxiter": 10000},
    )
    print(sol)

do_plotting = True
if do_plotting:
    if run_fitting:
        m = make_orthorhombic_mineral_from_parameters(sol.x)
    else:
        m = make_orthorhombic_mineral_from_parameters(
            [
                1.00258391,
                0.99370699,
                0.99752955,
                1.12714466,
                0.37112938,
                0.44371596,
                -0.91751765,
                0.57212708,
                0.26245383,
                0.33363313,
                -5.73660271,
                13.61011798,
                0.50361811,
                0.77519018,
                -0.99654079,
                5.48884001,
                0.51512727,
                1.48952566,
                13.17545958,
                -8.79522221,
                0.31525521,
                0.65662634,
                -0.6299461,
                7.97883382,
                1.71310531,
                1.59569965,
                6.1824309,
                3.91279783,
                0.45100248,
                1.973719,
                -0.53564915,
                21.27700395,
                -3.15830421,
                2.62856492,
                -5.2336551,
                12.4248871,
                -6.05264045,
                1.61904353,
                -1.93685182,
                23.30968136,
                -2.23850942,
                2.7957687,
                1.79408332,
                3.61628653,
                -6.49969595,
                1.55704258,
                -1.38247656,
                12.99413311,
                -0.82259761,
                5.79739211,
                7.82201844,
                -4.90148767,
                -2.3374101,
                -0.12059504,
                0.47144612,
                -1.54276206,
                0.46036409,
                -0.19727193,
                8.82188458,
                8.52558995,
                -0.02754936,
                -0.10031719,
                0.25444599,
                -1.90038101,
                -0.69354727,
                -0.38390188,
                -4.39474523,
                12.20397732,
                0.06899027,
            ]
        )

    assert burnman.tools.eos.check_anisotropic_eos_consistency(m)

    print(
        "The following parameters were used for the volumetric part of "
        f'the isotropic model: $V_0$: {m.params["V_0"]*1.e6:.5f} cm$^3$/mol, '
        f'$K_0$: {m.params["K_0"]/1.e9:.5f} GPa, '
        f'$K\'_0$: {m.params["Kprime_0"]:.5f}, '
        f'$\Theta_0$: {m.params["Debye_0"]:.5f} K, '
        f'$\gamma_0$: {m.params["grueneisen_0"]:.5f}, '
        f'and $q_0$: {m.params["q_0"]:.5f}.'
    )

    print_table_for_mineral_constants(
        m, [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3), (2, 3)]
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
            ax[0].plot(
                temperatures,
                ol_SLB.evaluate(["alpha"], pressures, temperatures)[0] * 1.0e5 / 3.0,
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

        for i, (p, q) in enumerate(i_pq):
            ln = ax[i].plot(f, SN[:, p - 1, q - 1] / betaN, label=f"{T} K")
            # ln = ax[i].plot(f, ST[:, p-1, q-1]/betaT, linestyle=':', color=ln[0].get_color())
            ax[i].scatter(
                np.log(rho0 / 1000.0 / ol_data[data_mask, 2]),
                SNoverbetaN_data[data_mask, p - 1, q - 1],
                color=ln[0].get_color(),
            )
        """
        for i, (p, q) in enumerate(i_pq):
            ln = ax[i].plot(pressures/1.e9, C[:, p-1, q-1]/1.e9, label=f'{T} K')
            j = 4 + 2*i
            ax[i].scatter(T_data[:,1], T_data[:,j], color=ln[0].get_color())
            ax[i].errorbar(T_data[:,1], T_data[:,j], yerr=T_data[:,j+1],
                           linestyle='None', color=ln[0].get_color())
        """

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
