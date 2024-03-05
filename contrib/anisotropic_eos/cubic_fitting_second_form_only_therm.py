# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
cubic_fitting
-------------

This script creates an AnisotropicMineral object corresponding to
periclase (a cubic mineral). If run_fitting is set to True, the script uses
experimental data to find the optimal anisotropic parameters.
If set to False, it uses pre-optimized parameters.
The data is used only to optimize the anisotropic parameters;
the isotropic equation of state is taken from
Stixrude and Lithgow-Bertelloni (2011).

The script ends by making three plots; one with elastic moduli
at high pressure, one with the corresponding shear moduli,
and one with the elastic moduli at 1 bar.
"""

from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, differential_evolution
from tools import print_table_for_mineral_constants

import burnman

from burnman import AnisotropicMineral


run_fitting = False

per = burnman.minerals.SLB_2011.periclase()
V0 = per.params["V_0"]
K0 = per.params["K_0"]
Kprime0 = per.params["Kprime_0"]
gr0 = per.params["grueneisen_0"]
q0 = per.params["q_0"]


atherm_params = {}
atherm_params["a"] = np.array(
    [
        [0.47944268, -0.07305467, -0.07305467, 0.0, 0.0, 0.0],
        [-0.07305467, 0.47944268, -0.07305467, 0.0, 0.0, 0.0],
        [-0.07305467, -0.07305467, 0.47944268, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.04445275, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.04445275, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.04445275],
    ]
)
atherm_params["b_1"] = np.array(
    [
        [-0.77294538, 0.38647269, 0.38647269, 0.0, 0.0, 0.0],
        [0.38647269, -0.77294538, 0.38647269, 0.0, 0.0, 0.0],
        [0.38647269, 0.38647269, -0.77294538, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -0.17422384, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -0.17422384, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, -0.17422384],
    ]
)
atherm_params["c_1"] = np.array(
    [
        [5.81473511, 5.81473511, 5.81473511, 1.0, 1.0, 1.0],
        [5.81473511, 5.81473511, 5.81473511, 1.0, 1.0, 1.0],
        [5.81473511, 5.81473511, 5.81473511, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.13934342, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.13934342, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.13934342],
    ]
)
atherm_params["b_2"] = np.array(
    [
        [0.79772947, -0.39886474, -0.39886474, 0.0, 0.0, 0.0],
        [-0.39886474, 0.79772947, -0.39886474, 0.0, 0.0, 0.0],
        [-0.39886474, -0.39886474, 0.79772947, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -0.53002167, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, -0.53002167, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, -0.53002167],
    ]
)
atherm_params["c_2"] = np.array(
    [
        [5.81567299, 5.81567299, 5.81567299, 1.0, 1.0, 1.0],
        [5.81567299, 5.81567299, 5.81567299, 1.0, 1.0, 1.0],
        [5.81567299, 5.81567299, 5.81567299, 1.0, 1.0, 1.0],
        [1.0, 1.0, 1.0, -2.37105322, 1.0, 1.0],
        [1.0, 1.0, 1.0, 1.0, -2.37105322, 1.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, -2.37105322],
    ]
)


per.set_state(1.0e5, 300.0)
beta_RT = per.beta_T


def psi_func(f, Pth, params):
    dPsidf = (
        params["a"]
        + params["b_1"] * params["c_1"] * np.exp(params["c_1"] * f)
        + params["b_2"] * params["c_2"] * np.exp(params["c_2"] * f)
    )
    dPsidf = dPsidf + Pth / 1.0e9 * (
        params["b_3"] * params["c_3"] * np.exp(params["c_3"] * f)
        + params["b_4"] * params["c_4"] * np.exp(params["c_4"] * f)
    )

    Psi = (
        0.0
        + params["a"] * f
        + params["b_1"] * (np.exp(params["c_1"] * f) - 1.0)
        + params["b_2"] * (np.exp(params["c_2"] * f) - 1.0)
    )
    Psi = Psi + Pth / 1.0e9 * (
        params["b_3"] * np.exp(params["c_3"] * f)
        + params["b_4"] * np.exp(params["c_4"] * f)
    )

    dPsidPth = (
        params["b_3"] * np.exp(params["c_3"] * f)
        + params["b_4"] * np.exp(params["c_4"] * f)
    ) / 1.0e9
    return (Psi, dPsidf, dPsidPth)


def make_cubic_mineral_from_parameters(x):
    anisotropic_parameters = {
        "b_3": np.zeros((6, 6)),
        "c_3": np.ones((6, 6)),
        "b_4": np.zeros((6, 6)),
        "c_4": np.ones((6, 6)),
        **atherm_params,
    }

    per.params["V_0"] = x[0] * V0
    per.params["K_0"] = x[1] * K0
    per.params["Kprime_0"] = x[2] * Kprime0
    per.params["grueneisen_0"] = x[3] * gr0
    per.params["q_0"] = np.power(x[4], 2.0) * q0

    a = np.cbrt(per.params["V_0"])
    cell_parameters = np.array([a, a, a, 90, 90, 90])

    i = 5
    for p in ["b_3", "c_3", "b_4", "c_4"]:
        if p == "a":
            anisotropic_parameters[p][:3, :3] = (1.0 - 3.0 * x[i]) / 6.0
        elif p[0] == "b" or p == "d":
            anisotropic_parameters[p][:3, :3] = -3.0 * x[i] / 6.0
        else:
            anisotropic_parameters[p][:3, :3] = x[i]

        anisotropic_parameters[p][0, 0] = x[i]
        anisotropic_parameters[p][1, 1] = x[i]
        anisotropic_parameters[p][2, 2] = x[i]
        i = i + 1
        anisotropic_parameters[p][3, 3] = x[i]
        anisotropic_parameters[p][4, 4] = x[i]
        anisotropic_parameters[p][5, 5] = x[i]
        i = i + 1
    assert len(x) == i

    return AnisotropicMineral(
        per, cell_parameters, anisotropic_parameters, psi_func, orthotropic=True
    )


per_data = np.loadtxt("data/Karki_2000_periclase_CSijs.dat")
per_data_2 = np.loadtxt("data/isentropic_stiffness_tensor_periclase.dat")
parameters = []

if run_fitting:

    def cubic_misfit(x, imin):
        m = make_cubic_mineral_from_parameters(x)

        chisqr = 0.0
        try:
            for d in per_data[::1]:
                T, P, C11S, C12S, C44S, KT, V, betaTminusbetaS = d

                C11Serr = 0.01 * C11S
                C12Serr = 0.01 * C12S
                C44Serr = 0.01 * C44S
                Verr = 0.01 * V

                if T < 3500.0:
                    m.set_state(P, T)

                    chisqr += np.power((m.V - V) / Verr, 2.0)
                    chisqr += np.power(
                        (m.isentropic_stiffness_tensor[0, 0] - C11S) / C11Serr, 2.0
                    )
                    chisqr += np.power(
                        (m.isentropic_stiffness_tensor[0, 1] - C12S) / C12Serr, 2.0
                    )
                    chisqr += np.power(
                        (m.isentropic_stiffness_tensor[3, 3] - C44S) / C44Serr, 2.0
                    )

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

    i = 0
    min = 1.0e10
    guesses = np.array(
        [
            1.00827773e00,
            9.98856164e-01,
            1.08401459e00,
            1.10777917e00,
            7.94493025e-01,
            2.04771439e-03,
            8.76321298e-03,
            1.87121476e00,
            2.01075106e00,
            1.68704173e-04,
            -2.38606584e-03,
            2.39190918e01,
            5.33127227e00,
        ]
    )

    # sol = minimize(cubic_misfit, guesses, method='COBYLA', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 10000})
    sol = minimize(
        cubic_misfit,
        guesses,
        method="Powell",
        args=[[i, min]],
        options={"rhobeg": 0.2, "maxiter": 10000},
    )  # good for rapid convergence
    # sol = minimize(cubic_misfit, guesses, method='Nelder-Mead', args=[[i, min]], options={'rhobeg': 0.2, 'maxiter': 100000}) # fine tuning
    print(sol)
    parameters = sol.x

if not run_fitting:
    parameters = np.array(
        [
            1.00827773e00,
            9.98856164e-01,
            1.08401459e00,
            1.10777917e00,
            7.94493025e-01,
            2.04771439e-03,
            8.76321298e-03,
            1.87121476e00,
            2.01075106e00,
            1.68704173e-04,
            -2.38606584e-03,
            2.39190918e01,
            5.33127227e00,
        ]
    )

m = make_cubic_mineral_from_parameters(parameters)

assert burnman.tools.eos.check_anisotropic_eos_consistency(m)

m.set_state(1.0e5, 300.0)

fig = plt.figure(figsize=(4, 12))
fig2 = plt.figure(figsize=(4, 4))
ax = [fig.add_subplot(3, 1, i) for i in range(1, 4)]
ax.append(fig2.add_subplot(1, 1, 1))


temperatures = [300.0, 1000.0, 2000.0, 3000.0]
for T in temperatures:
    # T, P, C11S, C12S, C44S, KT, V, betaTminusbetaS
    # every 3rd data point
    T_data = (
        np.array([[d[1], d[2], d[3], d[4]] for d in per_data if np.abs(d[0] - T) < 1])
        / 1.0e9
    )[::3]

    Pmin = np.min(T_data[:, 0]) * 1.0e9
    Pmax = np.max(T_data[:, 0]) * 1.0e9
    pressures = np.linspace(Pmin, Pmax, 101)
    G_iso = np.empty_like(pressures)
    G_aniso = np.empty_like(pressures)
    C11 = np.empty_like(pressures)
    C12 = np.empty_like(pressures)
    C44 = np.empty_like(pressures)
    G_slb = np.empty_like(pressures)
    G_V = np.empty_like(pressures)
    G_R = np.empty_like(pressures)

    f = np.empty_like(pressures)
    dXdf = np.empty_like(pressures)

    for i, P in enumerate(pressures):
        per.set_state(P, T)
        m.set_state(P, T)
        C11[i] = m.isentropic_stiffness_tensor[0, 0]
        C12[i] = m.isentropic_stiffness_tensor[0, 1]
        C44[i] = m.isentropic_stiffness_tensor[3, 3]
        G_V[i] = m.isentropic_shear_modulus_voigt
        G_R[i] = m.isentropic_shear_modulus_reuss

        G_slb[i] = per.shear_modulus

    ax[0].plot(pressures / 1.0e9, C11 / 1.0e9, label=f"{T} K")
    ax[1].plot(pressures / 1.0e9, C12 / 1.0e9, label=f"{T} K")
    ax[2].plot(pressures / 1.0e9, C44 / 1.0e9, label=f"{T} K")

    ln = ax[3].plot(pressures / 1.0e9, G_R / 1.0e9)
    ax[3].plot(pressures / 1.0e9, G_V / 1.0e9, color=ln[0].get_color(), label=f"{T} K")

    ax[3].fill_between(
        pressures / 1.0e9, G_R / 1.0e9, G_V / 1.0e9, alpha=0.3, color=ln[0].get_color()
    )

    ax[3].plot(
        pressures / 1.0e9,
        G_slb / 1.0e9,
        label=f"{T} K (SLB2011)",
        linestyle="--",
        color=ln[0].get_color(),
    )

    ax[0].scatter(T_data[:, 0], T_data[:, 1])
    ax[1].scatter(T_data[:, 0], T_data[:, 2])
    ax[2].scatter(T_data[:, 0], T_data[:, 3])

for i in range(4):
    ax[i].set_xlabel("Pressure (GPa)")

ax[0].set_ylabel("$C_{N 11}$ (GPa)")
ax[1].set_ylabel("$C_{N 12}$ (GPa)")
ax[2].set_ylabel("$C_{N 44}$ (GPa)")
ax[3].set_ylabel("$G$ (GPa)")

for i in range(4):
    ax[i].legend()

fig.set_tight_layout(True)
fig.savefig("periclase_stiffness_tensor.pdf")
fig2.set_tight_layout(True)
fig2.savefig("periclase_shear_modulus.pdf")
plt.show()

fig = plt.figure(figsize=(4, 4))
ax = [fig.add_subplot(1, 1, 1)]

temperatures = np.linspace(10.0, 2000.0, 101)
P = 1.0e5
for i, T in enumerate(temperatures):
    m.set_state(P, T)
    per.set_state(P, T)
    C11[i] = m.isentropic_stiffness_tensor[0, 0]
    C12[i] = m.isentropic_stiffness_tensor[0, 1]
    C44[i] = m.isentropic_stiffness_tensor[3, 3]
    G_V[i] = m.isentropic_shear_modulus_voigt
    G_R[i] = m.isentropic_shear_modulus_reuss

    G_slb[i] = per.shear_modulus
ax[0].plot(temperatures, C11 / 1.0e9, label="$C_{N 11}$")
ax[0].plot(temperatures, C12 / 1.0e9, label="$C_{N 12}$")
ax[0].plot(temperatures, C44 / 1.0e9, label="$C_{44}$")
ln = ax[0].plot(temperatures, G_R / 1.0e9, color=ln[0].get_color(), label=f"$G$")
ax[0].plot(temperatures, G_V / 1.0e9, color=ln[0].get_color())
ax[0].fill_between(
    temperatures, G_R / 1.0e9, G_V / 1.0e9, alpha=0.3, color=ln[0].get_color(), zorder=2
)
ax[0].plot(
    temperatures,
    G_slb / 1.0e9,
    label="$G$ (SLB2011)",
    color=ln[0].get_color(),
    linestyle="--",
    zorder=1,
)

# T, PGPa, Perr, rho, rhoerr, C11S, C11Serr, C12S, C12Serr, C44S, C44Serr
LP_data = np.array([[d[0], d[5], d[7], d[9]] for d in per_data_2 if d[1] < 0.1])

ax[0].scatter(LP_data[:, 0], LP_data[:, 1])
ax[0].scatter(LP_data[:, 0], LP_data[:, 2])
ax[0].scatter(LP_data[:, 0], LP_data[:, 3])

ax[0].set_xlabel("Temperature (K)")
ax[0].set_ylabel("Elastic modulus (GPa)")

# print_table_for_mineral_constants(m, [(1, 1), (1, 2), (4, 4)])

ax[0].legend()

fig.set_tight_layout(True)
fig.savefig("periclase_properties_1bar.pdf")
plt.show()
