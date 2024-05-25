from __future__ import absolute_import

from burnman import Composite, equilibrate
from burnman.constants import gas_constant
from burnman.tools.polytope import simplify_composite_with_composition
from burnman.tools.eos import check_eos_consistency
from burnman.tools.unitcell import molar_volume_from_unit_cell_volume
from burnman.minerals import SLB_2024
from burnman.minerals.HP_2011_fluids import O2
from burnman.eos.property_modifiers import calculate_property_modifications
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


def check_iron_bearing_mineral_property_values():
    print("Checking iron bearing mineral properties...")
    m_names = np.genfromtxt("slb_2024_Fe_endmembers.dat", dtype=str)[:, 1]
    minerals = [eval("SLB_2024." + m)() for m in m_names]
    d = np.genfromtxt("slb_2024_Fe_endmembers.dat", dtype=float)[:, 2:]
    for i, m in enumerate(minerals):
        P, T, _, entropy, Cp, gibbs, volume, _, _, _, _ = d[i]
        m.set_state(P * 1.0e4, T)
        delta_gibbs = gibbs * 1.0e3 - m.gibbs
        if np.abs(delta_gibbs) > 1.0:
            print(m_names[i])
            print("delta Gibbs:", delta_gibbs)
            print("delta S:", entropy - m.S)
            print("delta Cp:", Cp - m.C_p)
            print("delta V:", volume * 1.0e-6 - m.V)
        else:
            rel_error = max(
                [
                    (entropy - m.S) / entropy,
                    (Cp - m.C_p) / Cp,
                    (volume - m.V * 1.0e6) / volume,
                ]
            )
            print(
                f"{m_names[i]} in agreement with HeFESTo (max rel error {rel_error*100.:.1g} %)"
            )


def check_bcc_iron_consistency():
    print("")
    fe_bcc = SLB_2024.alpha_bcc_iron()
    assert check_eos_consistency(fe_bcc, verbose=True)


def check_fper_entropy():
    print("\nChecking fper configurational entropy...")
    fper = SLB_2024.ferropericlase()
    print("Endmember entropies (should all be zero):")
    print(np.abs(fper.solution_model.endmember_configurational_entropies))

    fper.set_composition([0.0, 0.5, 0.0, 0.0, 0.5])
    Sxs = -2 * gas_constant * np.log(0.5)
    print(
        f"Equimolar wu-mag: {fper.excess_entropy:.5f} J/K/mol (should be {Sxs:.5f} J/K/mol)"
    )


def check_fig_1_fper_relaxed():
    print("\nChecking Figure 1...")
    fper = SLB_2024.ferropericlase_relaxed()
    c = molar_volume_from_unit_cell_volume(1.0, SLB_2024.periclase().params["Z"])

    pressures = np.linspace(1.0e5, 140.0e9, 141)
    temperatures = pressures * 0.0 + 300.0

    fig = plt.figure(figsize=(12, 8))
    fig.suptitle("Figure 1 (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(2, 2, i) for i in range(1, 5)]

    fig1 = mpimg.imread("figures/SLB_2024_fper_V.png")
    fig2 = mpimg.imread("figures/SLB_2024_fper_K_T.png")
    fig3 = mpimg.imread("figures/SLB_2024_fper_spin_20_50_80.png")
    fig4 = mpimg.imread("figures/SLB_2024_fper_spin_percents.png")

    ax[0].imshow(fig1, extent=[0.0, 140.0, 50.0, 80.0], aspect="auto")
    ax[1].imshow(fig2, extent=[0.0, 140.0, 100.0, 700.0], aspect="auto")
    ax[2].imshow(fig3, extent=[0.0, 1.0, 0.0, 100.0], aspect="auto")
    ax[3].imshow(fig4, extent=[0.0, 200.0, 0.0, 4000.0], aspect="auto")

    for x_fe in [0.1, 0.3, 0.5, 1.0]:
        fper.set_composition([1.0 - x_fe, x_fe, 0.0, 0.0])
        Vs, K_Ss = fper.evaluate(["V", "K_S"], pressures, temperatures)

        ax[0].plot(pressures / 1.0e9, Vs / c, linestyle=":", linewidth=3.0)
        ax[1].plot(pressures / 1.0e9, K_Ss / 1.0e9, linestyle=":", linewidth=3.0)

    ax[0].set_ylim(50.0, 80.0)
    ax[1].set_ylim(100.0, 700.0)

    fper.set_composition([0.5, 0.5, 0.0, 0.0])
    fper.set_state(70.0e9, 300.0)
    print("Mg2Fe2O4 ferropericlase at 70 GPa, 300 K:")
    fstr = ", ".join([f"{f:.2f}" for f in fper.molar_fractions])
    print(f"Molar proportions: {fstr}")
    print(f"Volume: {fper.V*1.e6:.2f} cm^3/mol")

    fper = SLB_2024.ferropericlase()
    assemblage = simplify_composite_with_composition(
        Composite([fper]), {"Mg": 0.5, "Fe": 0.5, "O": 1.0}
    )
    assemblage.set_state(50.0e9, 300.0)
    fper = assemblage.phases[0]
    x_MgOs = np.linspace(1.0e-5, 0.9999, 101)
    pressures = np.empty((3, 101))
    for i, fLS in enumerate([0.2, 0.5, 0.8]):
        for j, x_MgO in enumerate(x_MgOs):
            composition = {"Mg": x_MgO, "Fe": 1.0 - x_MgO, "O": 1.0}
            fper.set_composition([x_MgO, (1.0 - x_MgO) / 2.0, (1.0 - x_MgO) / 2.0])
            equality_constraints = [
                ["T", 300.0],
                [
                    "phase_composition",
                    (
                        fper,
                        (
                            ["Fe_A", "Fels_A", "Fe_B", "Fels_B"],
                            [0.0, 1.0, 0.0, 1.0],
                            [1.0, 1.0, 1.0, 1.0],
                            fLS,
                        ),
                    ),
                ],
            ]
            equilibrate(composition, assemblage, equality_constraints)
            pressures[i, j] = assemblage.pressure
    ax[2].plot(
        1.0 - x_MgOs,
        (pressures[2] - pressures[0]) / 1.0e9,
        linestyle=":",
        linewidth=3.0,
    )
    ax[2].plot(1.0 - x_MgOs, (pressures[1]) / 1.0e9, linestyle=":", linewidth=3.0)

    composition = {"Mg": 0.75, "Fe": 0.25, "O": 1.0}
    fractions = np.linspace(0.1, 0.9, 9)
    temperatures = np.linspace(1.0, 4000.0, 21)
    equality_constraints = [
        ["T", temperatures],
        [
            "phase_composition",
            (
                fper,
                (
                    ["Fe_A", "Fels_A", "Fe_B", "Fels_B"],
                    [0.0, 1.0, 0.0, 1.0],
                    [1.0, 1.0, 1.0, 1.0],
                    fractions,
                ),
            ),
        ],
    ]
    solss = equilibrate(composition, assemblage, equality_constraints)[0]
    pressures = np.array([[sol.assemblage.pressure for sol in sols] for sols in solss])
    for i, f in enumerate(fractions):
        ax[3].plot(pressures[:, i] / 1.0e9, temperatures, linestyle=":", linewidth=3.0)
    ax[3].set_xlim(0.0, 200.0)
    ax[3].set_ylim(0.0, 4000.0)

    plt.show()


def check_fig_3_fcc_ferric_fper():
    print("\nChecking Figure 3...")

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle("Figure 3 (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    fig1 = mpimg.imread("figures/SLB_2024_Fe_O_phase_diagram.png")
    fig2 = mpimg.imread("figures/SLB_2024_fper_fO2.png")

    ax[0].imshow(fig1, extent=[-1, -0.75, 700.0, 1800.0], aspect="auto")
    ax[1].imshow(fig2, extent=[0.0, 1.0, 0.0, 0.12], aspect="auto")

    # Subplot a
    bcc = SLB_2024.alpha_bcc_iron()
    fcc = SLB_2024.gamma_fcc_iron()
    fper = SLB_2024.ferropericlase()
    mag = SLB_2024.smag()

    assemblage = simplify_composite_with_composition(
        Composite([fper]), {"Fe": 1.0, "O": 1.1}
    )
    fper = assemblage.phases[0]

    composition = {"Fe": 1.0, "O": 1.0}
    assemblage = Composite([bcc, fper, mag], [0.5, 0.49, 0.01])
    fper.set_composition([0.05, 0.0, 0.95])
    equality_constraints = [["P", 1.0e5], ["phase_fraction", (mag, 0.0)]]
    sol = equilibrate(composition, assemblage, equality_constraints, tol=1.0e-5)
    T_a_wu_mag = sol[0].assemblage.temperature
    print(f"BCC-fper-mag triple point at 1 bar: {T_a_wu_mag:.2f} K")

    assemblage = Composite([bcc, fcc], [0.5, 0.5])
    equality_constraints = [["P", 1.0e5], ["phase_fraction", (bcc, 0.0)]]
    sol = equilibrate(composition, assemblage, equality_constraints, tol=1.0e-5)
    T_a_g = sol[0].assemblage.temperature

    assemblage = Composite([bcc, fper], [0.5, 0.5])
    temperatures = np.linspace(T_a_wu_mag, T_a_g, 31)
    equality_constraints = [["P", 1.0e5], ["T", temperatures]]
    sol = equilibrate(composition, assemblage, equality_constraints, tol=1.0e-5)
    xs = np.empty_like(temperatures)
    for i, s in enumerate(sol[0]):
        xs[i] = -s.assemblage.phases[1].formula["Fe"] / 4.0
    ax[0].plot(xs, temperatures, linestyle=":", linewidth=3.0)

    assemblage = Composite([fcc, fper], [0.5, 0.5])
    temperatures = np.linspace(T_a_g, 1800.0, 31)
    equality_constraints = [["P", 1.0e5], ["T", temperatures]]
    sol = equilibrate(composition, assemblage, equality_constraints, tol=1.0e-5)
    xs = np.empty_like(temperatures)
    for i, s in enumerate(sol[0]):
        xs[i] = -s.assemblage.phases[1].formula["Fe"] / 4.0
    ax[0].plot(xs, temperatures, linestyle=":", linewidth=3.0)

    assemblage = Composite([mag, fper], [0.5, 0.5])
    temperatures = np.linspace(T_a_wu_mag, 1800.0, 31)
    equality_constraints = [["P", 1.0e5], ["T", temperatures]]
    sol = equilibrate(composition, assemblage, equality_constraints, tol=1.0e-5)
    xs = np.empty_like(temperatures)
    for i, s in enumerate(sol[0]):
        xs[i] = -s.assemblage.phases[1].formula["Fe"] / 4.0
    ax[0].plot(xs, temperatures, linestyle=":", linewidth=3.0)

    # Subplot b
    fcc = SLB_2024.gamma_fcc_iron()
    fper = SLB_2024.ferropericlase()

    x_MgOs = np.linspace(0.01, 0.99, 33)
    x_Mgs = np.empty_like(x_MgOs) * np.nan
    x_Fe3oversumFes = np.empty_like(x_MgOs) * np.nan
    composition = {"Mg": 0.01, "Fe": 2.0, "O": 1.0}
    assemblage = Composite([fper, fcc], [0.5, 0.5])
    x_MgO = 0.01
    fper.set_composition([x_MgO, 0.95 * (1.0 - x_MgO), 0.0, 0.0, 0.05 * (1.0 - x_MgO)])
    assemblage = simplify_composite_with_composition(assemblage, composition)

    for equality_constraints in [
        [["P", 1.0e5], ["T", 1473.0]],
        [["P", 21.0e9], ["T", 1673.0]],
    ]:
        for i, x_MgO in enumerate(x_MgOs):
            composition = {"Mg": x_MgO, "Fe": 2.0, "O": 1.0}
            try:
                sol = equilibrate(
                    composition, assemblage, equality_constraints, tol=1.0e-5
                )
                if sol[0].success:
                    f = assemblage.phases[0].formula
                    x_Mgs[i] = f["Mg"] / 4.0
                    x_Fe3oversumFes[i] = 2.0 * ((f["O"] - f["Mg"]) / f["Fe"] - 1.0)
            except AssertionError:
                pass
        ax[1].plot(x_Mgs, x_Fe3oversumFes, linestyle=":", linewidth=3.0)
    plt.show()


def check_fig_6a_iron_Cp_V():
    print("\nChecking Figure 6a...")
    fe_bcc = SLB_2024.alpha_bcc_iron()
    fe_fcc = SLB_2024.gamma_fcc_iron()

    temperatures = np.linspace(1.0, 2000.0, 101)
    pressures = temperatures * 0.0 + 1.0e5

    fig = plt.figure(figsize=(10, 5))
    fig.suptitle("Figure 6a (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    fig1 = mpimg.imread("figures/SLB_2024_Fe_Cp_V.png")
    ax[0].imshow(fig1, extent=[0.0, 2000.0, -15.0, 70.0], aspect="auto")
    ax[1].imshow(fig1, extent=[0.0, 2000.0, 6.9, 8.4], aspect="auto")

    for m in [fe_bcc, fe_fcc]:
        Cp, V = m.evaluate(["C_p", "V"], pressures, temperatures)
        ax[0].plot(temperatures, Cp, linestyle=":", linewidth=3.0)
        ax[1].plot(temperatures, V * 1.0e6, linestyle=":", linewidth=3.0)

    fe_bcc.set_state(1.0e5, 1000.0)
    print(
        f"BCC iron Cp, V at 1 bar, 1000 K: {fe_bcc.C_p:.2f} J/K/mol, {fe_bcc.V*1.e6:.2f} cm^3/mol"
    )
    ax[0].set_ylim(0.0, 70.0)
    ax[1].set_ylim(6.9, 7.8)
    plt.show()


def check_fig_6c_iron_phase_diagram():
    print("\nChecking Figure 6c...")
    fe_bcc = SLB_2024.alpha_bcc_iron()
    fe_fcc = SLB_2024.gamma_fcc_iron()
    fe_hcp = SLB_2024.epsilon_hcp_iron()

    fig = plt.figure(figsize=(8, 5))
    fig.suptitle("Figure 6c (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]
    fig1 = mpimg.imread("figures/SLB_2024_Fe_phase_diagram.png")
    ax[0].imshow(fig1, extent=[0.0, 250.0, 300.0, 5800.0], aspect="auto")

    # Calculate the invariant point
    irons = Composite([fe_bcc, fe_fcc, fe_hcp], [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0])
    irons.set_state(5.0e9, 600.0)
    equilibrate(
        {"Fe": 1.0},
        irons,
        equality_constraints=[
            ["phase_fraction", [fe_bcc, 0.0]],
            ["phase_fraction", [fe_hcp, 0.0]],
        ],
    )
    P_invariant = irons.pressure
    T_invariant = irons.temperature

    print(
        f"The bcc-fcc-hcp invariant is at {P_invariant/1.e9:.2f} GPa, {T_invariant:.1f} K"
    )

    for Tmin, Tmax, m1, m2 in [
        [300.0, T_invariant, fe_bcc, fe_hcp],
        [T_invariant, 2000.0, fe_bcc, fe_fcc],
        [T_invariant, 4800.0, fe_fcc, fe_hcp],
    ]:

        temperatures = np.linspace(Tmin, Tmax, 101)
        irons = Composite([m1, m2], [1.0 / 2.0, 1.0 / 2.0])
        sols = equilibrate(
            {"Fe": 1.0},
            irons,
            equality_constraints=[["phase_fraction", [m1, 0.0]], ["T", temperatures]],
        )[0]
        Ps, Ts = np.array(
            [
                [sol.assemblage.pressure, sol.assemblage.temperature]
                for sol in sols
                if sol.success
            ]
        ).T
        ax[0].plot(Ps / 1.0e9, Ts, linestyle=":", linewidth=3.0)

    plt.show()


def check_fig_7_fO2():
    print("\nChecking Figure 7...")
    hem = SLB_2024.hematite()
    mag = SLB_2024.smag()
    hmag = SLB_2024.high_pressure_magnetit()
    wu = SLB_2024.ferropericlase()
    wu.set_composition([0.0, 0.97, 0.005, 0.0, 0.025])
    wu = simplify_composite_with_composition(
        Composite([wu]), {"Fe": 1, "O": 1.01}
    ).phases[0]

    hepv = SLB_2024.hepv()
    hppv = SLB_2024.hppv()
    Fe_bcc = SLB_2024.alpha_bcc_iron()
    Fe_fcc = SLB_2024.gamma_fcc_iron()
    Fe_hcp = SLB_2024.epsilon_hcp_iron()
    q = SLB_2024.quartz()
    fa = SLB_2024.fayalite()
    O2_gas = O2()

    O2_gas.set_state(1.0e5, 298.15)
    f = O2_gas.S * 298.15

    fig = plt.figure(figsize=(8, 5))
    fig.suptitle("Figure 7 (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]

    fig1 = mpimg.imread("figures/SLB_2024_Fe_O_fO2_T.png")
    ax[0].imshow(fig1, extent=[700.0, 1500.0, -32.0, 0.0], aspect="auto")

    print("fO2 values at 1500 K:")
    for phases in [[fa, Fe_fcc, q], [fa, Fe_bcc, q], [fa, mag, q], [mag, hem, hem]]:
        assemblage = Composite(phases, [0.2, 0.3, 0.5])
        temperatures = np.linspace(700.0, 1500.0)
        logfO2 = np.empty_like(temperatures)

        for i, T in enumerate(temperatures):
            assemblage.set_state(1.0e5, T)
            mu_O2 = assemblage.chemical_potential([{"O": 2.0}])[0]
            O2_gas.set_state(1.0e5, T)
            O2_gibbs = O2_gas.gibbs + f
            logfO2[i] = (mu_O2 - O2_gibbs) / (gas_constant * T) / np.log(10)
        print([ph.name for ph in phases])
        print(f"log10fO2 = {logfO2[-1]:.3f}")
        ax[0].plot(temperatures, logfO2, linestyle=":", linewidth=3.0)

    ax[0].set_xlim(700.0, 1500.0)
    ax[0].set_ylim(-32, 0)
    plt.show()

    P_bounds = [[1.0e5, 100.0e9] for i in range(9)]

    a1 = Composite([mag, hmag], [0.5, 0.5])
    a2 = Composite([hem, hepv], [0.5, 0.5])
    a3 = Composite([hepv, hppv], [0.5, 0.5])
    a4 = Composite([hmag, hppv, wu], [0.3, 0.2, 0.5])
    a5 = Composite([Fe_fcc, Fe_hcp], [0.5, 0.5])

    for i, a in enumerate([a1, a2, a3, a4, "null", a1, "null", a5]):
        if a != "null":
            equality_constraints = [
                ["phase_fraction", [a.phases[0], 0.0]],
                ["T", 1500.0],
            ]
            equilibrate(a.phases[0].formula, a, equality_constraints, tol=1.0e-5)
            P_bounds[i][1] = a.pressure
            P_bounds[i + 1][0] = a.pressure

    P_bounds[6][1] = P_bounds[3][1]

    a1 = Composite([hem, mag], [0.5, 0.5])
    a2 = Composite([hem, hmag], [0.5, 0.5])
    a3 = Composite([hepv, hmag], [0.5, 0.5])
    a4 = Composite([hppv, hmag], [0.5, 0.5])
    a5 = Composite([hppv, wu], [0.5, 0.5])
    a6 = Composite([mag, wu], [0.5, 0.5])
    a7 = Composite([hmag, wu], [0.5, 0.5])
    a8 = Composite([Fe_fcc, wu], [0.5, 0.5])
    a9 = Composite([Fe_hcp, wu], [0.5, 0.5])

    fig = plt.figure(figsize=(8, 5))
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]

    fig1 = mpimg.imread("figures/SLB_2024_Fe_O_fO2.png")
    ax[0].imshow(fig1, extent=[0, 100.0, -14.0, 22.0], aspect="auto")

    for i, assemblage in enumerate([a1, a2, a3, a4, a5, a6, a7, a8, a9]):
        pressures = np.linspace(P_bounds[i][0], P_bounds[i][1], 11)
        logfO2 = np.empty_like(pressures) * np.nan
        T = 1500.0
        for i, P in enumerate(pressures):
            assemblage.set_state(P, T)
            try:
                sol = equilibrate(
                    assemblage.formula, assemblage, [["P", P], ["T", T]], tol=1.0e-7
                )[0]
                if sol.success:
                    mu_O2 = assemblage.chemical_potential([{"O": 2.0}])[0]
                    O2_gas.set_state(1.0e5, T)
                    O2_gibbs = O2_gas.gibbs + f
                    logfO2[i] = (mu_O2 - O2_gibbs) / (gas_constant * T) / np.log(10.0)
            except AssertionError:
                pass
        ax[0].plot(pressures / 1.0e9, logfO2, linestyle=":", linewidth=3.0)

    ax[0].set_ylim(-14, 22)
    plt.show()


def check_fig_a1_fe2o3():
    print("\nChecking Figure A1...")
    hem = SLB_2024.hematite()
    hppv = SLB_2024.hppv()
    bdg = SLB_2024.bridgmanite_relaxed()
    bdg.set_composition([0.0, 0.0, 0.0, 1.0, 0.0, 0.0])

    fig = plt.figure(figsize=(8, 5))
    fig.suptitle("Figure A1 (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]

    fig1 = mpimg.imread("figures/SLB_2024_Fe2O3_V.png")
    ax[0].imshow(fig1, extent=[0.0, 120.0, 19.0, 31.0], aspect="auto")

    pressures = np.linspace(1.0e5, 120.0e9, 101)
    temperatures = 300.0 + 0.0 * pressures
    for phase in [hem, hppv, bdg]:
        volumes = phase.evaluate(["V"], pressures, temperatures)[0]
        ax[0].plot(pressures / 1.0e9, volumes * 1.0e6, linestyle=":", linewidth=3.0)

    bdg.set_composition([0.5, 0.0, 0.0, 0.5, 0.0, 0.0])
    volumes = bdg.evaluate(["V"], pressures, temperatures)[0]
    ax[0].plot(pressures / 1.0e9, volumes * 1.0e6, linestyle=":", linewidth=3.0)

    bdg.set_state(60.0e9, 300.0)
    print(f"(Mg0.5Fe0.5)(Fe0.5Si0.5)O3 V at 60 GPa, 300 K: {bdg.V*1.e6:.3f} cm^3/mol")
    plt.show()


def check_fig_a4_stishovite():
    print("\nChecking Figure A4...")
    stv = SLB_2024.stishovite()

    fig = plt.figure(figsize=(8, 5))
    fig.suptitle("Figure A4 (Stixrude and Lithgow-Bertelloni, 2024)")
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]

    fig1 = mpimg.imread("figures/SLB_2024_stv_VS.png")
    ax[0].imshow(fig1, extent=[0.0, 100.0, 0.0, 9.0], aspect="auto")

    for T in [300, 1000, 2000]:
        pressures = np.linspace(1.0e5, 100.0e9, 101)
        temperatures = pressures * 0.0 + T
        Vs = stv.evaluate(["shear_wave_velocity"], pressures, temperatures)[0]
        ax[0].plot(pressures / 1.0e9, Vs / 1.0e3, linestyle=":", linewidth=3.0)

    stv.set_state(1.0e10, 2000.0)
    print(f"Stv Vs at 10 GPa, 2000 K: {stv.shear_wave_velocity/1.e3:.3f} km/s")
    plt.show()


if __name__ == "__main__":
    check_iron_bearing_mineral_property_values()
    check_bcc_iron_consistency()
    check_fper_entropy()
    check_fig_1_fper_relaxed()
    check_fig_3_fcc_ferric_fper()
    check_fig_6a_iron_Cp_V()
    check_fig_6c_iron_phase_diagram()
    check_fig_7_fO2()
    check_fig_a1_fe2o3()
    check_fig_a4_stishovite()
