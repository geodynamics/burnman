# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_optimal_thermobarometry
-------------------------------

This example script is intended to demonstrate optimal thermobarometry
using BurnMan. The technique is based on the work of Powell and Holland
(1994; American Mineralogist 79 (1-2): 120-133).
We cover importing BurnMan modules, creating a composite
material representing a mineral assemblage, fitting the compositions of
the constituent minerals to their solution models, and estimating the
pressure and temperature conditions of equilibration based on the
mineral compositions and their uncertainties. Finally, we print
the estimated conditions along with their uncertainties and correlations.

*Uses:*

* :doc:`mineral_database`
* :func:`burnman.optimize.composition_fitting.fit_composition_to_solution`
* :func:`burnman.optimize.thermobarometry.estimate_conditions`


*Demonstrates:*

* creating mineral assemblages
* fitting mineral compositions to solution models
* estimating pressure and temperature conditions of equilibration
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings

from burnman import Composite
from burnman.minerals import mp50MnNCKFMASHTO, HP_2011_ds62, SLB_2024
from burnman.optimize.composition_fitting import fit_composition_to_solution
from burnman.tools.thermobarometry import estimate_conditions
from burnman.tools.thermobarometry import (
    get_reaction_matrix,
    assemblage_set_state_from_params,
    assemblage_affinity_misfit,
)
from burnman.utils.chemistry import formula_to_string
from burnman import equilibrate
from burnman.optimize.nonlinear_fitting import plot_cov_ellipse

if __name__ == "__main__":

    # Example 1: Sillimanite-Kyanite equilibrium
    # In this example, we will estimate the temperature of
    # the sillimanite-kyanite equilibrium at known pressure,
    # taking into account uncertainties in the
    # thermodynamic dataset.
    sill = HP_2011_ds62.sill()
    ky = HP_2011_ds62.ky()
    assemblage = Composite([sill, ky])

    temperatures = np.linspace(600, 800, 5)
    equality_constraints = [["T", temperatures], ["phase_fraction", (sill, 0.0)]]
    sols, prm = equilibrate(
        ky.formula, assemblage, equality_constraints=equality_constraints
    )
    Ps = np.array([sol.assemblage.pressure for sol in sols])
    Ts = np.array([sol.assemblage.temperature for sol in sols])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = estimate_conditions(
            assemblage,
            dataset_covariances=HP_2011_ds62.cov(),
            guessed_conditions=np.array([0.5e9, 500.0]),
            pressure_bounds=[0.1e9, 0.1e9],
        )

    print("Estimated conditions for sillimanite-kyanite equilibrium:")
    print(f"    Temperature: {res.x[1]:.0f} +/- {np.sqrt(res.xcov[1, 1]):.0f} K")
    print(f"    (fixed pressure at {res.x[0]/1.e9} GPa)")

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(Ps / 1.0e9, Ts, label="Equilibrium sill+ky")
    ax.scatter(res.x[0] / 1.0e9, res.x[1], color="red", label="Estimated conditions")
    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    plt.show()

    # Example 2: Olivine-Wadsleyite equilibrium with compositional fitting
    # In this example, we will estimate the pressure and temperature of
    # the olivine-wadsleyite equilibrium based on fitted compositions
    # of the two minerals, taking into account uncertainties in the
    # thermodynamic dataset as well as uncertainties in the fitted
    # compositions, and also a prior on the conditions themselves.
    ol = SLB_2024.olivine()
    wad = SLB_2024.wadsleyite()
    assemblage = Composite([ol, wad])

    f_fwd = 0.18
    ol.set_composition([0.9, 0.1])
    wad.set_composition([1.0 - f_fwd, f_fwd])

    ol.compositional_covariances = np.diag(np.array([0.00001, 0.00001]) ** 2)
    wad.compositional_covariances = np.diag(np.array([0.01, 0.01]) ** 2)

    # The covariance matrix here is not very realistic - the correlations
    # between the endmembers of olivine and wadsleyite are not negligible.
    # Here, we apply 70 J/mol uncorrelated uncertainty to each endmember
    # to illustrate the method.
    assemblage.state_priors = np.array([13.0e9, 1400.0])
    assemblage.state_inverse_covariances = np.array(
        [[1.0e-9**2, 0.0], [0.0, 1.0e-2**2]]
    )
    dataset_covariances = {
        "endmember_names": [*ol.endmember_names, *wad.endmember_names],
        "covariance_matrix": np.eye(4) * 70.0 * 70.0,
    }
    res = estimate_conditions(
        assemblage,
        dataset_covariances=dataset_covariances,
        include_state_misfit=True,
        guessed_conditions=np.array([14.0e9, 1400.0]),
    )
    print("\nEstimated conditions for olivine-wadsleyite equilibrium:")
    print(
        f"    Pressure: {res.x[0]/1.e9:.1f} +/- {np.sqrt(res.xcov[0, 0])/1.e9:.1f} GPa"
    )
    print(
        f"    Temperature: {int(round(res.x[1], -1))} +/- {np.sqrt(res.xcov[1, 1]):.0f} K"
    )

    f_fwds = np.linspace(f_fwd - 0.02, f_fwd + 0.02, 5)
    assemblage.set_state(10.0e9, 1400.0)
    equality_constraints = [
        ["phase_fraction", (ol, 1.0)],
        [
            "phase_composition",
            (wad, (["Mg_A", "Fe_A"], [0.0, 1.0], [1.0, 1.0], f_fwds)),
        ],
    ]
    sols, prm = equilibrate(ol.formula, assemblage, equality_constraints)

    Ps = np.array([sol.assemblage.pressure for sol in sols])
    Ts = np.array([sol.assemblage.temperature for sol in sols])

    # Plot the uncertainty ellipse
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Plot the prior uncertainty
    prior_cov = np.linalg.inv(assemblage.state_inverse_covariances)
    prior_cov[0, :] /= 1.0e9
    prior_cov[:, 0] /= 1.0e9
    plot_cov_ellipse(
        prior_cov,
        assemblage.state_priors / np.array([1.0e9, 1.0]),
        nstd=2,
        ax=ax,
        facecolor="none",
        edgecolor="black",
        linestyle="--",
        label="2-sigma prior uncertainty",
    )

    # Plot the equilibrium curve and estimated conditions
    ax.scatter(Ps[2] / 1.0e9, Ts[2], s=40, color="purple")
    ax.plot(Ps / 1.0e9, Ts, color="purple")
    ax.scatter(
        Ps / 1.0e9,
        Ts,
        s=10,
        color="purple",
        label="equilibrium (wd $x_{{{Fe}}}$=0.16-0.20)",
    )

    # Plot the estimated conditions and uncertainty
    ax.scatter(res.x[0] / 1.0e9, res.x[1], color="red", label="estimated conditions")
    cov = res.xcov.copy()
    cov[0, :] /= 1.0e9
    cov[:, 0] /= 1.0e9
    plot_cov_ellipse(
        cov,
        res.x / np.array([1.0e9, 1.0]),
        nstd=2,
        ax=ax,
        facecolor="none",
        edgecolor="red",
        label="2-sigma uncertainty",
    )

    ax.set_xlim(10, 16)
    ax.set_ylim(800, 2000)
    ax.set_xlabel("Pressure (GPa)")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    plt.show()

    # Example 3: Complex assemblage with multiple solid solutions
    # In this example, we will estimate the pressure and temperature of
    # equilibration of a complex mineral assemblage based on fitted
    # compositions of the constituent minerals, taking into account
    # uncertainties in the thermodynamic dataset as well as uncertainties
    # in the fitted compositions.

    # 1) Define the observed mineral assemblage
    mu = mp50MnNCKFMASHTO.mu()
    bi = mp50MnNCKFMASHTO.bi()
    g = mp50MnNCKFMASHTO.g()
    ilmm = mp50MnNCKFMASHTO.ilmm()
    st = mp50MnNCKFMASHTO.st()
    q = HP_2011_ds62.q()

    assemblage = Composite([mu, bi, g, ilmm, st, q])

    # 2) Fit the measured compositions (in molar amounts) to the solution
    # models for each mineral in the assemblage.
    # The amounts do not need to sum to any particular number,
    # as the fitting function will normalize them appropriately.
    # The compositions correspond approximately to the expected
    # compositions at about 4 kbar and 873 K.

    # a) Muscovite
    fitted_species = ["Na", "Ca", "K", "Fe", "Mg", "Al", "Si"]
    species_amounts = np.array([0.40, 0.01, 0.55, 0.01, 0.01, 3.00, 3.00])
    species_covariances = np.diag(
        np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) ** 2
    )
    popt, pcov, res = fit_composition_to_solution(
        mu, fitted_species, species_amounts, species_covariances
    )
    mu.set_composition(popt)
    mu.compositional_covariances = pcov

    # b) Biotite (requires parameter constraining order state of
    # Mg and Fe on the first two sites).
    # Here, we will allow the Mg-Fe ordering to vary to improve the fit
    # by adding a free compositional vector.
    fitted_species = ["Mn", "Fe", "Mg", "Al", "Si", "Ti", "Mg_M3"]
    species_amounts = np.array([0.01, 1.50, 1.00, 1.65, 2.65, 0.20, 0.10])
    species_covariances = np.diag(
        np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1]) ** 2
    )
    variable_conversions = {"Mg_M3": {"Mgmthree_A": 1.0}}
    popt, pcov, res = fit_composition_to_solution(
        bi, fitted_species, species_amounts, species_covariances, variable_conversions
    )
    bi.set_composition(popt)
    bi.compositional_covariances = pcov
    bi.free_compositional_vectors = Composite([bi]).reaction_basis

    # c) Garnet
    fitted_species = ["Mn", "Fe", "Mg", "Ca", "Al", "Si"]
    species_amounts = np.array([0.25, 2.30, 0.40, 0.05, 2.00, 3.00])
    species_covariances = np.diag(np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) ** 2)
    popt, pcov, res = fit_composition_to_solution(
        g, fitted_species, species_amounts, species_covariances
    )
    g.set_composition(popt)
    g.compositional_covariances = pcov

    # d) Ilmenite (requires parameter constraining order state of Fe and Ti)
    # Here, we will allow the Fe-Ti ordering to vary to improve the fit
    # by adding a free compositional vector.
    fitted_species = ["Mn", "Fe", "Ti", "Mg", "Fe2+_A"]
    species_amounts = np.array([0.05, 1.0, 0.90, 0.05, 0.4])
    species_covariances = np.diag(np.array([0.01, 0.01, 0.01, 0.01, 0.2]) ** 2)
    species_conversions = {"Fe2+_A": {"Fea_A": 1.0}}
    popt, pcov, res = fit_composition_to_solution(
        ilmm, fitted_species, species_amounts, species_covariances, species_conversions
    )
    ilmm.set_composition(popt)
    ilmm.compositional_covariances = pcov
    ilmm.free_compositional_vectors = Composite([ilmm]).reaction_basis

    # e) Staurolite
    fitted_species = ["Mn", "Fe", "Mg", "Al", "Si", "Ti"]
    species_amounts = np.array([0.05, 3.30, 0.72, 17.78, 7.50, 0.12])
    species_covariances = np.diag(np.array([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]) ** 2)
    popt, pcov, res = fit_composition_to_solution(
        st, fitted_species, species_amounts, species_covariances
    )
    st.set_composition(popt)
    st.compositional_covariances = pcov

    # f) Quartz (no fitting needed)

    # 3) Estimate the pressure and temperature conditions of equilibration
    #    based on the fitted compositions and their uncertainties, and
    #    also the endmember covariances provided by the underlying dataset.
    res = estimate_conditions(
        assemblage,
        dataset_covariances=HP_2011_ds62.cov(),
        guessed_conditions=np.array([0.5e9, 500.0]),
        small_fraction_tol=0.01,
    )

    # 4) Print the estimated conditions along with their uncertainties
    #    and correlations.
    assemblage_name = "-".join([phase.name for phase in assemblage.phases])
    print(f"\nEstimated conditions for {assemblage_name} assemblage:")
    print(
        f"    Pressure: {assemblage.pressure/1.e9:.2f} "
        f"+/- {np.sqrt(res.xcov[0, 0])/1.e9:.2f} GPa"
    )
    print(
        f"    Temperature: {assemblage.temperature:.0f} "
        f"+/- {np.sqrt(res.xcov[1, 1]):.0f} K"
    )
    print(f"    Correlation between P and T: {res.xcorr[0][1]:.2f}")
    print(f"    Number of Reactions: {res.n_reactions}")
    print(f"    Number of Parameters: {res.n_params}")
    print(f"    Degrees of Freedom: {res.degrees_of_freedom}")
    print(f"    Reduced Chi-squared: {res.reduced_chisqr:.2f}")
    print(f"    Fit (sqrt reduced chi-squared): {res.fit:.2f}")

    # The names of the sites in these models are a bit cumbersome,
    # so we will do some string replacements to make them more
    # readable in the output.
    replace_strings = [
        ["Fethree", "Fe3+"],
        ["monetwo", ""],
        ["mtwoa", ""],
        ["mtwob", ""],
        ["mthree", ""],
        ["x", ""],
        ["y", ""],
        ["tone", ""],
        ["t", ""],
        ["Ka", "K"],
        ["Naa", "Na"],
        ["Caa", "Ca"],
        ["Fea", "Fe"],
        ["Mga", "Mg"],
        ["Mna", "Mn"],
        ["Tia", "Ti"],
        ["Fe3+a", "Fe3+"],
        ["Oh", "(OH)"],
        ["b", ""],
        ["v", ""],
    ]

    print("\n    Mineral compositions and site occupancies:")
    for phase in assemblage.phases:
        np.set_printoptions(precision=2, formatter={"all": lambda x: f"{x:0.2f}"})
        print(f"    {phase.name}:")

        if hasattr(phase, "free_compositional_vectors"):
            free_site_occupancy_vectors = phase.free_compositional_vectors.dot(
                phase.solution_model.endmember_occupancies
            )
            for v in free_site_occupancy_vectors:
                site_formula = phase.site_formula(
                    precision=2, site_occupancies=v, print_absent_species=False
                )
                for old, new in replace_strings:
                    site_formula = site_formula.replace(old, new)
                print(f"        free site occupancy vector: {site_formula}")

        print(
            f"        composition: {formula_to_string(phase.formula, use_fractions=False)}"
        )
        if hasattr(phase, "molar_fractions"):
            site_formula = phase.site_formula(2)
            for old, new in replace_strings:
                site_formula = site_formula.replace(old, new)
            print(f"        site occupancies: {site_formula}")
    print("")

    # Let's also have a look at how the misfit changes when varying the
    # compositional degrees of freedom.
    R = get_reaction_matrix(assemblage, small_fraction_tol=0.01)

    x = np.linspace(-0.1, 0.1, 51)
    fig, ax = plt.subplots()
    labels = ["bi Fe-Mg order", "ilmm Fe-Ti order"]
    print(
        "    Reduced chi-squared misfit variation with compositional degrees of freedom:"
    )
    for i in range(res.n_params - 2):
        chisqr_vals = []
        for xi in x:
            params = res.x.copy()
            params[i + 2] += xi
            assemblage_set_state_from_params(assemblage, params)
            chisqr = assemblage_affinity_misfit(
                assemblage, R, dataset_covariances=HP_2011_ds62.cov()
            )
            chisqr_vals.append(chisqr)
        red_chisqr = np.array(chisqr_vals) / (res.n_reactions - res.n_params)
        ax.plot(x, red_chisqr, label=f"{labels[i]}")
        print(
            f"    {labels[i]}:\n"
            f"        chisq(x={x[0]})={red_chisqr[0]:.2f},\n"
            f"        chisq(x={x[25]})={red_chisqr[25]:.2f},\n"
            f"        chisq(x={x[-1]})={red_chisqr[-1]:.2f}\n"
        )

    ax.set_xlabel("Change in order relative to optimal value")
    ax.set_ylabel("Reduced chi-squared misfit")
    ax.legend()
    plt.show()

    # Uncomment to see the weighted reaction affinities
    # print("Weighted reaction affinities:")
    # np.set_printoptions(precision=2)
    # print(res.weighted_affinities)
