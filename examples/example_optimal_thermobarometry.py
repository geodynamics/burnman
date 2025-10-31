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
from burnman import Composite
from burnman.minerals import mp50MnNCKFMASHTO, HP_2011_ds62
from burnman.optimize.composition_fitting import fit_composition_to_solution
from burnman.tools.thermobarometry import estimate_conditions

if __name__ == "__main__":

    # 1) Define observed mineral assemblage
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
    # Mg and Fe on the first two sites)
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
    fitted_species = ["Mn", "Fe", "Ti", "Mg", "Fe2+_A"]
    species_amounts = np.array([0.05, 1.0, 0.90, 0.05, 0.4])
    species_covariances = np.diag(np.array([0.01, 0.01, 0.01, 0.01, 0.2]) ** 2)
    species_conversions = {"Fe2+_A": {"Fea_A": 1.0}}
    popt, pcov, res = fit_composition_to_solution(
        ilmm, fitted_species, species_amounts, species_covariances, species_conversions
    )
    ilmm.set_composition(popt)
    ilmm.compositional_covariances = pcov

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
    print(
        f"Estimated Pressure: {assemblage.pressure/1.e9:.2f} "
        f"+/- {np.sqrt(res.xcov[0, 0])/1.e9:.2f} GPa"
    )
    print(
        f"Estimated Temperature: {assemblage.temperature:.2f} "
        f"+/- {np.sqrt(res.xcov[1, 1]):.2f} K"
    )
    print(f"Correlation between P and T: {res.xcorr:.4f}")
    print(f"Number of Reactions: {res.n_reactions}")
    print(f"Number of Parameters: {res.n_params}")
    print(f"Degrees of Freedom: {res.degrees_of_freedom}")
    print(f"Reduced Chi-squared: {res.reduced_chisqr:.4f}")
    print(f"Fit (sqrt reduced chi-squared): {res.fit:.4f}")
    print("Weighted reaction affinities:")
    np.set_printoptions(precision=2)
    print(res.weighted_affinities)
