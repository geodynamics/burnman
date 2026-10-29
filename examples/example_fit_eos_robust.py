# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_fit_eos_robust
----------------------

This example demonstrates some of BurnMan's more advanced
least squares fitting capabilities for an equation of state (EoS).

A much more detailed example is provided in example_fit_eos.py.
This script focuses on what to do when an initial attempt at fitting
fails, when parameters are poorly known a-priori, or when there is a
risk of converging to a local minimum instead of the global minimum.

teaches:
- least squares fitting

"""

import numpy as np

import burnman

from burnman.utils.unitcell import molar_volume_from_unit_cell_volume
from burnman.utils.misc import pretty_print_values

if __name__ == "__main__":

    def make_stv():
        return burnman.Mineral(
            {
                "name": "stishovite",
                "formula": {"Si": 1, "O": 2},
                "molar_mass": 0.0600843,
                "V_0": 46.5e-6,
                "K_0": 300.0e9,
                "Kprime_0": 4.0,
                "equation_of_state": "bm3",
            }
        )

    print("Least squares equation of state fitting to unknown EoS\n")

    # Now let's fit some EoS data
    # Let's just take a bit of data from Andrault et al. (2003) on stishovite
    PV = np.array(
        [
            [0.0001, 46.5126, 0.0061],
            [5.037, 45.7902, 0.0053],
            [11.69, 44.947, 0.002],
            [17.67, 44.264, 0.002],
            [22.38, 43.776, 0.003],
            [29.38, 43.073, 0.009],
            [37.71, 42.278, 0.008],
            [46.03, 41.544, 0.017],
            [52.73, 40.999, 0.009],
        ]
    )

    Z = 2.0  # number of formula units per unit cell in stishovite
    PTV = np.array(
        [
            PV[:, 0] * 1.0e9,
            298.15 * np.ones_like(PV[:, 0]),
            molar_volume_from_unit_cell_volume(PV[:, 1], Z),
        ]
    ).T

    # Here, we assume that the pressure uncertainties are equal to 3% of the total pressure,
    # that the temperature uncertainties are negligible, and take the unit cell volume
    # uncertainties from the paper.
    # We also assume that the uncertainties in pressure and volume are uncorrelated.
    nul = np.zeros_like(PTV[:, 0])
    PTV_covariances = np.array(
        [
            [0.03 * PTV[:, 0], nul, nul],
            [nul, nul, nul],
            [nul, nul, molar_volume_from_unit_cell_volume(PV[:, 2], Z)],
        ]
    ).T
    PTV_covariances = np.power(PTV_covariances, 2.0)

    # Here's where we fit the data
    # The mineral parameters are automatically updated during fitting
    stv = make_stv()
    params = ["V_0", "K_0", "Kprime_0"]

    # Fit #1
    # No bounds are set for this initial fit.
    # In this case, the inversion will fail, because the
    # starting parameters are too far from the true values.

    # We here catch the exception.
    print("Fit #1 (no bounds)")
    try:
        fitted_eos = burnman.eos_fitting.fit_PTV_data(
            stv, params, PTV, PTV_covariances, verbose=False
        )
        raise ValueError("Fitting succeeded unexpectedly.")
    except Exception:
        print(
            "Fitting failed, because the starting "
            "parameters were too far from the true values. "
            "This led to one or more of the EoS parameters "
            "leaving their valid range.\n"
        )
        pass

    # Fit #2 (with priors)
    # Sometimes, we might have some prior knowledge about the parameters we are trying to fit.
    # We can use that prior knowledge to help the fitting process, but
    # we need to be aware that the priors will bias the fit towards the prior values.
    # In this case, we apply extremely weak priors on V_0,
    # and stronger priors on K_0 and Kprime_0.
    print("Fit #2 (with priors)")
    stv = make_stv()
    param_priors = np.array([46.5e-6, 300.0e9, 4.0])
    param_prior_inv_cov_matrix = np.diag(
        np.square(np.array([1.0 / 1.0e-3, 1.0 / 1.0e10, 1.0]))
    )
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        stv,
        params,
        PTV,
        PTV_covariances,
        param_priors=param_priors,
        param_prior_inv_cov_matrix=param_prior_inv_cov_matrix,
        verbose=False,
    )

    print("Optimized equation of state for stishovite:")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("")

    # Fit #3
    # Now we set bounds for the parameters, rather than using priors.
    print("Fit #3 (with bounds)")
    stv = make_stv()
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        stv,
        params,
        PTV,
        PTV_covariances,
        bounds=[(1.0e-10, np.inf), (1.0e9, np.inf), (2.0, 12.0)],
        verbose=False,
    )

    # Print the optimized parameters
    print("Optimized equation of state for stishovite:")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("")

    # Fit #4
    # You may also encounter situations where you don't know whether you have found the global minimum of the misfit function.
    # In this case, one effective strategy is to use differential evolution.

    # make the bounds 0.1 and 10. times the optimized parameters from the previous fit
    print("Fit #4 (with differential evolution)")
    bounds = [
        (0.1 * fitted_eos.popt[i], 10.0 * fitted_eos.popt[i])
        for i in range(len(fitted_eos.popt))
    ]
    stv = make_stv()
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        stv,
        params,
        PTV,
        PTV_covariances,
        bounds=bounds,
        first_run_differential_evolution=True,
        differential_evolution_seed=42,
        verbose=False,
    )

    # Print the optimized parameters
    print("Optimized equation of state for stishovite:")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("")
