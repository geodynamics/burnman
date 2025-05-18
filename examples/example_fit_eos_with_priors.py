# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_fit_eos_with_priors
---------------------------

This example demonstrates BurnMan's functionality to fit data to
an EoS of the user's choice with (potentially correlated)
Gaussian priors on the parameters.

This provides an alternative to artificially fixing the parameters.

teaches:
- least squares fitting with Gaussian priors on the parameter values.

"""
import numpy as np
import matplotlib.pyplot as plt

import burnman

from burnman.utils.unitcell import molar_volume_from_unit_cell_volume
from burnman.utils.misc import attribute_function
from burnman.utils.misc import pretty_print_values

if __name__ == "__main__":
    print("Least squares equation of state fitting with Gaussian priors\n")

    # Let's fit some of the high pressure stishovite data from Andrault et al. (2003)
    PV = np.array(
        [
            [0.0001, 46.5126, 0.0061],
            [1.168, 46.3429, 0.0053],
            [2.299, 46.1756, 0.0043],
            [3.137, 46.0550, 0.0051],
            [4.252, 45.8969, 0.0045],
            [5.037, 45.7902, 0.0053],
            [5.851, 45.6721, 0.0038],
            [6.613, 45.5715, 0.0050],
            [7.504, 45.4536, 0.0041],
            [8.264, 45.3609, 0.0056],
            [9.635, 45.1885, 0.0042],
            [11.69, 44.947, 0.002],
            [17.67, 44.264, 0.002],
            [22.38, 43.776, 0.003],
            [29.38, 43.073, 0.009],
            [37.71, 42.278, 0.008],
            [46.03, 41.544, 0.017],
            [52.73, 40.999, 0.009],
            [26.32, 43.164, 0.006],
            [30.98, 42.772, 0.005],
            [34.21, 42.407, 0.003],
            [38.45, 42.093, 0.004],
            [43.37, 41.610, 0.004],
            [47.49, 41.280, 0.007],
        ]
    )

    PV[:, 0] = PV[:, 0] * 1.0e9
    PV[:, 1] = molar_volume_from_unit_cell_volume(PV[:, 1], 2.0)
    PV[:, 2] = molar_volume_from_unit_cell_volume(PV[:, 2], 2.0)

    # Use only data points where P > 15 GPa in the fitting
    mask = PV[:, 0] > 15.0e9
    PTV = np.array(
        [
            PV[mask, 0],
            PV[mask, 0] * 0.0 + 298.15,
            PV[mask, 1],
        ]
    ).T

    nul = 0.0 * PTV.T[0]
    PTV_covariances = np.array(
        [
            [0.03 * PTV.T[0], nul, nul],
            [nul, nul, nul],
            [nul, nul, PV[mask, 2]],
        ]
    ).T
    PTV_covariances = np.power(PTV_covariances, 2.0)

    # The mineral parameters are automatically updated during fitting
    stv = burnman.minerals.HP_2011_ds62.stv()
    params = ["V_0", "K_0", "Kprime_0"]
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        stv, params, PTV, PTV_covariances, verbose=False
    )

    # Evaluate some volumes for plotting later
    pressures = np.linspace(1.0e5, 100.0e9, 101)
    temperatures = 300.0 + 0.0 * pressures
    V_fit = stv.evaluate(["V"], pressures, temperatures)[0]
    confidence_interval = 0.95
    cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(
        model=fitted_eos,
        x_array=np.array([pressures, temperatures, V_fit]).T,
        confidence_interval=0.95,
        f=attribute_function(stv, "V"),
        flag="V",
    )

    # Print the optimized parameters
    print("Optimized equation of state for stishovite with no priors on parameters:")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print(f"Weighted sum of squares: {fitted_eos.WSS:.2f}")
    print(f"Goodness of fit: {fitted_eos.goodness_of_fit:.2f}")
    print("")

    # Now define a weak prior on K_0' while leaving V_0 and K_0 free
    prior_values = np.array([0.0, 0.0, 4.0])
    prior_cov = np.power(np.diag([np.inf, np.inf, 1.0]), 2.0)
    prior_inv_cov = np.linalg.inv(prior_cov)

    # Refit
    fitted_eos = burnman.eos_fitting.fit_PTV_data(
        stv,
        params,
        PTV,
        PTV_covariances,
        param_priors=prior_values,
        param_prior_inv_cov_matrix=prior_inv_cov,
        verbose=False,
    )
    V_fit_with_priors = stv.evaluate(["V"], pressures, temperatures)[0]

    print("Optimized equation of state for stishovite with weak prior; K' = 4(1):")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print(f"Weighted sum of squares: {fitted_eos.WSS:.2f}")
    print(f"Goodness of fit: {fitted_eos.goodness_of_fit:.2f}")
    print("")

    cp_bands_with_priors = burnman.nonlinear_fitting.confidence_prediction_bands(
        model=fitted_eos,
        x_array=np.array([pressures, temperatures, V_fit_with_priors]).T,
        confidence_interval=0.95,
        f=attribute_function(stv, "V"),
        flag="V",
    )

    fig = plt.figure()
    ax = [fig.add_subplot(1, 1, i) for i in range(1, 2)]

    ax[0].scatter(
        PV[~mask, 0] / 1.0e9,
        PV[~mask, 1] * 1.0e6,
        label="Andrault et al., 2003 (not fitted)",
    )
    ax[0].scatter(
        PV[mask, 0] / 1.0e9, PV[mask, 1] * 1.0e6, label="Andrault et al., 2003 (fitted)"
    )
    ax[0].errorbar(
        PV[:, 0] / 1.0e9,
        PV[:, 1] * 1.0e6,
        xerr=0.03 * PV[:, 0] / 1.0e9,
        yerr=PV[:, 2] * 1.0e6,
        linestyle="",
    )

    ax[0].plot(
        pressures / 1.0e9, V_fit * 1.0e6, label="Fit without priors", color="orange"
    )
    ax[0].fill_between(
        pressures / 1.0e9,
        cp_bands[0] * 1.0e6,
        cp_bands[1] * 1.0e6,
        alpha=0.15,
        label="95% confidence interval",
        color="orange",
    )

    ax[0].plot(
        pressures / 1.0e9,
        V_fit_with_priors * 1.0e6,
        label="Fit with weak prior (K'=4$\\pm$1)",
        color="blue",
    )
    ax[0].fill_between(
        pressures / 1.0e9,
        cp_bands_with_priors[0] * 1.0e6,
        cp_bands_with_priors[1] * 1.0e6,
        alpha=0.15,
        label="95% confidence interval",
        color="blue",
    )

    ax[0].set_xlabel("Pressure (GPa)")
    ax[0].set_ylabel("Volume (cm$^3$/mol)")
    ax[0].legend()
    plt.show()
