# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_fit_eos_with_independent_volumes
----------------------------------------

This example demonstrates BurnMan's functionality to fit data to
an EoS of the user's choice using volumes as independent variables.

The first example deals with simple P(V, T) fitting.
The first example deals with E(V, T) fitting.

teaches:
- least squares fitting

"""
import numpy as np
import matplotlib.pyplot as plt

import burnman

from burnman.utils.unitcell import molar_volume_from_unit_cell_volume
from burnman.utils.misc import attribute_function
from burnman.utils.misc import pretty_print_values

if __name__ == "__main__":
    print(
        "Least squares equation of state fitting with volume as "
        "an independent variable.\n"
    )

    print("1) Fitting to stishovite room temperature P(V) data\n")
    # Let's fit an EoS to stishovite data from Andrault et al. (2003)
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

    Z = 2.0  # number of formula units per unit cell in stishovite
    VTP = np.array(
        [
            molar_volume_from_unit_cell_volume(PV[:, 1], Z),
            298.15 * np.ones_like(PV[:, 0]),
            PV[:, 0] * 1.0e9,
        ]
    ).T

    # Here, we assume that the pressure uncertainties are equal to 3% of the total pressure,
    # that the temperature uncertainties are negligible, and take the unit cell volume
    # uncertainties from the paper.
    # We also assume that the uncertainties in pressure and volume are uncorrelated.
    nul = np.zeros_like(VTP[:, 0])
    VTP_covariances = np.array(
        [
            [molar_volume_from_unit_cell_volume(PV[:, 2], Z), nul, nul],
            [nul, nul, nul],
            [nul, nul, 0.03 * VTP[:, 2]],
        ]
    ).T
    VTP_covariances = np.power(VTP_covariances, 2.0)

    # Here's where we fit the data
    # The mineral parameters are automatically updated during fitting
    stv = burnman.minerals.HP_2011_ds62.stv()
    params = ["V_0", "K_0", "Kprime_0"]
    fitted_eos = burnman.eos_fitting.fit_VTP_data(
        stv, params, VTP, VTP_covariances, verbose=False
    )

    # Print the optimized parameters
    print("Optimized equation of state for stishovite (HP):")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("")

    # Create a corner plot of the covariances
    fig = burnman.nonlinear_fitting.corner_plot(
        fitted_eos.popt, fitted_eos.pcov, params
    )
    plt.show()

    # Finally, let's plot our equation of state
    T = 298.15
    pressures = np.linspace(1.0e5, 60.0e9, 101)
    volumes = np.empty_like(pressures)

    VTPs = np.empty((len(pressures), 3))
    for i, P in enumerate(pressures):
        stv.set_state(P, T)
        VTPs[i] = [stv.V, stv.temperature, stv.pressure]

    # Plot the 95% confidence and prediction bands
    cp_bands = burnman.nonlinear_fitting.confidence_prediction_bands(
        model=fitted_eos,
        x_array=VTPs,
        confidence_interval=0.95,
        f=attribute_function(stv, "P", using_volume=True),
        flag="P",
    )

    plt.fill_between(
        VTPs[:, 0] * 1.0e6,
        cp_bands[0] / 1.0e9,
        cp_bands[1] / 1.0e9,
        alpha=0.3,
        linewidth=0.0,
        color="r",
        label="95% confidence bands",
    )
    plt.fill_between(
        VTPs[:, 0] * 1.0e6,
        cp_bands[2] / 1.0e9,
        cp_bands[3] / 1.0e9,
        alpha=0.3,
        linewidth=0.0,
        color="b",
        label="95% prediction bands",
    )

    plt.plot(
        VTPs[:, 0] * 1.0e6, VTPs[:, 2] / 1.0e9, label="Optimized fit for stishovite"
    )
    plt.errorbar(
        VTP[:, 0] * 1.0e6,
        VTP[:, 2] / 1.0e9,
        xerr=np.sqrt(VTP_covariances.T[0][0]) * 1.0e6,
        yerr=np.sqrt(VTP_covariances.T[2][2]) / 1.0e9,
        linestyle="None",
        marker="o",
        label="Andrault et al. (2003)",
    )

    plt.xlabel("Volume (cm^3/mol)")
    plt.ylabel("Pressure (GPa)")
    plt.legend(loc="upper right")
    plt.title("Stishovite EoS (room temperature)")
    plt.show()

    # Now let's fit the same parameters for the SLB EoS.
    stv_SLB = burnman.minerals.SLB_2011.stishovite()

    # Removing the property modifiers makes fitting the SLB EoS
    # a lot faster as we avoid an iterative solve
    stv_SLB.property_modifiers = []

    fitted_eos_SLB = burnman.eos_fitting.fit_VTP_data(
        stv_SLB, params, VTP, VTP_covariances, verbose=False
    )
    # Print the optimized parameters
    print("Optimized equation of state for stishovite (SLB):")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("")

    # 2) Now let's fit some F(V, T) data
    print("2) Fitting to synthetic F(V, T) data\n")

    # Create the mineral instance that we'll optimise
    per_SLB = burnman.minerals.SLB_2011.periclase()

    V_0 = per_SLB.params["V_0"]
    volumes = np.linspace(0.7 * V_0, V_0, 4)
    temperatures = np.linspace(500.0, 2000.0, 4)
    VV, TT = np.meshgrid(volumes, temperatures)
    V_data = VV.flatten()
    T_data = TT.flatten()
    F_data = per_SLB.evaluate_with_volumes(["helmholtz"], V_data, T_data)[0]

    # Add some random noise
    np.random.seed(100)
    F_stdev = 10.0
    F_data = F_data + np.random.normal(scale=F_stdev, size=F_data.shape)

    # Collect all the objects we need to supply as arguments to the
    # fitting function
    VTF_data = np.array([V_data, T_data, F_data]).T
    nul = 0.0 * V_data
    Fcov = F_stdev * F_stdev * np.ones_like(nul)
    VTF_covariances = np.array([[nul, nul, nul], [nul, nul, nul], [nul, nul, Fcov]]).T
    flags = ["helmholtz"] * len(F_data)
    fit_params = ["V_0", "K_0", "Kprime_0", "grueneisen_0", "q_0", "Debye_0", "F_0"]

    print("Actual equation of state:")
    for p in fit_params:
        print(f"{p}: {per_SLB.params[p]:.4e}")
    print("")

    # Let's change some parameters in our equation of state
    # so that fitting is not trivial
    per_SLB.params["F_0"] = per_SLB.params["F_0"] + 2300.0
    per_SLB.params["V_0"] = per_SLB.params["V_0"] * 0.9
    per_SLB.params["K_0"] = per_SLB.params["K_0"] * 0.9

    # And now we can fit our data!
    fitted_eos = burnman.eos_fitting.fit_VTp_data(
        mineral=per_SLB,
        flags=flags,
        fit_params=fit_params,
        data=VTF_data,
        data_covariances=VTF_covariances,
        verbose=False,
    )

    # Print the optimized parameters
    print("Optimized equation of state:")
    pretty_print_values(fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params)
    print("\nGoodness of fit:")
    print(fitted_eos.goodness_of_fit)
    print("")
