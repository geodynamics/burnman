# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_fit_data
----------------

This example demonstrates BurnMan's functionality to fit various
mineral physics data to an EoS of the user's choice.

Please note also the separate file example_fit_eos.py, which can be viewed
as a more advanced example in the same general field.

teaches:
- least squares fitting

"""
import numpy as np
import matplotlib.pyplot as plt

import burnman
import warnings


if __name__ == "__main__":
    # 1) Fitting shear modulus and its derivative to shear wave velocity data
    print("1) Fitting shear modulus and its derivative to shear wave velocity data\n")

    # First, read in the data from file and convert to SI units.
    PTp_data = np.loadtxt("../burnman/data/input_minphys/Murakami_perovskite.txt")
    PTp_data[:, 0] = PTp_data[:, 0] * 1.0e9
    PTp_data[:, 2] = PTp_data[:, 2] * 1.0e3

    # Make the test mineral
    mg_perovskite_test = burnman.Mineral()
    mg_perovskite_test.params = {
        "V_0": 24.45e-6,
        "K_0": 281.0e9,
        "Kprime_0": 4.1,
        "molar_mass": 0.10,
        "G_0": 200.0e9,
        "Gprime_0": 2.0,
    }

    def best_fit():
        return burnman.eos_fitting.fit_PTp_data(
            mineral=mg_perovskite_test,
            flags="shear_wave_velocity",
            fit_params=["G_0", "Gprime_0"],
            data=PTp_data,
            verbose=False,
        )

    pressures = np.linspace(1.0e5, 150.0e9, 101)
    temperatures = pressures * 0.0 + 300.0

    # Fit to the second order Birch-Murnaghan EoS
    mg_perovskite_test.set_method("bm3shear2")
    fitted_eos = best_fit()
    print("2nd order fit:")
    burnman.utils.misc.pretty_print_values(
        fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params
    )
    model_vs_2nd_order_correct = mg_perovskite_test.evaluate(
        ["shear_wave_velocity"], pressures, temperatures
    )[0]

    with warnings.catch_warnings(record=True) as w:
        mg_perovskite_test.set_method("bm3")
        print(w[-1].message)
    model_vs_2nd_order_incorrect = mg_perovskite_test.evaluate(
        ["shear_wave_velocity"], pressures, temperatures
    )[0]
    print("")

    # Fit to the third order Birch-Murnaghan EoS
    mg_perovskite_test.set_method("bm3")
    fitted_eos = best_fit()
    print("3rd order fit:")
    burnman.utils.misc.pretty_print_values(
        fitted_eos.popt, fitted_eos.pcov, fitted_eos.fit_params
    )
    model_vs_3rd_order_correct = mg_perovskite_test.evaluate(
        ["shear_wave_velocity"], pressures, temperatures
    )[0]

    with warnings.catch_warnings(record=True) as w:
        mg_perovskite_test.set_method("bm3shear2")
        print(w[-1].message)
    model_vs_3rd_order_incorrect = mg_perovskite_test.evaluate(
        ["shear_wave_velocity"], pressures, temperatures
    )[0]
    print("")

    plt.plot(
        pressures / 1.0e9,
        model_vs_2nd_order_correct / 1000.0,
        color="r",
        linestyle="-",
        linewidth=2,
        label="Correct 2nd order fit",
    )
    plt.plot(
        pressures / 1.0e9,
        model_vs_2nd_order_incorrect / 1000.0,
        color="r",
        linestyle="-.",
        linewidth=2,
        label="Incorrect 2nd order fit",
    )
    plt.plot(
        pressures / 1.0e9,
        model_vs_3rd_order_correct / 1000.0,
        color="b",
        linestyle="-",
        linewidth=2,
        label="Correct 3rd order fit",
    )
    plt.plot(
        pressures / 1.0e9,
        model_vs_3rd_order_incorrect / 1000.0,
        color="b",
        linestyle="-.",
        linewidth=2,
        label="Incorrect 3rd order fit",
    )
    plt.scatter(PTp_data[:, 0] / 1.0e9, PTp_data[:, 2] / 1.0e3)
    plt.ylim([6.55, 8])
    plt.xlim([25.0, 135.0])
    plt.ylabel("Shear velocity (km/s)")
    plt.xlabel("Pressure (GPa)")
    plt.legend(loc="lower right", prop={"size": 12}, frameon=False)
    plt.savefig("output_figures/example_fit_data1.png")
    plt.show()

    # 2) Fitting standard enthalpy and heat capacity to enthalpy data
    print("2) Fitting standard enthalpy and heat capacity to enthalpy data\n")

    per_SLB = burnman.minerals.SLB_2011.periclase()
    per_HP = burnman.minerals.HP_2011_ds62.per()

    # this is the mineral we'll optimise
    per_opt = burnman.minerals.HP_2011_ds62.per()

    # Load some example enthalpy data
    TH_data = np.loadtxt(
        "../burnman/data/input_fitting/" "Victor_Douglas_1963_deltaH_MgO.dat"
    )
    per_HP.set_state(1.0e5, 298.15)
    PTH_data = np.array(
        [TH_data[:, 0] * 0.0 + 1.0e5, TH_data[:, 0], TH_data[:, 2] * 4.184 + per_HP.H]
    ).T
    nul = TH_data[:, 0] * 0.0
    PTH_covariances = np.array(
        [
            [nul, nul, nul],
            [nul, TH_data[:, 1], nul],
            [nul, nul, np.power(TH_data[:, 2] * 4.184 * 0.0004, 2.0)],
        ]
    ).T

    per_opt.params["S_0"] = 6.439 * 4.184
    model = burnman.eos_fitting.fit_PTp_data(
        mineral=per_opt,
        flags="H",
        fit_params=["H_0", "Cp"],
        data=PTH_data,
        data_covariances=PTH_covariances,
        max_lm_iterations=10,
        verbose=False,
    )

    print("Optimised values:")
    params = ["H_0", "Cp_a", "Cp_b", "Cp_c", "Cp_d"]
    burnman.utils.misc.pretty_print_values(model.popt, model.pcov, params)
    print("")

    # Corner plot
    fig = burnman.nonlinear_fitting.corner_plot(model.popt, model.pcov, params)
    plt.savefig("output_figures/example_fit_data2.png")
    plt.show()

    # Plot models
    temperatures = np.linspace(200.0, 2000.0, 101)
    pressures = np.array([298.15] * len(temperatures))
    plt.plot(
        temperatures,
        per_HP.evaluate(["molar_heat_capacity_p"], pressures, temperatures)[0],
        linestyle="--",
        label="HP",
    )
    plt.plot(
        temperatures,
        per_SLB.evaluate(["molar_heat_capacity_p"], pressures, temperatures)[0],
        linestyle="--",
        label="SLB",
    )
    plt.plot(
        temperatures,
        per_opt.evaluate(["molar_heat_capacity_p"], pressures, temperatures)[0],
        label="Optimised fit",
    )

    plt.legend(loc="lower right")
    plt.xlim(0.0, temperatures[-1])
    plt.xlabel("Temperature (K)")
    plt.ylabel("Heat capacity (J/K/mol)")
    plt.legend(loc="lower right", prop={"size": 12}, frameon=False)
    plt.savefig("output_figures/example_fit_data3.png")
    plt.show()
