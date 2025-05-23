# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_optimize_pv
-------------------

Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM. For more extensive comments on this setup, see tutorial/step_2.py

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Composite`
* :class:`burnman.seismic.PREM`
* :func:`burnman.geotherm.brown_shankland`
* :func:`burnman.Material.evaluate`
* :func:`burnman.utils.math.l2_norm_profiles`

*Demonstrates:*

* compare errors between models
* loops over models

"""
import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman import minerals


if __name__ == "__main__":
    # Define reference model and depth to evaluate
    seismic_model = burnman.seismic.PREM()
    number_of_points = 20
    depths = np.linspace(700e3, 2800e3, number_of_points)
    (
        seis_p,
        seis_rho,
        seis_vp,
        seis_vs,
        seis_vphi,
        seis_K,
        seis_G,
    ) = seismic_model.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi", "K", "G"], depths
    )

    # Define geotherm
    temperature = burnman.geotherm.brown_shankland(depths)

    # Define solid solutions
    perovskite = minerals.SLB_2011.mg_fe_perovskite()
    # Set molar_fraction of fe_perovskite and al_perovskite:
    perovskite.set_composition([0.94, 0.06, 0.0])
    ferropericlase = minerals.SLB_2011.ferropericlase()
    # Set molar_fraction of MgO and FeO:
    ferropericlase.set_composition([0.8, 0.2])

    def material_error(amount_perovskite):
        # Define composite using the values
        rock = burnman.Composite(
            [perovskite, ferropericlase], [amount_perovskite, 1.0 - amount_perovskite]
        )

        # Compute velocities
        mat_rho, mat_vp, mat_vs, mat_vphi, mat_K, mat_G = rock.evaluate(
            ["density", "v_p", "v_s", "v_phi", "K_S", "G"], seis_p, temperature
        )

        print("Calculations are done for:")
        rock.debug_print()
        # Calculate errors
        [vs_err, vphi_err, rho_err, K_err, G_err] = np.square(
            burnman.utils.math.l2_norm_profiles(
                depths,
                [mat_vs, mat_vphi, mat_rho, mat_K, mat_G],
                [seis_vs, seis_vphi, seis_rho, seis_K, seis_G],
            )
        )
        # Normalize errors
        vs_err = vs_err / np.mean(seis_vs) ** 2.0
        vphi_err = vphi_err / np.mean(seis_vphi) ** 2.0
        rho_err = rho_err / np.mean(seis_rho) ** 2.0
        K_err = K_err / np.mean(seis_K) ** 2.0
        G_err = G_err / np.mean(seis_G) ** 2.0
        return vs_err, vphi_err, rho_err, K_err, G_err

    # Run through fractions of perovskite
    xx = np.linspace(0.0, 1.0, 40)
    errs = np.array([material_error(x) for x in xx])

    # Plot results
    yy_vs = errs[:, 0]
    yy_vphi = errs[:, 1]
    yy_rho = errs[:, 2]
    yy_K = errs[:, 3]
    yy_G = errs[:, 4]
    plt.plot(xx, yy_vs, "r-x", label=("vs error"))
    plt.plot(xx, yy_vphi, "b-x", label=("vphi error"))
    plt.plot(xx, yy_rho, "m-x", label=("rho error"))
    plt.plot(xx, yy_K, "g-x", label=("K error"))
    plt.plot(xx, yy_G, "y-x", label=("G error"))
    plt.yscale("log")
    plt.xlabel("% Perovskite")
    plt.ylabel("Error")
    plt.legend()
    plt.show()
