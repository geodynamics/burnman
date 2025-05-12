# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_perplex
---------------

This minimal example demonstrates how burnman can be used
to read and interrogate a Perple_X tab file, and how to smooth
thermodynamic properties on a given P-T grid..
For a workflow showing how to create Perple_X tables and save
smoothed properties along an adiabat in an ASPECT-compatible
format, see burnman/contrib/perplex/.

*Uses:*

* :doc:`PerplexMaterial`
* :func:`burnman.Material.evaluate`
* :func:`burnman.utils.math.smooth_array`

*Demonstrates:*

* Use of PerplexMaterial
* Smoothing gridded properties


"""
import numpy as np
import matplotlib.pyplot as plt

import burnman
from burnman.utils.math import smooth_array


if __name__ == "__main__":
    rock = burnman.PerplexMaterial("../burnman/data/input_perplex/in23_1.tab")
    P = 1.0e9
    T = 1650.0
    rock.set_state(P, T)
    print(
        "P: {0:.1f} GPa, T: {1:.1f} K, density: {2:.1f} kg/m^3".format(
            P / 1.0e9, T, rock.rho
        )
    )

    pressures = np.linspace(10.0e9, 25.0e9, 151)
    temperatures = [T] * len(pressures)
    densities = rock.evaluate(["rho"], pressures, temperatures)[0]
    plt.plot(pressures / 1.0e9, densities)
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Density (kg/m^3)")
    plt.show()

    pressures = np.linspace(10.0e9, 25.0e9, 151)
    temperatures = np.linspace(1600.0, 1800.0, 3)

    T = 1650.0
    entropies = rock.evaluate(["S"], pressures, np.array([T] * len(pressures)))[0]

    smoothed_entropies = smooth_array(
        array=entropies,
        grid_spacing=np.array([pressures[1] - pressures[0]]),
        gaussian_rms_widths=np.array([5.0e8]),
    )

    plt.plot(pressures / 1.0e9, entropies, label="entropies")
    plt.plot(pressures / 1.0e9, smoothed_entropies, label="smoothed entropies")
    plt.xlabel("Pressure (GPa)")
    plt.ylabel("Entropy (J/K/mol)")
    plt.legend(loc="upper right")
    plt.show()
