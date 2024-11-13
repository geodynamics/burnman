from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from burnman import calibrants
from burnman.tools.unitcell import molar_volume_from_unit_cell_volume


def make_VPT_figure(calibrant, temperatures, figure, figure_extent, plot_extent):
    name = str(calibrant).split(" ")[0][1:]
    print(f"Checking {name}...")

    fig = plt.figure(figsize=(6, 4))
    fig.suptitle(f"{name}")

    ax = [fig.add_subplot(1, 1, 1)]

    fig1 = mpimg.imread(figure)
    ax[0].imshow(fig1, extent=figure_extent, aspect="auto")

    pressures = np.linspace(plot_extent[0] * 1.0e9, plot_extent[1] * 1.0e9, 101)
    volumes = np.empty_like(pressures)
    for T in temperatures:
        for i, P in enumerate(pressures):
            volumes[i] = calibrant.volume(P, T) / molar_volume_from_unit_cell_volume(
                1.0, calibrant.params["Z"]
            )

        plt.plot(pressures / 1.0e9, volumes)

    ax[0].set_xlim(plot_extent[0], plot_extent[1])
    ax[0].set_ylim(plot_extent[2], plot_extent[3])
    plt.show()


def check_figures():
    make_VPT_figure(
        calibrants.Fei_2007.Au(),
        [300.0, 1473.0, 2173.0],
        "figures/Fei_2007_Au.png",
        [0, 139.5, 50, 68],
        [0, 139.5, 50, 68],
    )

    make_VPT_figure(
        calibrants.Fei_2007.Pt(),
        [300.0, 1473.0, 1873.0],
        "figures/Fei_2007_Pt.png",
        [-2, 125, 46.98, 61.02],
        [0, 125, 47, 61],
    )

    V_0 = calibrants.Holmes_1989.Pt().params[
        "V_0"
    ] / molar_volume_from_unit_cell_volume(1.0, 4)
    make_VPT_figure(
        calibrants.Holmes_1989.Pt(),
        [300.0],
        "figures/Holmes_1989_Pt.png",
        [0, 600, 0.6 * V_0, 1.1 * V_0],
        [0, 600, 0.6 * V_0, 1.1 * V_0],
    )


def check_Decker_1971():
    NaCl = calibrants.Decker_1971.NaCl_B1()
    V_0 = NaCl.params["V_0"]

    # First column in table
    V = V_0

    print(
        "Pressures from Decker 1971 calibrant "
        "vs. tabulated data in original paper (given in GPa)"
    )
    print(f"\nV={V:.4e} m^3/mol (standard state volume):")
    T_C = [25.0, 100.0, 200.0, 300.0, 500.0, 800.0]
    P_kbar = [0.00, 2.13, 5.00, 7.89, 13.72, 22.48]
    for i, T in enumerate(T_C):
        print(f"{T} C: {NaCl.pressure(V, T+273.15)/1.e9:.2f}, {P_kbar[i]/10.:.2f}")

    # Middle column in table
    V = (1.0 - 0.1904) * V_0
    print(f"\nV={V:.4e} m^3/mol (middle row):")
    T_C = [0.0, 25.0, 100.0, 200.0, 300.0, 500.0, 800.0]
    P_kbar = [83.24, 83.93, 86.02, 88.89, 91.79, 97.65, 106.52]
    for i, T in enumerate(T_C):
        print(f"{T} C: {NaCl.pressure(V, T+273.15)/1.e9:.2f}, {P_kbar[i]/10.:.2f}")

    # Last column in table
    V = (1.0 - 0.2950) * V_0
    print(f"\nV={V:.4e} m^3/mol (last row):")
    P_kbar = [193.45, 194.12, 196.18, 199.03, 201.93, 207.81, 216.73]
    for i, T in enumerate(T_C):
        print(f"{T} C: {NaCl.pressure(V, T+273.15)/1.e9:.2f}, {P_kbar[i]/10.:.2f}")
    print("")


if __name__ == "__main__":
    check_Decker_1971()
    check_figures()
