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


if __name__ == "__main__":
    check_figures()
