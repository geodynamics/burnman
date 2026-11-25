# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences.
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
Benchmark SSC, Voigt and Reuss moduli against published Echoes results
(https://jfbarthelemy.github.io/echoes/4_schemes/schemes.html).

The original online figure is an underlay, and its digitized SC curve and
printed numerical results are independent references.
"""

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

from burnman.classes.averaging_schemes import Reuss, SelfConsistent, Voigt

# Pixel calibration from the original figure's ticks, not calculated curves.
# Entries are (x at phi=0, pixels per unit phi, y at modulus=0, pixels/GPa).
CALIBRATION = [(145.0, 655.0, 609.0, 7.29), (997.5, 655.0, 609.0, 98.5 / 6.0)]
PANELS = [(113, 833, 58, 635), (966, 1685, 58, 635)]


def calculate(porosities, scheme=None):
    """Use the same public scalar interface as Composite/VoigtReussHill."""
    if scheme is None:
        scheme = SelfConsistent()
    bulk, shear = np.array([72.0, 1.0e-6]), np.array([32.0, 1.0e-6])
    return np.array(
        [
            [
                scheme.average_bulk_moduli([1.0 - phi, phi], bulk, shear),
                scheme.average_shear_moduli([1.0 - phi, phi], bulk, shear),
            ]
            for phi in porosities
        ]
    )


if __name__ == "__main__":
    figures = Path(__file__).resolve().parent / "figures"
    reference = np.array(
        [
            [0, 210, 218.0],
            [0, 243, 281.0],
            [0, 276, 339.0],
            [0, 309, 394.0],
            [0, 342, 445.0],
            [0, 374, 491.0],
            [0, 407, 535.0],
            [0, 440, 574.0],
            [1, 1063, 185.5],
            [1, 1096, 237.0],
            [1, 1128, 287.0],
            [1, 1161, 339.0],
            [1, 1194, 392.0],
            [1, 1227, 446.0],
            [1, 1260, 500.0],
            [1, 1292, 554.0],
        ]
    )
    for phi, expected, tolerance in [
        (0.2, [37.156, 19.629], 5.1e-4),
        (0.3, [22.689245, 13.264364], 2.0e-5),
    ]:
        actual = calculate([phi])[0]
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=tolerance)
        print(
            f"phi={phi:.1f}: K={actual[0]:.6f}, G={actual[1]:.6f} GPa; published values passed"
        )

    for panel, name in enumerate(["Bulk", "Shear"]):
        pixels = reference[reference[:, 0] == panel, 1:]
        x0, dx, y0, dy = CALIBRATION[panel]
        phi = (pixels[:, 0] - x0) / dx
        observed = (y0 - pixels[:, 1]) / dy
        actual = calculate(phi)[:, panel]
        error_pixels = np.abs(actual - observed) * dy
        np.testing.assert_allclose(actual, observed, rtol=0.0, atol=2.0 / dy)
        print(
            f"{name}: {len(phi)} digitized points; maximum error {max(error_pixels):.3f} pixels (limit 2)"
        )

    image = mpimg.imread(figures / "averaging_schemes_echoes_Fig_12_4.png")
    porosities = np.linspace(0.0, 1.0, 101)
    curves = [
        ("SSC", calculate(porosities), "#007c91", "--", "o", (0, 5)),
        ("Voigt", calculate(porosities, Voigt()), "#c55a00", "-.", "^", (5, 10)),
        ("Reuss", calculate(porosities, Reuss()), "#7954af", ":", "s", (2, 10)),
    ]
    # Independent closed-form bounds, Echoes equations 12.3 and 12.4.
    # Check the tiny Reuss values numerically: they look like zero at the
    # scale of the reference figure and are printed as 0.000 on the page.
    fractions = np.column_stack((1.0 - porosities, porosities))
    phase_moduli = np.array([[72.0, 32.0], [1.0e-6, 1.0e-6]])
    voigt_reference = fractions @ phase_moduli
    reuss_reference = 1.0 / (fractions @ (1.0 / phase_moduli))
    for name, actual, expected in [
        ("Voigt", curves[1][1], voigt_reference),
        ("Reuss", curves[2][1], reuss_reference),
    ]:
        np.testing.assert_allclose(actual, expected, rtol=1.0e-12, atol=0.0)
        print(f"{name}: 101 porosities; analytical bound passed")
    for name, scheme, expected in [
        ("Voigt", Voigt(), [57.600, 25.600]),
        ("Reuss", Reuss(), [0.000, 0.000]),
    ]:
        actual = calculate([0.2], scheme)[0]
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=5.1e-4)
        print(
            f"{name} phi=0.2: K={actual[0]:.6g}, G={actual[1]:.6g} GPa; published values passed"
        )
    np.testing.assert_array_less(curves[2][1], curves[0][1] + 1.0e-12)
    np.testing.assert_array_less(curves[0][1], curves[1][1] + 1.0e-12)
    print("SSC lies between Reuss and Voigt at all 101 porosities")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8))
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.23, top=0.86, wspace=0.20)
    for panel, (ax, modulus_name) in enumerate(zip(axes, ["Bulk", "Shear"])):
        x0, dx, y0, dy = CALIBRATION[panel]
        left, right, top, bottom = PANELS[panel]
        extent = [
            (left - x0) / dx,
            (right - x0) / dx,
            (y0 - bottom) / dy,
            (y0 - top) / dy,
        ]
        ax.imshow(
            image[top:bottom, left:right], extent=extent, aspect="auto", alpha=0.45
        )
        for name, moduli, color, linestyle, marker, markevery in curves:
            ax.plot(
                porosities,
                moduli[:, panel],
                color=color,
                linestyle=linestyle,
                linewidth=1.6,
                marker=marker,
                markerfacecolor="none",
                markersize=5,
                markevery=markevery,
                label=f"BurnMan {name}",
            )
        pixels = reference[reference[:, 0] == panel, 1:]
        ax.plot(
            (pixels[:, 0] - x0) / dx,
            (y0 - pixels[:, 1]) / dy,
            "+",
            color="black",
            markersize=7,
            label="Digitized Echoes SC",
        )
        ax.set_xlim(extent[:2])
        ax.set_ylim(extent[2:])
        ax.set_xlabel(r"Porosity $\varphi$", fontsize=13)
        ax.set_ylabel(f"{modulus_name} modulus (GPa)", fontsize=13)
        ax.tick_params(labelsize=11)
    fig.suptitle(
        "Elastic averaging benchmark: SSC, Voigt and Reuss", fontsize=17, y=0.98
    )
    fig.text(
        0.5,
        0.91,
        "Solid: K = 72, G = 32 GPa; pore phase: K = G = 1 kPa",
        ha="center",
        fontsize=12,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.085),
        ncol=4,
        frameon=False,
        fontsize=12,
    )
    fig.text(
        0.5,
        0.035,
        "Faded underlay: Echoes Fig. 12.4, Jean-François Barthélémy. CC BY-SA 4.0\nhttps://jfbarthelemy.github.io/echoes/4_schemes/schemes.html",
        ha="center",
        fontsize=9,
    )
    plt.show()
