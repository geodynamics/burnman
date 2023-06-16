# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import matplotlib.pyplot as plt


def pretty_plot():
    """
    Makes pretty plots. Overwrites the matplotlib default settings to allow
    for better fonts. Slows down plotting
    """
    import matplotlib.pyplot as plt

    plt.rc("text", usetex=True)
    plt.rcParams["text.latex.preamble"] = "\\usepackage{relsize}"
    plt.rc("font", family="sanserif")


def plot_projected_elastic_properties(
    mineral, plot_types, axes, n_zenith=31, n_azimuth=91, n_divs=100
):
    """
    :param mineral: Mineral object on which calculations should be done
    :type mineral: :class:`burnman.Mineral`

    :param plot_types: Plot types must be one of the following
        * 'vp' - V$_{P}$ (km/s)
        * 'vs1' - 'V$_{S1}$ (km/s)
        * 'vs2' - V$_{S2}$ (km/s)
        * 'vp/vs1' - V$_{P}$/V$_{S1}$
        * 'vp/vs2' - V$_{P}$/V$_{S2}$
        * 's anisotropy' - S-wave anisotropy (%), 200(vs1s - vs2s)/(vs1s + vs2s))
        * 'linear compressibility' - Linear compressibility (GPa$^{-1}$)
        * 'youngs modulus' - Youngs Modulus (GPa)
    :type plot_types: list of str

    :param axes: axes objects to be modified.
        Must be initialised with projection='polar'.
    :type axes: matplotlib.pyplot.axes objects

    :param n_zenith: Number of zeniths (plot resolution).
    :type n_zenith: int

    :param n_azimuth: Number of azimuths (plot resolution).
    :type n_azimuth: int

    :param n_divs: Number of divisions for the color scale.
    :type n_divs: int
    """

    assert len(plot_types) == len(axes)

    zeniths = np.linspace(np.pi / 2.0, np.pi, n_zenith)
    azimuths = np.linspace(0.0, 2.0 * np.pi, n_azimuth)
    Rs = np.sin(zeniths) / (1.0 - np.cos(zeniths))
    r, theta = np.meshgrid(Rs, azimuths)

    vps = np.empty_like(r)
    vs1s = np.empty_like(r)
    vs2s = np.empty_like(r)
    betas = np.empty_like(r)
    Es = np.empty_like(r)
    for i, az in enumerate(azimuths):
        for j, phi in enumerate(zeniths):
            d = np.array(
                [np.cos(az) * np.sin(phi), np.sin(az) * np.sin(phi), -np.cos(phi)]
            )  # change_hemispheres
            velocities = mineral.wave_velocities(d)
            betas[i][j] = mineral.isentropic_linear_compressibility(d)
            Es[i][j] = mineral.isentropic_youngs_modulus(d)
            vps[i][j] = velocities[0][0]
            vs1s[i][j] = velocities[0][1]
            vs2s[i][j] = velocities[0][2]

    prps = []
    for type in plot_types:
        if type == "vp":
            prps.append(("V$_{P}$ (km/s)", vps / 1000.0))
        elif type == "vs1":
            prps.append(("V$_{S1}$ (km/s)", vs1s / 1000.0))
        elif type == "vs2":
            prps.append(("V$_{S2}$ (km/s)", vs2s / 1000.0))
        elif type == "vp/vs1":
            prps.append(("V$_{P}$/V$_{S1}$", vps / vs1s))
        elif type == "vp/vs2":
            prps.append(("V$_{P}$/V$_{S2}$", vps / vs2s))
        elif type == "s anisotropy":
            prps.append(
                ("S-wave anisotropy (%)", 200.0 * (vs1s - vs2s) / (vs1s + vs2s))
            )
        elif type == "linear compressibility":
            prps.append(("Linear compressibility (GPa$^{-1}$)", betas * 1.0e9))
        elif type == "youngs modulus":
            prps.append(("Youngs Modulus (GPa)", Es / 1.0e9))
        else:
            raise Exception("plot_type not recognised.")

    contour_sets = []
    ticks = []
    lines = []
    for i, prp in enumerate(prps):
        title, item = prp

        axes[i].set_title(title)

        vmin = np.min(item)
        vmax = np.max(item)
        spacing = np.power(10.0, np.floor(np.log10(vmax - vmin)))
        nt = int((vmax - vmin - vmax % spacing + vmin % spacing) / spacing)
        if nt == 1:
            spacing = spacing / 4.0
        elif nt < 4:
            spacing = spacing / 2.0
        elif nt > 8:
            spacing = spacing * 2.0

        tmin = vmin + (spacing - vmin % spacing)
        tmax = vmax - vmax % spacing
        nt = int((tmax - tmin) / spacing + 1)

        ticks.append(np.linspace(tmin, tmax, nt))
        contour_sets.append(
            axes[i].contourf(
                theta, r, item, n_divs, cmap=plt.cm.jet_r, vmin=vmin, vmax=vmax
            )
        )
        lines.append(
            axes[i].contour(
                theta, r, item, ticks[-1], colors=("black",), linewidths=(1,)
            )
        )
        axes[i].set_yticks([100])

    return contour_sets, ticks, lines
