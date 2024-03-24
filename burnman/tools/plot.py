# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R


def pretty_plot():
    """
    Makes pretty plots. Overwrites the matplotlib default settings to allow
    for better fonts. Slows down plotting
    """
    import matplotlib.pyplot as plt

    plt.rc("text", usetex=True)
    plt.rcParams["text.latex.preamble"] = "\\usepackage{relsize}"
    plt.rc("font", family="sanserif")


def _to_polar(vecs):
    """
    Convenience function converting x, y, z vectors to
    normalised positions on a polar plot.

    :param vecs: Array of vectors to convert
    :type vecs: numpy array
    :return: Arrays of azimuths and distances from the origin
    :rtype: Two numpy arrays
    """
    shape = vecs.shape
    vecs = vecs.reshape(shape[0] * shape[1], 3)
    Rs = np.linalg.norm(vecs, axis=-1)

    zeniths = np.arccos(-vecs[:, 2] / Rs)
    thetas = (np.arctan2(vecs[:, 1], vecs[:, 0])).reshape(shape[0], shape[1])
    r = (np.sin(zeniths) / (1.0 - np.cos(zeniths))).reshape(shape[0], shape[1])
    return thetas, r


def _positions(
    n_zenith,
    n_azimuth,
    min_zenith,
    max_zenith,
    latitude_direction=False,
    azimuth_midpoints=False,
):
    """
    Convenience function creating a set of equally spaced
    x, y, z vectors and theta and r values on a polar plot.

    :param n_zenith: Number of zeniths
    :type n_zenith: integer
    :param n_azimuth: Number of azimuths
    :type n_azimuth: integer
    :param min_zenith: _description_
    :type min_zenith: float
    :param max_zenith: _description_
    :type max_zenith: float
    :param latitude_direction: Whether to return both positions and vectors
        along local lines of latitude,
        defaults to False
    :type latitude_direction: bool, optional
    :param azimuth_midpoints: Whether to use azimuthal midpoints,
        rather than starting at 0 and ending at 2*pi. defaults to False
    :type azimuth_midpoints: bool, optional
    :return: Either a full basis or just the positional x, y, z vectors,
        azimuth values and r values for plotting on a polar plot.
    :rtype: Three or five numpy arrays
    """
    # Upper hemisphere plot, zenith starts at pi/2, increases to pi.
    zeniths = np.linspace(min_zenith, max_zenith, n_zenith)
    if azimuth_midpoints:
        azimuths = np.linspace(
            (np.pi / (n_azimuth + 1)),
            2.0 * np.pi - (np.pi / (n_azimuth + 1)),
            n_azimuth,
        )
    else:
        azimuths = np.linspace(0.0, 2.0 * np.pi, n_azimuth)

    Rs = np.sin(zeniths) / (1.0 - np.cos(zeniths))
    r, theta = np.meshgrid(Rs, azimuths)

    cos_az = np.cos(azimuths)
    sin_phi = np.sin(zeniths)
    sin_az = np.sin(azimuths)
    cos_phi = np.cos(zeniths)

    cos_phi, cos_az = np.meshgrid(cos_phi, cos_az)
    sin_phi, sin_az = np.meshgrid(sin_phi, sin_az)

    d = np.moveaxis(np.array([cos_az * sin_phi, sin_az * sin_phi, -cos_phi]), 0, 2)

    if latitude_direction:
        d1 = np.moveaxis(np.array([cos_az * cos_phi, sin_az * cos_phi, sin_phi]), 0, 2)
        return d, d1, theta, r
    else:
        return d, theta, r


def _rotated_vector_sets(rotation_axes, vector_initial, n_rots):
    """
    Convenience function to calculate a set of
    vectors rotated between 0 and 180 degrees
    from their initial positions about an axis.

    :param rotation_axes: rotation axes
    :type rotation_axes: numpy array
    :param vector_1_initial: vectors to be rotated
    :type vector_1_initial: numpy array
    :param vector_2_initial: another set of vectors to be rotated
    :type vector_2_initial: numpy array
    :param n_rots: number of rotations between 0 and 180 degrees.
    :type n_rots: integer
    :return: A new set of vectors
    :rtype: Two numpy arrays
    """
    # now rotate d1 and d2 about d in 1 degree increments
    n_azimuth = vector_initial.shape[0]
    n_zenith = vector_initial.shape[1]
    rotations = np.linspace(0.0, np.pi, n_rots)
    axial_directions = np.empty((n_azimuth, n_zenith, n_rots, 3))
    lateral_directions = np.empty((n_azimuth, n_zenith, n_rots, 3))
    for i in range(n_rots):
        axial_directions[:, :, i, :] = rotation_axes
    for i in range(n_azimuth):
        for j in range(n_zenith):
            rotvecs = np.einsum("i, l", rotations, rotation_axes[i, j])
            rot = R.from_rotvec(rotvecs)
            lateral_directions[i, j] = rot.apply(vector_initial[i, j])
    return axial_directions, lateral_directions


def _make_lines(positions, vectors, scales):
    """
    Convenience function to make line segments for polar plots.

    :param positions: Center position (x, y, z)
    :type positions: numpy array
    :param vectors: Vector away from the center position (also Cartesian)
    :type vectors: numpy array
    :param scales: Scaling fractions for the length of line.
    :type scales: numpy array
    :return: The requested line segments
    :rtype: numpy array
    """
    vec_start = positions - np.einsum("ijk, ij->ijk", vectors, scales)
    vec_end = positions + np.einsum("ijk, ij->ijk", vectors, scales)
    n_lines = positions.shape[0] * positions.shape[1]

    vlines = np.array([_to_polar(vec_start), _to_polar(vec_end)]).reshape(2, 2, n_lines)
    vlines = np.moveaxis(vlines, (0, 2), (2, 0))
    return vlines


def _make_poisson_lines(mineral, n_rots, vector_fraction, plot_minimum=False):
    """
    Convenience function to make line segments
    for the polar plot to indicate the orientations of
    the minimum or maximum Poisson ratio eigenvectors

    :param mineral: Anisotropic mineral to plot.
    :type mineral: AnisotropicMineral
    :param n_rots: Number of rotations to find the minimum
    :type n_rots: integer
    :param vector_fraction: Fractional length of plotted unit vector.
    :type vector_fraction: float
    :param plot_minimum: Whether to plot the minimum (or maximum),
        defaults to False
    :type plot_minimum: bool, optional

    :return: The requested line segments.
    :rtype: numpy array
    """
    v = []
    for i in range(3):
        n_zenith = 2
        n_azimuth = 36 - 16 * i  # 36 - 16*2 = 4
        min_zenith = (7.0 + 2 * i) * np.pi / 13.0
        max_zenith = (8.0 + 2 * i) * np.pi / 13.0
        d0, d01, _, r = _positions(
            n_zenith, n_azimuth, min_zenith, max_zenith, True, True
        )

        axial, lateral = _rotated_vector_sets(d0, d01, n_rots)
        poisson0 = mineral.isentropic_poissons_ratio(axial, lateral)

        if plot_minimum:
            inds0 = np.ix_(*[np.arange(i) for i in poisson0.shape[:-1]]) + (
                np.argmin(poisson0, axis=-1),
            )
            v.extend(
                list(_make_lines(d0, lateral[inds0], vector_fraction / np.sqrt(r)))
            )
        else:
            inds0 = np.ix_(*[np.arange(i) for i in poisson0.shape[:-1]]) + (
                np.argmax(poisson0, axis=-1),
            )
            v.extend(
                list(_make_lines(d0, lateral[inds0], vector_fraction / np.sqrt(r)))
            )
    return np.array(v)


def _make_velocity_lines(mineral, vector_fraction):
    """
    Convenience function to make line segments
    for the polar plot to indicate the orientations of
    the VS1 and VS2 eigenvectors

    :param mineral: Anisotropic mineral to plot.
    :type mineral: AnisotropicMineral
    :param vector_fraction: Fractional length of plotted unit vector.
    :type vector_fraction: float

    :return: The requested line segments.
    :rtype: Two 2D numpy arrays
    """
    v_vs1 = []
    v_vs2 = []
    for i in range(3):
        n_zenith = 2
        n_azimuth = 36 - 16 * i  # 36 - 16*2 = 4
        min_zenith = (7.0 + 2 * i) * np.pi / 13.0
        max_zenith = (8.0 + 2 * i) * np.pi / 13.0
        d0, _, r = _positions(n_zenith, n_azimuth, min_zenith, max_zenith, False, True)

        _, eigenvectors = mineral.wave_velocities(d0)

        v_vs1.extend(
            list(_make_lines(d0, eigenvectors[:, :, 1], vector_fraction / np.sqrt(r)))
        )
        v_vs2.extend(
            list(_make_lines(d0, eigenvectors[:, :, 2], vector_fraction / np.sqrt(r)))
        )
    return np.array(v_vs1), np.array(v_vs2)


def plot_projected_elastic_properties(
    mineral,
    plot_types,
    axes,
    n_zenith=31,
    n_azimuth=91,
    n_divs=100,
    n_rots=181,
    plot_vectors=True,
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
        * 'minimum poisson ratio' - Minimum Poisson ratio
        * 'maximum poisson ratio' - Maximum Poisson ratio
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

    :param n_rots: Number of along-vector rotations to
        find minimum and maximum Poisson ratios.
    :type n_rots: int

    :param plot_vectors: Whether or not to plot vectors
        for $V_S$ and Poisson axial directions
    :type plot_vectors: bool
    """

    assert len(plot_types) == len(axes)

    # Upper hemisphere plot, zenith starts at pi/2, increases to pi.
    d, d1, theta, r = _positions(n_zenith, n_azimuth, np.pi / 2.0, np.pi, True)

    velocities, _ = mineral.wave_velocities(d)
    vps, vs1s, vs2s = np.moveaxis(velocities, -1, 0)
    v_vs1s, v_vs2s = _make_velocity_lines(mineral, 0.02)

    if "minimum poisson ratio" in plot_types or "maximum poisson ratio" in plot_types:
        axial_directions, lateral_directions = _rotated_vector_sets(d, d1, n_rots)
        poisson = mineral.isentropic_poissons_ratio(
            axial_directions, lateral_directions
        )

    prps = []
    for type in plot_types:
        if type == "vp":
            prps.append(("V$_{P}$ (km/s)", vps / 1000.0))
        elif type == "vs1":
            prps.append(("V$_{S1}$ (km/s)", (vs1s / 1000.0, v_vs1s)))
        elif type == "vs2":
            prps.append(("V$_{S2}$ (km/s)", (vs2s / 1000.0, v_vs2s)))
        elif type == "vp/vs1":
            prps.append(("V$_{P}$/V$_{S1}$", vps / vs1s))
        elif type == "vp/vs2":
            prps.append(("V$_{P}$/V$_{S2}$", vps / vs2s))
        elif type == "s anisotropy":
            prps.append(
                ("S-wave anisotropy (%)", 200.0 * (vs1s - vs2s) / (vs1s + vs2s))
            )
        elif type == "linear compressibility":
            betas = mineral.isentropic_linear_compressibility(d)
            prps.append(("Linear compressibility (GPa$^{-1}$)", betas * 1.0e9))
        elif type == "youngs modulus":
            Es = mineral.isentropic_youngs_modulus(d)
            prps.append(("Youngs Modulus (GPa)", Es / 1.0e9))
        elif type == "minimum poisson ratio":
            inds = np.ix_(*[np.arange(i) for i in poisson.shape[:-1]]) + (
                np.argmin(poisson, axis=-1),
            )
            vlines = _make_poisson_lines(mineral, n_rots, 0.02, plot_minimum=True)
            min_poisson = poisson[inds]
            prps.append(("Minimum Poisson ratio", (min_poisson, vlines)))
        elif type == "maximum poisson ratio":
            inds = np.ix_(*[np.arange(i) for i in poisson.shape[:-1]]) + (
                np.argmax(poisson, axis=-1),
            )
            max_poisson = poisson[inds]
            vlines = _make_poisson_lines(mineral, n_rots, 0.02, plot_minimum=False)
            prps.append(("Maximum Poisson ratio", (max_poisson, vlines)))
        else:
            raise Exception("plot_type not recognised.")

    contour_sets = []
    ticks = []
    lines = []
    both_vs = "vs1" in plot_types and "vs2" in plot_types
    both_poisson = (
        "minimum poisson ratio" in plot_types and "maximum poisson ratio" in plot_types
    )
    if both_vs:
        vsmin = np.min(np.array([vs1s, vs2s])) / 1.0e3
        vsmax = np.max(np.array([vs1s, vs2s])) / 1.0e3
    if both_poisson:
        poissonmin = np.min(min_poisson)
        poissonmax = np.max(max_poisson)
    for i, prp in enumerate(prps):
        title, values = prp
        plot_lines = False
        if len(values) == 2:
            if plot_vectors is True:
                plot_lines = True
            vlines = values[1]
            values = values[0]

        axes[i].set_title(title)

        vmin = np.min(values)
        vmax = np.max(values)
        vmin2 = np.min(values)
        vmax2 = np.max(values)

        if (title == "V$_{S1}$ (km/s)" or title == "V$_{S2}$ (km/s)") and both_vs:
            vmin = vsmin
            vmax = vsmax
        if (
            title == "Minimum Poisson ratio" or title == "Maximum Poisson ratio"
        ) and both_poisson:
            vmin = poissonmin
            vmax = poissonmax

        spacing = np.power(10.0, np.floor(np.log10(vmax2 - vmin2)))
        nt = int((vmax2 - vmin2 - vmax2 % spacing + vmin2 % spacing) / spacing)
        if nt == 1:
            spacing = spacing / 4.0
        elif nt < 4:
            spacing = spacing / 2.0
        elif nt > 8:
            spacing = spacing * 2.0

        tmin = vmin2 + (spacing - vmin2 % spacing)
        tmax = vmax2 - vmax2 % spacing
        nt = int((tmax - tmin) / spacing + 1)

        ticks.append(np.linspace(tmin, tmax, nt))
        contour_sets.append(
            axes[i].contourf(
                theta, r, values, n_divs, cmap=plt.cm.jet_r, vmin=vmin, vmax=vmax
            )
        )
        lines.append(
            axes[i].contour(
                theta, r, values, ticks[-1], colors=("black",), linewidths=(1,)
            )
        )
        if plot_lines:
            for line in vlines:
                axes[i].plot(line[0], line[1], color="white", linewidth=1)
            axes[i].set_ylim(0.0, 1.0)

        axes[i].set_yticks([])

    return contour_sets, ticks, lines
