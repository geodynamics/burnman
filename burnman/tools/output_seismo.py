# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


import os
import numpy as np
import matplotlib.pyplot as plt
import pkgutil

from ..classes.planet import Planet
from ..classes.layer import Layer


def write_tvel_file(planet_or_layer, modelname="burnmanmodel", background_model=None):
    """
    Function to write input file for obspy travel time calculations.
    Helper function which uses tvel_formatted_data_and_header()

    :param planet_or_layer: Planet or layer to write out to tvel file
    :type planet_or_layer: :class:`burnman.Planet` or :class:`burnman.Layer`.
    :param modelname: Basename of file to write to (suffix .tvel added by function).
    :type modelname: str
    :param background_model:  1D seismic model to fill in parts of planet
        (likely to be an earth model) that aren't defined by layer
        (only need when using :class:`burnman.Layer`)
    :type background_model: :class:`burnman.seismic.Seismic1DModel`
    """
    data, header = tvel_formatted_data_and_header(planet_or_layer, background_model)

    filename = os.devnull
    if modelname:
        filename = modelname + ".tvel"

    with open(filename, "wb") as f:
        np.savetxt(f, data, header=header, fmt="%5.2f", delimiter="\t")


def write_axisem_input(
    layers,
    modelname="burnmanmodel_foraxisem",
    axisem_ref="axisem_prem_ani_noocean.txt",
    plotting=False,
    verbose=True,
):
    """
    Function to write input file for AXISEM (www.axisem.info).
    Helper function which uses axisem_formatted_data_and_reference()

    The input can be a single layer, or a list of layers taken from a planet
    (planet.layers).
    Currently this function will implement explicit discontinuities between
    layers in the seismic model.
    Currently this function is only set for Earth.

    :param layers: List of layers to put in AXISEM file.
    :type layers: list of one or more :class:`burnman.Layer`
    :param modelname: Name of model, appears in name of output file.
    :type modelname: str
    :param axisem_ref: Reference file, used to copy the header
        and for the rest of the planet, in the case of a :class:`burnman.Layer`.
    :type axisem_ref: str
    :param plotting: Choose whether to show plot of the old model and replaced model.
    :type plotting: bool
    """

    table, lines = axisem_formatted_data_and_reference(layers, axisem_ref, plotting)

    filename = os.devnull
    if modelname:
        filename = "axisem_" + modelname + ".txt"
    else:
        modelname = ""

    if verbose:
        print("Writing " + filename + " ...")

    with open(filename, "w") as f:
        f.write(
            f"# Input file {modelname} for AXISEM created using BurnMan, "
            f"replacing {axisem_ref} between\n"
        )
        for layer in layers:
            f.write(
                f"# {str(np.round(layer.inner_radius / 1.e3))} and "
                f"{str(np.round(layer.outer_radius / 1.e3))} km\n"
            )

        formats = [
            "8.0f",
            "9.2f",
            "9.2f",
            "9.2f",
            "9.1f",
            "9.1f",
            "9.2f",
            "9.2f",
            "9.5f",
        ]
        discontinuity = 0  # Number discontinuities
        f.write("NAME " + modelname + "\n")
        for line in lines[2:18]:
            f.write(line[:] + "\n")
        for i in range(len(table[:, 0])):
            if i > 0 and table[i, 0] == table[i - 1, 0]:
                discontinuity = discontinuity + 1
                f.write(
                    f"#          Discontinuity   {str(discontinuity)}, "
                    f"depth:    {str(np.round((6371.e3 - table[i, 0]) / 1.e3, decimals=2))}"
                    " km \n"
                )

            f.write(
                " ".join(f"{val:{fmt}}" for val, fmt in zip(table[i, :9], formats))
                + "\n"
            )


def write_mineos_input(
    layers,
    modelname="burnmanmodel_formineos",
    mineos_ref="mineos_prem_noocean.txt",
    plotting=False,
    verbose=True,
):
    """
    Function to write input file for Mineos (https://geodynamics.org/cig/software/mineos/).
    Helper function which uses mineos_formatted_data_and_reference()

    Note: currently, this function only honors the discontinuities already
    in the synthetic input file, so it is best to only replace
    certain layers with burnman values

    :param layers: List of layers to put in Mineos file.
    :type layers: list of one or more :class:`burnman.Layer`
    :param modelname: Name of model, appears in name of output file.
    :type modelname: str
    :param mineos_ref: Reference file, used to copy the header
        and for the rest of the planet, in the case of a :class:`burnman.Layer`.
    :type mineos_ref: str
    :param plotting: Choose whether to show plot of the old model and replaced model.
    :type plotting: bool
    """

    table, lines = mineos_formatted_data_and_reference(layers, mineos_ref, plotting)

    formats = ["8.0f", "9.2f", "9.2f", "9.2f", "9.1f", "9.1f", "9.2f", "9.2f", "9.5f"]
    rows = (
        " ".join(f"{val:{fmt}}" for val, fmt in zip(row, formats)) + "\n"
        for row in table[:, :9]
    )

    filename = os.devnull
    if modelname:
        filename = "mineos_" + modelname + ".txt"
    else:
        modelname = ""

    if verbose:
        print("Writing " + filename + " ...")

    with open(filename, "w") as f:
        f.write(lines[0][:-2] + " +  " + filename + "\n")
        for line in lines[1:3]:
            f.write(line[:-2] + "\n")
            f.writelines(rows)


def tvel_formatted_data_and_header(planet_or_layer, background_model=None):
    """
    Formats data and creates a header for obspy travel time calculations.
    Note: Because density isn't defined for most 1D seismic models, densities
    are output as zeroes.  The tvel format has a column for density,
    but this column is not used by obspy for travel time calculations.

    :param planet_or_layer: Planet or layer to write out to tvel file
    :type planet_or_layer: :class:`burnman.Planet` or :class:`burnman.Layer`.
    :param background_model:  1D seismic model to fill in parts of planet
        (likely to be an earth model) that aren't defined by layer
        (only need when using :class:`burnman.Layer`)
    :type background_model: :class:`burnman.seismic.Seismic1DModel`
    """

    if not isinstance(planet_or_layer, (Planet, Layer)):
        raise TypeError("Input must be a Planet() or Layer() object.")

    if isinstance(planet_or_layer, Layer):
        assert background_model
        layer = planet_or_layer
        depths = background_model.internal_depth_list()

        above_layer = np.where(depths < (np.max(depths) - layer.outer_radius))[-1]
        below_layer = np.where(depths > (np.max(depths) - layer.inner_radius))[0]

        data_above = list(
            zip(
                depths[above_layer] / 1.0e3,
                background_model.v_p(depths[above_layer]) / 1.0e3,
                background_model.v_s(depths[above_layer]) / 1.0e3,
                np.zeros_like(depths[above_layer]),
            )
        )
        data_layer = list(
            zip(
                (np.max(depths) - layer.radii)[::-1] / 1.0e3,
                layer.v_p[::-1] / 1.0e3,
                layer.v_s[::-1] / 1.0e3,
                layer.density[::-1] / 1.0e3,
            )
        )
        data_below = list(
            zip(
                depths[below_layer] / 1.0e3,
                background_model.v_p(depths[below_layer]) / 1.0e3,
                background_model.v_s(depths[below_layer]) / 1.0e3,
                np.zeros_like(depths[below_layer]),
            )
        )

        data = data_above + data_layer + data_below

        header = (
            f"{layer.name}  model from BurnMan between a radius of "
            f"{str(layer.inner_radius)} and "
            f"{str(layer.outer_radius)} km \n"
            f"{background_model.__class__.__name__} "
            f"for the rest of the Earth"
        )

        return data, header

    if isinstance(planet_or_layer, Planet):
        planet = planet_or_layer
        data = list(
            zip(
                (planet.radius_planet - planet.radii)[::-1] / 1.0e3,
                planet.v_p[::-1] / 1.0e3,
                planet.v_s[::-1] / 1.0e3,
                planet.density[::-1] / 1.0e3,
            )
        )

        header = (
            f"{planet.name} model from BurnMan with a radius of "
            f"{str(planet.radius_planet)} km \n"
            f"Layers of planet are "
            f'{", ".join(layer.name for layer in planet.layers)}'
        )

        return data, header


def axisem_formatted_data_and_reference(
    layers,
    axisem_ref="axisem_prem_ani_noocean.txt",
    plotting=False,
):
    """
    Formats data for AXISEM (www.axisem.info).
    The input can be a single layer, or a list of layers taken from a planet
    (planet.layers).
    Currently this function will implement explicit discontinuities between
    layers in the seismic model.
    Currently this function is only set for Earth.

    :param layers: List of layers to put in AXISEM file.
    :type layers: list of one or more :class:`burnman.Layer`
    :param axisem_ref: Reference file, used to copy the header
        and for the rest of the planet, in the case of a :class:`burnman.Layer`.
    :type axisem_ref: str
    :param plotting: Choose whether to show plot of the old model and replaced model.
    :type plotting: bool
    """

    if not isinstance(layers[0], Layer):
        raise TypeError("Input must be a list of Layer()")
    # Load reference input
    datastream = pkgutil.get_data("burnman", "data/input_seismic/" + axisem_ref)
    reference_data = [
        line.strip() for line in datastream.decode("ascii").split("\n") if line.strip()
    ]
    table = []
    for line in reference_data[18:]:
        line = line.strip()
        if not line or line[0] in "#%":
            continue  # skip empty or comment lines
        numbers = np.array(line.split(), dtype=float)
        if numbers.size > 0:
            table.append(numbers)
    table = np.array(table)
    # format is
    # radius density vpv vsv Qk Qmu vph vsh eta

    if plotting:
        plt.figure(figsize=(12, 6))
        plt.plot(table[:, 0] / 1.0e3, table[:, 2] / 1.0e3, color="g", linestyle="--")
        plt.plot(table[:, 0] / 1.0e3, table[:, 3] / 1.0e3, color="b", linestyle="--")
        plt.plot(table[:, 0] / 1.0e3, table[:, 1] / 1.0e3, color="r", linestyle="--")

    for layer in layers:
        # Cutting out range to input in Axisem reference file, and adding in
        # model values at the top and bottom of Layer.
        i_min = next(x[0] for x in enumerate(table[:, 0]) if x[1] <= layer.outer_radius)
        if table[i_min, 0] - layer.outer_radius < 0:
            table = np.insert(table, i_min, None, axis=0)
            table[i_min, 0] = layer.outer_radius

        i_max = next(x[0] for x in enumerate(table[:, 0]) if x[1] <= layer.inner_radius)

        if table[i_max, 0] - layer.inner_radius < 0:
            table = np.insert(table, i_max, None, axis=0)
            table[i_max, 0] = layer.inner_radius

        lvp, lvs, lrho = layer.evaluate(
            ["v_p", "v_s", "density"],
            radlist=table[i_min:i_max, 0],
            radius_planet=np.max(table[:, 0]),
        )

        table[i_min:i_max, 2] = lvp
        table[i_min:i_max, 6] = lvp
        table[i_min:i_max, 3] = lvs
        table[i_min:i_max, 7] = lvs
        table[i_min:i_max, 1] = lrho

    if plotting:
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 2] / 1.0e3,
            color="g",
            linestyle="-",
            label="V_p",
        )
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 3] / 1.0e3,
            color="b",
            linestyle="-",
            label="V_s",
        )
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 1] / 1.0e3,
            color="r",
            linestyle="-",
            label="density",
        )

        plt.title(
            f"{axisem_ref} replaced between "
            f"{str(layer.inner_radius / 1.e3)} and "
            f"{str(layer.outer_radius / 1.e3)} km"
        )
        plt.legend(loc="lower right")
        plt.show()

    return table, reference_data


def mineos_formatted_data_and_reference(
    layers,
    mineos_ref="mineos_prem_noocean.txt",
    plotting=False,
):
    """
    Formats data and reference data for
    Mineos (https://geodynamics.org/cig/software/mineos/).

    Note: currently, this function only honors the discontinuities already
    in the synthetic input file, so it is best to only replace
    certain layers with burnman values

    :param layers: List of layers to put in Mineos file.
    :type layers: list of one or more :class:`burnman.Layer`
    :param mineos_ref: Reference file, used to copy the header
        and for the rest of the planet, in the case of a :class:`burnman.Layer`.
    :type mineos_ref: str
    :param plotting: Choose whether to show plot of the old model and replaced model.
    :type plotting: bool
    """

    if not isinstance(layers[0], Layer):
        raise TypeError("Input must be a list of Layer()")

    # Load reference input
    datastream = pkgutil.get_data("burnman", "data/input_seismic/" + mineos_ref)
    reference_data = [
        line.strip() for line in datastream.decode("ascii").split("\n") if line.strip()
    ]
    table = []
    for line in reference_data[3:]:
        numbers = np.fromstring(line, sep=" ")
        table.append(numbers)
    table = np.array(table)

    if plotting:
        plt.figure(figsize=(12, 6))
        plt.plot(table[:, 0] / 1.0e3, table[:, 2] / 1.0e3, color="g", linestyle="--")
        plt.plot(table[:, 0] / 1.0e3, table[:, 3] / 1.0e3, color="b", linestyle="--")
        plt.plot(table[:, 0] / 1.0e3, table[:, 1] / 1.0e3, color="r", linestyle="--")

    for layer in layers:
        i_min = next(x[0] for x in enumerate(table[:, 0]) if x[1] >= layer.inner_radius)
        if table[i_min, 0] - layer.inner_radius > 0:
            table[i_min, 0] = layer.inner_radius

        i_max = next(x[0] for x in enumerate(table[:, 0]) if x[1] >= layer.outer_radius)
        if table[i_max, 0] - layer.outer_radius > 0:
            table[i_max, 0] = layer.outer_radius

        lvp, lvs, lrho = layer.evaluate(
            ["v_p", "v_s", "density"],
            radlist=table[i_min:i_max, 0],
            radius_planet=np.max(table[:, 0]),
        )

        table[i_min:i_max, 2] = lvp
        table[i_min:i_max, 6] = lvp
        table[i_min:i_max, 3] = lvs
        table[i_min:i_max, 7] = lvs
        table[i_min:i_max, 1] = lrho

    if plotting:
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 2] / 1.0e3,
            color="g",
            linestyle="-",
            label="V_p",
        )
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 3] / 1.0e3,
            color="b",
            linestyle="-",
            label="V_s",
        )
        plt.plot(
            table[:, 0] / 1.0e3,
            table[:, 1] / 1.0e3,
            color="r",
            linestyle="-",
            label="density",
        )

        plt.title(
            f"{mineos_ref} replaced between "
            f"{str(layer.inner_radius / 1.e3)} and "
            f"{str(layer.outer_radius / 1.e3)} km"
        )
        plt.legend(loc="lower right")
        plt.show()

    return table, reference_data
