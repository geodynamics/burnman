# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


from __future__ import absolute_import

import numpy as np
import matplotlib.pyplot as plt
import pkgutil

from ..classes.planet import Planet
from ..classes.layer import Layer


def write_tvel_file(planet_or_layer, modelname='burnmanmodel',
                    background_model=None):
    """
    Function to write input file for obspy travel time calculations.
    Note: Because density isn't defined for most 1D seismic models, densities
    are output as zeroes.  The tvel format has a column for density,
    but this column is not used by obspy for travel time calculations.
    Parameters
    ----------
    planet_or_layer  :  burnman.Planet() or burnman.Layer()
        Planet or layer to write out to tvel file
    filename : string
        Filename to read to
    background_model : burnman.seismic.Seismic1DModel()
        1D seismic model to fill in parts of planet
        (likely to be an earth model) that aren't defined by layer
        (only need when using Layer())
    """

    if not isinstance(planet_or_layer, (Planet, Layer)):
        raise TypeError("Input must be a Planet() or Layer() ")

    if isinstance(planet_or_layer, Layer):
        assert(background_model)
        layer = planet_or_layer
        depths = background_model.internal_depth_list()

        above_layer = np.where(
            depths < (np.max(depths) - layer.outer_radius))[-1]
        below_layer = np.where(
            depths > (np.max(depths) - layer.inner_radius))[0]

        data_above = list(zip(depths[above_layer] / 1.e3,
                              background_model.v_p(depths[above_layer]) / 1.e3,
                              background_model.v_s(depths[above_layer]) / 1.e3,
                              np.zeros_like(depths[above_layer])))
        data_layer = list(zip((np.max(depths) - layer.radii)[::- 1] / 1.e3,
                              layer.v_p[::-1] / 1.e3,
                              layer.v_s[::-1] / 1.e3,
                              layer.density[::-1] / 1.e3))
        data_below = list(zip(depths[below_layer] / 1.e3,
                              background_model.v_p(depths[below_layer]) / 1.e3,
                              background_model.v_s(depths[below_layer]) / 1.e3,
                              np.zeros_like(depths[below_layer])))

        data = data_above + data_layer + data_below

        header = (f'{layer.name}  model from BurnMan between a radius of '
                  f'{str(layer.inner_radius)} and '
                  f'{str(layer.outer_radius)} km \n'
                  f'{background_model.__class__.__name__} '
                  f'for the rest of the Earth')
        with open(modelname + '.tvel', 'wb') as f:
            np.savetxt(f, data, header=header, fmt='%5.2f', delimiter='\t')

    if isinstance(planet_or_layer, Planet):
        planet = planet_or_layer
        data = list(zip((planet.radius_planet - planet.radii)[::-1] / 1.e3,
                        planet.v_p[::-1] / 1.e3,
                        planet.v_s[::-1] / 1.e3,
                        planet.density[::- 1] / 1.e3))

        header = (f'{planet.name} model from BurnMan with a radius of '
                  f'{str(planet.radius_planet)} km \n'
                  f'Layers of planet are '
                  f'{", ".join(layer.name for layer in planet.layers)}')
        with open(modelname + '.tvel', 'wb') as f:
            np.savetxt(f, data, header=header, fmt='%5.2f', delimiter='\t')


def write_axisem_input(layers, modelname='burnmanmodel_foraxisem',
                       axisem_ref='axisem_prem_ani_noocean.txt',
                       plotting=False):
    """
    Writing velocities and densities to AXISEM (www.axisem.info) input file.
    The input can be a single layer, or a list of layers taken from a planet (planet.layers).
    Currently this function will implement explicit discontinuities between layers in the seismic model.
    Currently this function is only set for Earth.


    Parameters
    ----------
    layers : list of one or more burnman.Layer()
        List of layers to put in axisem file
    modelname : string
        Name of model, appears in name of output file
    axisem_ref : string
        Reference file, used to copy the header and for the rest of the planet, in the case of a Layer(), default = 'axisem_prem_ani_noocean.txt'
    plotting: Boolean
        True means plot of the old model and replaced model will be shown (default = False)
    """

    if not isinstance(layers[0], Layer):
        raise TypeError("Input must be a list of Layer()")
    # Load reference input
    datastream = pkgutil.get_data(
        'burnman', 'data/input_seismic/' + axisem_ref)
    lines = [line.strip()
             for line in datastream.decode('ascii').split('\n') if line.strip()]
    table = []
    for line in lines[18:]:
        numbers = np.fromstring(line, sep=' ')
        if len(numbers) > 0:
            if line[0] != "#" and line[0] != "%":
                table.append(numbers)
    table = np.array(table)
    # format is
    # radius density vpv vsv Qk Qmu vph vsh eta

    if plotting:
        plt.figure(figsize=(12, 6))
        plt.plot(table[:, 0] / 1.e3, table[:, 2] / 1.e3,
                 color='g', linestyle='--')
        plt.plot(table[:, 0] / 1.e3, table[:, 3] / 1.e3,
                 color='b', linestyle='--')
        plt.plot(table[:, 0] / 1.e3, table[:, 1] / 1.e3,
                 color='r', linestyle='--')

    for layer in layers:
        # Cutting out range to input in Axisem reference file, and adding in
        # model values at the top and bottom of Layer.
        i_min = next(x[0] for x in enumerate(table[:, 0])
                     if x[1] <= layer.outer_radius)
        if table[i_min, 0] - layer.outer_radius < 0:
            table = np.insert(table, i_min, None, axis=0)
            table[i_min, 0] = layer.outer_radius

        i_max = next(x[0] for x in enumerate(table[:, 0])
                     if x[1] <= layer.inner_radius)

        if table[i_max, 0] - layer.inner_radius < 0:
            table = np.insert(table, i_max, None, axis=0)
            table[i_max, 0] = layer.inner_radius

        layer_vp, layer_vs, layer_density = layer.evaluate(
            ['v_p', 'v_s', 'density'], radlist=table[i_min:i_max, 0], radius_planet=np.max(table[:, 0]))

        table[i_min:i_max, 2] = layer_vp
        table[i_min:i_max, 6] = layer_vp
        table[i_min:i_max, 3] = layer_vs
        table[i_min:i_max, 7] = layer_vs
        table[i_min:i_max, 1] = layer_density

    # WRITE OUT FILE
    filename = 'axisem_' + modelname + '.txt'
    f = open(filename, 'w')
    print('Writing ' + filename + ' ...')
    f.write(f'# Input file {modelname} for AXISEM created using BurnMan, '
            f'replacing {axisem_ref} between '
            f'{str(np.round(layer.inner_radius / 1.e3))} and '
            f'{str(np.round(layer.outer_radius / 1.e3))} km\n')

    discontinuity = 0  # Number discontinuities
    f.write('NAME ' + modelname + '\n')
    for line in lines[2:18]:
        f.write(line[:] + '\n')
    for i in range(len(table[:, 0])):
        if i > 0 and table[i, 0] == table[i - 1, 0]:
            discontinuity = discontinuity + 1
            f.write(f'#          Discontinuity   {str(discontinuity)}, '
                    f'depth:    {str(np.round((6371.e3 - table[i, 0]) / 1.e3, decimals=2))} km \n')

        f.write(f'{table[i, 0]:8.0f} {table[i, 1]:9.2f} {table[i, 2]:9.2f} '
                f'{table[i, 3]:9.2f} {table[i, 4]:9.1f} {table[i, 5]:9.1f} '
                f'{table[i, 6]:9.2f} {table[i, 7]:9.2f} {table[i, 8]:9.5f} \n')

    f.close()

    if plotting:

        plt.plot(table[:, 0] / 1.e3, table[:, 2] / 1.e3,
                 color='g', linestyle='-', label='V_p')
        plt.plot(table[:, 0] / 1.e3, table[:, 3] / 1.e3,
                 color='b', linestyle='-', label='V_s')
        plt.plot(table[:, 0] / 1.e3, table[:, 1] / 1.e3,
                 color='r', linestyle='-', label='density')

        plt.title(f'{filename} = {axisem_ref} replaced between '
                  f'{str(layer.inner_radius / 1.e3)} and '
                  f'{str(layer.outer_radius / 1.e3)} km')
        plt.legend(loc='lower right')
        plt.show()


def write_mineos_input(layers, modelname='burnmanmodel_formineos',
                       mineos_ref='mineos_prem_noocean.txt', plotting=False):
    """
    Writing velocities and densities to
    Mineos (https://geodynamics.org/cig/software/mineos/) input file
    Note:
        - Currently, this function only honors the discontinuities already
        in the synthetic input file, so it is best to only replace
        certain layers with burnman values

    Parameters
    ----------
    layers : list of one or more burnman.Layer()
        List of layers to put in axisem file
    modelname : string
        Name of model, appears in name of output file
    mineos_ref : string
        Reference file, used to copy the header and for the rest of the planet,
        in the case of a Layer(), default = 'mineos_prem_noocean.txt'
    plotting: Boolean
        True means plot of the old model and replaced model will be shown
        (default = False)

    """

    if not isinstance(layers[0], Layer):
        raise TypeError("Input must be a list of Layer()")

    # Load reference input
    datastream = pkgutil.get_data(
        'burnman', 'data/input_seismic/' + mineos_ref)
    lines = [line.strip()
             for line in datastream.decode('ascii').split('\n') if line.strip()]
    table = []
    for line in lines[3:]:
        numbers = np.fromstring(line, sep=' ')
        table.append(numbers)
    table = np.array(table)

    if plotting:
        plt.figure(figsize=(12, 6))
        plt.plot(table[:, 0] / 1.e3, table[:, 2] / 1.e3,
                 color='g', linestyle='--')
        plt.plot(table[:, 0] / 1.e3, table[:, 3] / 1.e3,
                 color='b', linestyle='--')
        plt.plot(table[:, 0] / 1.e3, table[:, 1] / 1.e3,
                 color='r', linestyle='--')

    for layer in layers:
        i_min = next(x[0] for x in enumerate(table[:, 0])
                     if x[1] >= layer.inner_radius)
        if table[i_min, 0] - layer.inner_radius > 0:
            table[i_min, 0] = layer.inner_radius

        i_max = next(x[0] for x in enumerate(table[:, 0])
                     if x[1] >= layer.outer_radius)
        if table[i_max, 0] - layer.outer_radius > 0:
            table[i_max, 0] = layer.outer_radius

        layer_vp, layer_vs, layer_density = layer.evaluate(
            ['v_p', 'v_s', 'density'], radlist=table[i_min:i_max, 0], radius_planet=np.max(table[:, 0]))

        table[i_min:i_max, 2] = layer_vp
        table[i_min:i_max, 6] = layer_vp
        table[i_min:i_max, 3] = layer_vs
        table[i_min:i_max, 7] = layer_vs
        table[i_min:i_max, 1] = layer_density

    # WRITE OUT FILE
    filename = 'mineos_' + modelname + '.txt'
    f = open(filename, 'w')
    print('Writing ' + filename + ' ...')
    f.write(lines[0][:-2] + ' +  ' + filename + '\n')
    for line in lines[1:3]:
        f.write(line[:-2] + '\n')
    for i in range(len(table[:, 0])):
        f.write('%8.0f %9.2f %9.2f %9.2f %9.1f %9.1f %9.2f %9.2f %9.5f \n' % (
            table[i, 0], table[i, 1], table[i, 2], table[i, 3], table[i, 4], table[i, 5], table[i, 6], table[i, 7], table[i, 8]))

    f.close()

    if plotting:

        plt.plot(table[:, 0] / 1.e3, table[:, 2] / 1.e3,
                 color='g', linestyle='-', label='V_p')
        plt.plot(table[:, 0] / 1.e3, table[:, 3] / 1.e3,
                 color='b', linestyle='-', label='V_s')
        plt.plot(table[:, 0] / 1.e3, table[:, 1] / 1.e3,
                 color='r', linestyle='-', label='density')

        plt.title(f'{filename} = {mineos_ref} replaced between '
                  f'{str(layer.inner_radius / 1.e3)} and '
                  f'{str(layer.outer_radius / 1.e3)} km')
        plt.legend(loc='lower right')
        plt.show()
