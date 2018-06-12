# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.
"""
This example script produces an input files for travel time calculations in ObsPy, either by replacing one layer with a composition defined in BurnMan (see example_layer) or with an entire planet (see example_planet). This example shows how to plot predicted travel times and ray paths and plots those in a PREM earth for reference.
Requires Obspy, see www.obspy.org and
L. Krischer, T. Megies, R. Barsch, M. Beyreuther, T. Lecocq, C. Caudron, J. Wassermann (2015)
    ObsPy: a bridge for seismology into the scientific Python ecosystem
    Computational Science & Discovery, 8(1), 014003
    DOI: 10.1088/1749-4699/8/1/014003

To find out more about the specific routines in this example see
    https://docs.obspy.org/packages/obspy.taup.html
"""
# Import supporting libraries
# Imports to be compatible with Python2 and Python3
from __future__ import absolute_import
from __future__ import print_function
import os
import sys  # Library used to interact with your operating system
import numpy as np  # Library used for general array
import matplotlib.pyplot as plt  # Library used for plotting
# Import BurnMan
sys.path.insert(1, os.path.abspath('../'))  # add path to burnman
import burnman
from burnman import minerals  # import mineral library seperately


# This example relies heavily  on the ObsPy, a python seismology toolkit
import obspy
from obspy.taup import TauPyModel




def plot_rays_and_times(modelname):
    """
    Calls obspy routines to plot ray paths and travel times for the model built in BurnMan and for a seismic reference model for comparison. 
    
    Parameters
    ----------
    modelname:
    Name for BurnMan model (*.tvel file must be present)
    """
    
    
    # Arrivals to plot, some random examples of phase names to use ("P","S", "PcP", "Sdiff", "SKS", "PKIKP")
    # Phase naming in obspy.taup is explained at
    # https://docs.obspy.org/packages/obspy.taup.html
    phase_list = ["P", "PKP", "PKIKP"]
    source_depth = 10 # in km
    min_degrees = 60 # minimum distance for ray paths
    max_degrees = 300 # maximum distance for ray paths
    npoints = 9 # number of distances to plot ray paths
    ref_model = 'prem' # choice of models available in obpsy, or input an npz file name
    
    # Build a taup_model for Obspy
    obspy.taup.taup_create.build_taup_model(
                "./" + modelname + ".tvel", ".")

    # Time to plot some predictions using routines from Obspy
    plt.figure(figsize=[9, 7])
    ax = plt.subplot(2, 2, 1)
    # plotting predicted travel times at all distances
    obspy.taup.plot_travel_times(
        ax=ax,
        model='./' +
        modelname +
        '.npz',
        source_depth=source_depth,
        phase_list=phase_list,
        show=False)
    plt.title(modelname)
    # plotting the same for PREM for reference
    ax = plt.subplot(2, 2, 2)
    obspy.taup.plot_travel_times(
        ax=ax,
        model=ref_model,
        source_depth=source_depth,
        phase_list=phase_list,
        show=False)
    # not sure why the grid dissapears on this subplot, reactivate here...
    ax.grid()
    plt.title(ref_model)
    # plotting predicted ray paths every 30 degrees between 60 and 300
    # degrees
    ax = plt.subplot(2, 2, 3, polar=True)
    obspy.taup.plot_ray_paths(
        ax=ax,
        model='./' +
        modelname +
        '.npz',
        source_depth=source_depth,
        min_degrees=min_degrees,
        max_degrees=max_degrees,
        npoints=npoints,
        phase_list=phase_list,
        verbose=True,
        show=False)
    # plotting the same for PREM for reference
    ax = plt.subplot(2, 2, 4, polar=True)
    obspy.taup.plot_ray_paths(
        ax=ax,
        model=ref_model,
        source_depth=source_depth,
        min_degrees=min_degrees,
        max_degrees=max_degrees,
        npoints=npoints,
        phase_list=phase_list,
        verbose=True)



if __name__ == "__main__":
    # Two examples available
    example_layer = True
    example_planet = True


    # First example: replacing the lower mantle with a composition from BurnMan
    if example_layer:
        modelname = 'perovskitic_mantle'
        # This is the first actual work done in this example.  We define
        # composite object and name it "rock".
        mg_fe_perovskite = minerals.SLB_2011.mg_fe_perovskite()
        mg_fe_perovskite.set_composition(
            [0.9, 0.1, 0])  # frac_mg, frac_fe, frac_al
        rock = burnman.Composite([mg_fe_perovskite], [1.])

        # We create an array of 20 depths at which we want to evaluate the
        # layer at
        depths = np.linspace(2890e3, 670e3, 20)
        # Here we define the lower mantle as a Layer(). The layer needs various
        # parameters to set a depth array and radius array.
        lower_mantle = burnman.Layer(
            name='Perovskitic Lower Mantle',
            radii=6371.e3 - depths)
        # Here we set the composition of the layer as the above defined 'rock'.
        lower_mantle.set_material(rock)

        # Now we set the temperature mode of the layer.
        # Here we use an adiabatic temperature and set the temperature at the
        # top of the layer
        lower_mantle.set_temperature_mode(
            temperature_mode='adiabatic',
            temperature_top=1900.)

        # And we set a self-consistent pressure. The pressure at the top of the layer and
        # gravity at the bottom of the layer are given by the PREM.
        pressure, gravity = burnman.seismic.PREM().evaluate(
            ['pressure', 'gravity'], depths)
        lower_mantle.set_pressure_mode(pressure_mode='self-consistent',
                                       pressure_top=pressure[-1], gravity_bottom=gravity[0])
        lower_mantle.make()

        # Constructing the tvel file for obspy. Here we use PREM to fill in the
        # rest of the planet
        burnman.output_seismo.write_tvel_file(
            lower_mantle,
            filename=modelname + '.tvel',
            background_model=burnman.seismic.PREM())
        
        # Plot ray paths and travel times
        plot_rays_and_times(modelname)

    # Second example implementing an entire planet
    if example_planet:
        modelname = 'planetzog'

        # A layer is defined by 4 parameters: Name, min_depth, max_depth,and number of slices within the layer.
        # Separately the composition and the temperature_mode need to set.
        # Radii of different layers
        radius_planet = 6371.e3
        radius_ic = 1220.e3
        radius_oc = 3580.e3
        radius_lm = 5711.e3

        # inner_core
        inner_core = burnman.Layer(
            'inner core', radii=np.linspace(
                0., radius_ic, 10))
        inner_core.set_material(burnman.minerals.other.Fe_Dewaele())

        # The minerals that make up our core do not currently implement the
        # thermal equation of state, so we will set the temperature at 300 K.
        inner_core.set_temperature_mode('user-defined',
                                        300. * np.ones_like(inner_core.radii))

        # outer_core
        outer_core = burnman.Layer(
            'outer core', radii=np.linspace(
                radius_ic, radius_oc, 10))
        outer_core.set_material(burnman.minerals.other.Liquid_Fe_Anderson())
        # The minerals that make up our core do not currently implement the
        # thermal equation of state, so we will define the temperature at 300
        # K.
        outer_core.set_temperature_mode(
            'user-defined',
            300. *
            np.ones_like(
                outer_core.radii))

        # Next the Mantle.
        lower_mantle = burnman.Layer(
            'lower mantle', radii=np.linspace(
                radius_oc, radius_lm, 10))
        lower_mantle.set_material(burnman.minerals.SLB_2011.mg_bridgmanite())
        lower_mantle.set_temperature_mode('adiabatic')
        upper_mantle = burnman.Layer(
            'upper mantle', radii=np.linspace(
                radius_lm, radius_planet, 10))
        upper_mantle.set_material(burnman.minerals.SLB_2011.forsterite())
        upper_mantle.set_temperature_mode('adiabatic', temperature_top=1200.)

        # Now we calculate the planet.
        planet_zog = burnman.Planet(
            'Planet Zog', [
                inner_core, outer_core, lower_mantle, upper_mantle], verbose=True)

        # Here we compute its state. Go BurnMan Go!
        # (If we were to change composition of one of the layers, we would have to
        # recompute the state)
        planet_zog.make()

        # Constructing the tvel file for obspy. Here we use PREM to fill in the
        # rest of the planet
        burnman.output_seismo.write_tvel_file(
            planet_zog,
            filename=modelname + '.tvel',
            background_model=burnman.seismic.PREM())

            
        # Plot ray paths and travel times
        plot_rays_and_times(modelname)
