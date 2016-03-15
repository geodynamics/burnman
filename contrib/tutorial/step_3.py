# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU GPL v2 or later.
# Ian Rose ian.r.rose@gmail.com

"""

CIDER 2014 BurnMan Tutorial --- step 3
--------------------------------------

In the previous two steps of the tutorial we tried to find a very simple
mineralogical model that best fit the 1D seismic model PREM.  But we
know that there is consideral uncertainty in many of the mineral
physical parameters that control how the elastic properties of minerals
change with pressure and temperature.  In this step we explore how uncertainties
in these parameters might affect the conclusions you draw.

The strategy here is to make many different "realizations" of the rock
that you determined was the closest fit to PREM, where each realization
has its mineral physical parameters perturbed by a small amount, hopefully
related to the uncertainty in that parameter.  In particular, we will
look at how perturbations to :math:`K_{0}^{'}` and :math:`G_{0}^{'}` (the pressure derivatives of the
bulk and shear modulus, respectively) change the calculated 1D seismic
profiles.

This script may be run by typing

    python step_3.py

"""
from __future__ import absolute_import
from __future__ import print_function

# We import a couple extra modules for this step than in the previous
# ones.  In particular, we will use a function for a normal distribution
# from numpy.random, as well as a class for doing linear 1D interpolation
# from scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt

import numpy.random
from scipy import interpolate

# We also import these modules for some more advanced plotting.
import numpy.ma as ma
from matplotlib.colors import LinearSegmentedColormap

# hack to allow scripts to be placed in subdirectories next to burnman:
import os
import sys
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

# The BurnMan imports, however, are the same.
import burnman
from burnman import minerals


if __name__ == '__main__':
    """
    Here we define a function realize_mineral() which takes a mineral
    as a paramter, perturbs its :math:`K_{0}^{'}` and :math:`G_{0}^{'}` parameters, then returns
    the modified mineral.  Typical values for :math:`K_{0}^{'}` are around 4 and typical
    values for :math:`G_{0}^{'}` are around 2 (both parameters are dimensionless).
    We perturb the mineral by drawing from a normal ditribution with a
    standard deviation given by K_prime_std_dev and G_prime_std_dev.
    These parameters are currently set to very small values but you will
    want to put in realistic values (perhaps 0.1-0.2)
    """

    def realize_mineral(mineral):
        K_prime_std_dev = 0.0001  # <----------------- One sigma uncertainty in K prime
        G_prime_std_dev = 0.0001  # <----------------- One sigma uncertainty in G prime

        # Perturb the G' and K' values by a random number drawn from a normal
        # distribution
        mineral.params['Kprime_0'] = mineral.params[
            'Kprime_0'] + numpy.random.normal(scale=K_prime_std_dev)
        mineral.params['Gprime_0'] = mineral.params[
            'Gprime_0'] + numpy.random.normal(scale=G_prime_std_dev)

    """
    Here we define a function realize_rock(), which creates the realization of
    our mantle model.  First it makes a mineral, then calls realize_mineral()
    on that mineral to perturb its values, and then makes the rock out of
    those perturbed minerals.  Put in your preferred mantle model from steps
    1 and 2 here.
    """

    def realize_rock():

        phase_1_fraction = 0.5
        phase_2_fraction = 1.0 - phase_1_fraction

        # Setup the minerals for the two phase composite.  This is different
        # from how we did it in step 1 and 2.  Instead, we create the two phases
        # then call realize_mineral() on them to perturb their properties.
        phase_1 = minerals.SLB_2011.stishovite()
        realize_mineral(phase_1)
        phase_2 = minerals.SLB_2011.wuestite()
        realize_mineral(phase_2)

        # Set up the rock with the now-perturbed mineral phases
        mantle_rock = burnman.Composite(
            [phase_1, phase_2], [phase_1_fraction, phase_2_fraction])
        mantle_rock.set_method('slb3')

        # Give back the realization of the rock with the perturbed phases.
        return mantle_rock

    # Again, we set up the seismic model and ask for its properties
    # in the lower mantle.  Basically the same thing as in steps 1 and 2.
    seismic_model = burnman.seismic.PREM()
    min_depth = 850.e3
    max_depth = 2800.e3
    n_depths = 10
    depths = np.linspace(min_depth, max_depth, n_depths)
    pressure, seis_rho, seis_vphi, seis_vs = seismic_model.evaluate(
        ['pressure', 'density', 'v_phi', 'v_s'], depths)
    pressures_sampled = np.linspace(
        pressure[0], pressure[-1], 20 * len(pressure))

    temperature = burnman.geotherm.brown_shankland(pressure)

    """
    Finally, we make 1000 realizations of the rock and save the seismic profiles
    for later visualization.
    """

    n_realizations = 1000
    outfile = open('uncertainty.dat', 'wb')

    for i in range(n_realizations):

        print("realization", i + 1)
        try:
            # We call the realize_rock() to create the ith model
            rock = realize_rock()

            # Calculate the wavespeed profiles, just as before.
            rho, vphi, vs = rock.evaluate(
                ['rho', 'v_phi', 'v_s'], pressure, temperature)

            # This block of code interpolates the resulting densities and wavespeeds
            # to a higher resolution line to make the plot look nicer.
            func_rho = interpolate.interp1d(pressure, rho)
            func_vs = interpolate.interp1d(pressure, vs)
            func_vphi = interpolate.interp1d(pressure, vphi)

            pressure_list = pressures_sampled
            density_list = func_rho(pressures_sampled)
            vs_list = func_vs(pressures_sampled)
            vphi_list = func_vphi(pressures_sampled)

            # Save the output to a file
            data = list(zip(pressure_list, vs_list, vphi_list, density_list))
            np.savetxt(outfile, data, fmt='%.10e', delimiter='\t')

        # It is possible for the Birch-Murnaghan equation of state to go unstable for
        # some values of the parameters, which can make it fail.  If that happens, we
        # simply disregard this realization of the rock.
        except ValueError:
            print("failed, skipping")

    """
    The rest of this script is concerned with plotting the results so that
    you may visualize the uncertainty space that we have sampled with our
    mineral model.  In different colors we see 2D histograms, where the color
    intensity corresponds to how many of the models went through that portion of
    the velocity/density space.  The dashed lines show the PREM values.  Tighter
    histograms indicate that the perturbations to the mineral physical values
    make less of a difference, more spread out ones show larger effects.

    In general, this plotting is more complex than the ones in step 1 and step 2
    and you may safely ignore it.
    """

    # Read in the data from uncertainty.dat
    outfile.close()
    data = np.loadtxt('uncertainty.dat')
    pressure_list = data[:, 0]
    vs_list = data[:, 1]
    vphi_list = data[:, 2]
    density_list = data[:, 3]

    # Create 2D histograms showing the spread of the data
    density_hist, rho_xedge, rho_yedge = np.histogram2d(
        pressure_list, density_list, bins=len(pressures_sampled), normed=True)
    vs_hist, vs_xedge, vs_yedge = np.histogram2d(
        pressure_list, vs_list, bins=len(pressures_sampled), normed=True)
    vphi_hist, vphi_xedge, vphi_yedge = np.histogram2d(
        pressure_list, vphi_list, bins=len(pressures_sampled), normed=True)

    vs_xedge /= 1.e9
    vphi_xedge /= 1.e9
    rho_xedge /= 1.e9
    vs_yedge /= 1.e3
    vphi_yedge /= 1.e3
    rho_yedge /= 1.e3

    left_edge = min(vs_xedge[0], vphi_xedge[0], rho_xedge[0])
    right_edge = max(vs_xedge[-1], vphi_xedge[-1], rho_xedge[-1])
    bottom_edge = 4.3
    top_edge = 11.3
    aspect_ratio = float((right_edge - left_edge) / (top_edge - bottom_edge))
    gamma = 0.5  # Mess with this to change intensity of colormaps near the edges

    plt.subplot(111, aspect='equal')
    plt.xlim(left_edge, right_edge)
    plt.ylim(bottom_edge, top_edge)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Wave Speed (km/s)')

    # plot density
    density_hist = ma.masked_where(density_hist <= 0.0, density_hist)
    c = LinearSegmentedColormap.from_list('vphi', [(0, '#ffffff'), (0.2, '#edf8e9'), (
        0.4, '#bae4b3'), (0.6, '#74c476'), (0.8, '#31a354'), (1.0, '#006d2c')], gamma=0.1)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        density_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[rho_xedge[0], rho_xedge[-1], rho_yedge[0], rho_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_rho / 1.e3, linestyle="--",
             color='g', linewidth=2.0, label='Density')

    # plot v_s
    vs_hist = ma.masked_where(vs_hist <= 0.0, vs_hist)
    c = LinearSegmentedColormap.from_list('vphi', [(0, '#ffffff'), (0.2, '#eff3ff'), (
        0.4, '#bdd7e7'), (0.6, '#6baed6'), (0.8, '#3182bd'), (1.0, '#08519c')], gamma=0.5)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        vs_hist.transpose(), origin='low', cmap=c,  interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[vs_xedge[0], vs_xedge[-1], vs_yedge[0], vs_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_vs / 1.e3, linestyle="--",
             color='b', linewidth=2.0, label='Vs')

    # plot v_phi
    vphi_hist = ma.masked_where(vphi_hist <= 0.0, vphi_hist)
    c = LinearSegmentedColormap.from_list('vphi', [(0, '#ffffff'), (0.2, '#fee5d9'), (
        0.4, '#fcae91'), (0.6, '#fb6a4a'), (0.8, '#de2d26'), (1.0, '#a50f15')], gamma=0.5)
    c.set_bad('w', alpha=1.0)
    plt.imshow(
        vphi_hist.transpose(), origin='low', cmap=c, interpolation='gaussian', alpha=.7,
        aspect=aspect_ratio, extent=[vphi_xedge[0], vphi_xedge[-1], vphi_yedge[0], vphi_yedge[-1]])
    plt.plot(pressure / 1.e9, seis_vphi / 1.e3,
             linestyle="--", color='r', linewidth=2.0, label='Vphi')

    plt.legend(loc='lower right')

    plt.show()
