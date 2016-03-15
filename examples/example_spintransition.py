# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_spintransition
----------------------

This example shows the different minerals that are implemented with a spin
transition.  Minerals with spin transition are implemented by defining two
separate minerals (one for the low and one for the high spin state).  Then a
third dynamic mineral is created that switches between the two previously
defined minerals by comparing the current pressure to the transition pressure.

*Specifically uses:*

* :func:`burnman.mineral_helpers.HelperSpinTransition`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase_HS`
* :func:`burnman.minerals.Murakami_etal_2012.fe_periclase_LS`


*Demonstrates:*

* implementation of spin transition in (Mg,Fe)O at user defined pressure
"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals

if __name__ == "__main__":
    # seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()
    number_of_points = 40  # set on how many depth slices the computations should be done
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    # alternatively, we could use the values where prem is defined:
    # depths = seismic_model.internal_depth_list(mindepth=700.e3,
    # maxdepth=2800.e3)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    # here we use the Brown & Shankland geotherm
    temperature = burnman.geotherm.brown_shankland(seis_p)

    # We create one mineral that contains spin transitions
    rock = minerals.Murakami_etal_2012.fe_periclase()

    # The mineral Murakami_fe_periclase is derived from minerals.helper_spin_transition
    # which contains the logic to switch between two other minerals based on the
    # current pressure. The mineral is implemented similar to the following lines:
    #
    #   class Murakami_fe_periclase(helper_spin_transition):
    #     def __init__(self):
    #       helper_spin_transition.__init__(self, 63.0e9, Murakami_fe_periclase_LS(), Murakami_fe_periclase_HS())
    #
    # Note the reference to the low spin and high spin minerals (_LS and _HS).

    # Now we calculate the velocities
    mat_rho, mat_vs, mat_vphi  = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    print("Calculations are done for:")
    rock.debug_print()

    # plot example 1
    plt.subplot(2, 2, 1)
    plt.plot(
        seis_p / 1.e9, mat_vs / 1.e3, color='b', linestyle='-', marker='o',
        markerfacecolor='b', markersize=4, label='Vs')
    plt.plot(
        seis_p / 1.e9, mat_vphi / 1.e3, color='r', linestyle='-', marker='o',
        markerfacecolor='r', markersize=4, label='Vp')
    plt.plot(
        seis_p / 1.e9, mat_rho / 1.e3, color='k', linestyle='-', marker='o',
        markerfacecolor='k', markersize=4, label='rho')
    plt.title("ferropericlase (Murakami et al. 2012)")
    plt.xlim(min(seis_p) / 1.e9, max(seis_p) / 1.e9)
    plt.ylim(5, 12)
    plt.legend(loc='upper left')

    # example 2: Here we show the effects of using purely High Spin or Low Spin

    rock = minerals.Murakami_etal_2012.fe_periclase_LS()

    mat_rho_LS, mat_vs_LS, mat_vphi_LS = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    rock = minerals.Murakami_etal_2012.fe_periclase_HS()
    mat_rho_HS, mat_vs_HS, mat_vphi_HS = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    rock = minerals.Murakami_etal_2012.fe_periclase()
    mat_rho_ON, mat_vs_ON, mat_vphi_ON = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    plt.subplot(2, 2, 2)
    plt.plot(
        seis_p / 1.e9, mat_vs_LS / 1.e3, color='b', linestyle='-', marker='.',
        markerfacecolor='b', markersize=4, label='Vs LS')
    plt.plot(
        seis_p / 1.e9, mat_vs_HS / 1.e3, color='r', linestyle='-', marker='.',
        markerfacecolor='b', markersize=4, label='Vs HS')
    plt.plot(
        seis_p / 1.e9, mat_vs_ON / 1.e3, color='g', linestyle='-', marker='o',
        markerfacecolor='b', markersize=4, label='Vs ON')
    plt.title("Murakami_fp")
    plt.xlim(min(seis_p) / 1.e9, max(seis_p) / 1.e9)
    # plt.ylim(300,800)
    plt.legend(loc='lower right')

    # Example 3: Periclase from Speziale et al. 2006
    # Here the compositions are implemented as fixed minerals.
    # For other options see example_composition.py
    rock = minerals.other.Speziale_fe_periclase()

    mat_rho, mat_vs, mat_vphi = \
        rock.evaluate(['density', 'v_s', 'v_phi'], seis_p, temperature)

    print("Calculations are done for:")
    rock.debug_print()

    # plot example 3
    plt.subplot(2, 2, 3)
    plt.plot(
        seis_p / 1.e9, mat_rho / 1.e3, color='k', linestyle='-', marker='o',
        markerfacecolor='k', markersize=4, label='rho')
    plt.title("ferropericlase (Speziale et al. 2007)")
    plt.xlim(min(seis_p) / 1.e9, max(seis_p) / 1.e9)
    plt.legend(loc='upper left')

    plt.subplot(2, 2, 4)
    plt.plot(
        seis_p / 1.e9, mat_vphi / 1.e3, color='b', linestyle='-', marker='o',
        markerfacecolor='b', markersize=4, label='Vphi')
    plt.title("ferropericlase (Speziale et al. 2007)")
    plt.xlim(min(seis_p) / 1.e9, max(seis_p) / 1.e9)
    plt.legend(loc='upper left')

    plt.savefig("output_figures/examples_spintransition.png")
    plt.show()
