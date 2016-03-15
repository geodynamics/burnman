# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_averaging
-----------------

This example shows the effect of different averaging schemes. Currently four
averaging schemes are available:

1. Voight-Reuss-Hill
2. Voight averaging
3. Reuss averaging
4. Hashin-Shtrikman averaging

See :cite:`Watt1976` Journal of Geophysics and Space Physics for explanations
of each averaging scheme.

*Specifically uses:*

* :class:`burnman.averaging_schemes.VoigtReussHill`
* :class:`burnman.averaging_schemes.Voigt`
* :class:`burnman.averaging_schemes.Reuss`
* :class:`burnman.averaging_schemes.HashinShtrikmanUpper`
* :class:`burnman.averaging_schemes.HashinShtrikmanLower`

*Demonstrates:*

* implemented averaging schemes

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
    """ choose 'slb2' (finite-strain 2nd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'slb3 (finite-strain 3rd order shear modulus,
        stixrude and lithgow-bertelloni, 2005)
    or 'mgd3' (mie-gruneisen-debeye 3rd order shear modulus,
        matas et al. 2007)
    or 'mgd2' (mie-gruneisen-debeye 2nd order shear modulus,
        matas et al. 2007)
    or 'bm2' (birch-murnaghan 2nd order, if you choose to ignore temperature
       (your choice in geotherm will not matter in this case))
       or 'bm3' (birch-murnaghan 3rd order, if you choose to ignore temperature
    (your choice in geotherm will not matter in this case))"""

    amount_perovskite = 0.6

    rock = burnman.Composite([minerals.SLB_2011.mg_perovskite(),
                              minerals.SLB_2011.periclase()],
                             [amount_perovskite, 1.0 - amount_perovskite])

    perovskitite = minerals.SLB_2011.mg_perovskite()

    periclasite = minerals.SLB_2011.periclase()

    # seismic model for comparison:
    # pick from .prem() .slow() .fast() (see burnman/seismic.py)
    seismic_model = burnman.seismic.PREM()
    # set on how many depth slices the computations should be done
    number_of_points = 20
    # we will do our computation and comparison at the following depth values:
    depths = np.linspace(700e3, 2800e3, number_of_points)
    # alternatively, we could use the values where prem is defined:
    # depths = seismic_model.internal_depth_list(mindepth=700.e3,
    # maxdepth=2800.e3)
    pressures, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    temperatures = burnman.geotherm.brown_shankland(pressures)

    print("Calculations are done for:")
    rock.debug_print()

    # calculate the seismic velocities of the rock using a whole battery of
    # averaging schemes:

    # do the end members, here averaging scheme does not matter (though it
    # defaults to Voigt-Reuss-Hill)
    model_pv = burnman.Model(
        perovskitite, pressures, temperatures, burnman.averaging_schemes.VoigtReussHill())
    model_fp = burnman.Model(
        periclasite, pressures, temperatures, burnman.averaging_schemes.VoigtReussHill())

    # Voigt Reuss Hill / Voigt / Reuss averaging
    model_vrh = burnman.Model(
        rock, pressures, temperatures, burnman.averaging_schemes.VoigtReussHill())
    model_v = burnman.Model(
        rock, pressures, temperatures, burnman.averaging_schemes.Voigt())
    model_r = burnman.Model(
        rock, pressures, temperatures, burnman.averaging_schemes.Reuss())

    # Upper/lower bound for Hashin-Shtrikman averaging
    model_hsu = burnman.Model(
        rock, pressures, temperatures, burnman.averaging_schemes.HashinShtrikmanUpper())
    model_hsl = burnman.Model(
        rock, pressures, temperatures, burnman.averaging_schemes.HashinShtrikmanLower())

    # PLOTTING
    # plot vs
    fig = plt.figure()

    plt.plot(
        pressures / 1.e9, model_v.v_s() / 1.e3, color='c', linestyle='-', marker='^',
        markersize=4, label='Voigt')
    plt.plot(
        pressures / 1.e9, model_r.v_s() / 1.e3, color='k', linestyle='-', marker='v',
        markersize=4, label='Reuss')
    plt.plot(
        pressures / 1.e9, model_vrh.v_s() / 1.e3, color='b', linestyle='-', marker='x',
        markersize=4, label='Voigt-Reuss-Hill')
    plt.plot(
        pressures / 1.e9, model_hsu.v_s() / 1.e3, color='r', linestyle='-', marker='x',
        markersize=4, label='Hashin-Shtrikman')
    plt.plot(
        pressures / 1.e9, model_hsl.v_s() / 1.e3, color='r', linestyle='-', marker='x',
        markersize=4)
    plt.plot(
        pressures / 1.e9, model_pv.v_s() / 1.e3, color='y', linestyle='-', marker='x',
        markersize=4, label='Mg Perovskite')
    plt.plot(
        pressures / 1.e9, model_fp.v_s() / 1.e3, color='g', linestyle='-', marker='x',
        markersize=4, label='Periclase')
    plt.xlim(min(pressures) / 1.e9, max(pressures) / 1.e9)
    plt.legend(loc='upper left', prop={'size': 11}, frameon=False)
    plt.xlabel('pressure (GPa)')
    plt.ylabel('Vs (km/s)')

    vs_pv_norm = (model_pv.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_fp_norm = (model_fp.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_vrh_norm = (model_vrh.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_v_norm = (model_v.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_r_norm = (model_r.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_hsu_norm = (model_hsu.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())
    vs_hsl_norm = (model_hsl.v_s() - model_fp.v_s()) / \
        (model_pv.v_s() - model_fp.v_s())

    ax = fig.add_axes([0.58, 0.18, 0.3, 0.3])
    plt.plot(pressures / 1.e9, vs_v_norm, color='c', linestyle='-', marker='^',
             markersize=4, label='Voigt')
    plt.plot(pressures / 1.e9, vs_r_norm, color='k', linestyle='-', marker='v',
             markersize=4, label='Reuss')
    plt.plot(
        pressures / 1.e9, vs_vrh_norm, color='b', linestyle='-', marker='x',
        markersize=4, label='Voigt-Reuss-Hill')
    plt.plot(
        pressures / 1.e9, vs_hsl_norm, color='r', linestyle='-', marker='x',
        markersize=4, label='Hashin-Shtrikman')
    plt.plot(
        pressures / 1.e9, vs_hsu_norm, color='r', linestyle='-', marker='x',
        markersize=4)
    plt.plot(
        pressures / 1.e9, vs_pv_norm, color='y', linestyle='-', marker='x',
        markersize=4, label='Mg Perovskite')
    plt.plot(
        pressures / 1.e9, vs_fp_norm, color='g', linestyle='-', marker='x',
        markersize=4, label='Periclase')
    ax.tick_params(labelsize=10)
    plt.title("normalized by mixture endmembers", fontsize=10)
    plt.xlim(min(pressures) / 1.e9, max(pressures) / 1.e9)
    plt.ylim(-0.005, 1.005)
    plt.xlabel('pressure (GPa)', fontsize=10)
    plt.ylabel('normalized Vs', fontsize=10)
    # plt.legend(loc='lower right')

    plt.savefig("output_figures/example_averaging_normalized.png")
    plt.show()
