# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

paper_opt_pv
------------

This script reproduces :cite:`Cottaar2014`, Figure 6.
Vary the amount perovskite vs. ferropericlase and compute the error in the
seismic data against PREM.

requires:
- creating minerals
- compute seismic velocities
- geotherms
- seismic models
- seismic comparison

teaches:
- compare errors between models
- loops over models

"""
from __future__ import absolute_import
from __future__ import print_function

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../../burnman'):
    sys.path.insert(1, os.path.abspath('../..'))

import burnman
from burnman import minerals
from misc.helper_solid_solution import *
import misc.colors as colors

if __name__ == "__main__":

    # figsize=(6,5)
    plt.figure(dpi=100, figsize=(12, 10))
    prop = {'size': 12}
    plt.rc('text', usetex=True)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{relsize}'
    plt.rc('font', family='sans-serif')
    figsize = (6, 5)

    dashstyle2 = (7, 3)
    dashstyle3 = (3, 2)

    # input variables ###
    #

    # INPUT for method
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
    method = 'slb3'

    seismic_model = burnman.seismic.PREM()
                                         # pick from .prem() .slow() .fast()
                                         # (see burnman/seismic.py)
    number_of_points = 20  # set on how many depth slices the computations should be done
    depths = np.linspace(850e3, 2700e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    print(seis_p[0], seis_p[-1])

    # temperature = burnman.geotherm.brown_shankland(seis_p)

    def eval_material(amount_perovskite):
        rock = burnman.Composite([SLB_2011_ZSB_2013_mg_fe_perovskite(0.07),
                                  other_ferropericlase(0.2)],
                                 [amount_perovskite, 1.0 - amount_perovskite])
        rock.set_method(method)
        temperature = burnman.geotherm.adiabatic(seis_p, 1900, rock)
        print("Calculations are done for:")
        rock.debug_print()

        mat_rho, mat_vs, mat_vphi = rock.evaluate(
            ['rho', 'v_s', 'v_phi'], seis_p, temperature)
        #[rho_err,vphi_err,vs_err]=burnman.compare_chifactor(mat_vs,mat_vphi,mat_rho,seis_vs,seis_vphi,seis_rho)

        return seis_p, mat_vs, mat_vphi, mat_rho

    def material_error(x):
        _, mat_vs, mat_vphi, mat_rho = eval_material(x)
        [vs_err, vphi_err, rho_err] = burnman.compare_l2(depths,
                                                         [mat_vs, mat_vphi,
                                                          mat_rho],
                                                         [seis_vs, seis_vphi, seis_rho])
        scale = 2700e3 - 850e3
        return vs_err / scale, vphi_err / scale

    xx = np.linspace(0.0, 1.0, 200)  # 200 for final image
    # speed up computation for the automatic tests:
    if "RUNNING_TESTS" in globals():
        xx = np.linspace(0.0, 1.0, 10)

    errs = np.array([material_error(x) for x in xx])
    yy_vs = errs[:, 0]
    yy_vphi = errs[:, 1]
    vs_average_prem = sum(seis_vs) / len(seis_vs)
    vphi_average_prem = sum(seis_vphi) / len(seis_vphi)
    print(vs_average_prem, vphi_average_prem)
    yy_vs /= vs_average_prem
    yy_vphi /= vphi_average_prem
    yy_sum = (yy_vs + yy_vphi)  # we scale by a factor so it fits in the plot
 #   plt.figure(dpi=100,figsize=figsize)
    plt.subplot(2, 2, 1)
    plt.plot(xx * 100, yy_vs, "-", color=colors.color(1),
             label=("$V_s$ error"), linewidth=1.5, dashes=dashstyle2)
    plt.plot(xx * 100, yy_vphi, "-", color=colors.color(3),
             label=("$V_\phi$ error"), linewidth=1.5)
    # plt.plot (xx*100,yy_vs+yy_vphi,"g--",label=("sum"),linewidth=1.5)
    plt.plot(xx * 100, yy_sum, "-", color=colors.color(4),
             label=("weighted sum"), linewidth=1.5, dashes=dashstyle3)

    ymin = 1e-2
    ymax = 1e2
    plt.ylim([ymin, ymax])

    print(xx[np.argmin(yy_vs)], xx[np.argmin(yy_vphi)], xx[np.argmin(yy_sum)])

    B = np.around(xx[np.argmin(yy_vs)], decimals=3)
    A = np.around(xx[np.argmin(yy_vphi)], decimals=3)
    C = np.around(xx[np.argmin(yy_sum)], decimals=3)

    plt.plot([A * 100., A * 100.], [ymin, ymax], color=colors.color(3),
             label='A (%g\%% pv)' % (A * 100), linewidth=1.5, linestyle='-')
    plt.plot([B * 100., B * 100.], [ymin, ymax], color=colors.color(1),
             label='B (%g\%% pv)' % (B * 100), linewidth=1.5, dashes=dashstyle2)
    plt.plot([C * 100., C * 100.], [ymin, ymax], color=colors.color(4),
             label='C (%g\%% pv)' % (C * 100), linewidth=1.5, dashes=dashstyle3)

    plt.yscale('log')
    plt.xlabel('\% Perovskite')
    plt.ylabel('Error')
    plt.legend(loc='lower left', prop=prop)
#    plt.tight_layout(pad=2)
#    plt.savefig("opt_pv_1.pdf",bbox_inches='tight')
#    plt.show()

    A_p, A_vs, A_vphi, _ = eval_material(A)
    B_p, B_vs, B_vphi, _ = eval_material(B)
    C_p, C_vs, C_vphi, _ = eval_material(C)

#    plt.figure(dpi=100,figsize=figsize)
    plt.subplot(2, 2, 3)
    plt.plot(seis_p / 1.e9, seis_vs / 1.e3, color='k', linestyle='-',
             linewidth=2.0, markersize=6, markerfacecolor='None', label='PREM')
    plt.plot(A_p / 1.e9, A_vs / 1.e3, color=colors.color(3), linestyle='-',
             label='A (%g\%% pv)' % (A * 100), linewidth=1.5, markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(3))
    plt.plot(B_p / 1.e9, B_vs / 1.e3, color=colors.color(1), dashes=dashstyle2,
             label='B (%g\%% pv)' % (B * 100), linewidth=1.5, markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(1))
    plt.plot(C_p / 1.e9, C_vs / 1.e3, color=colors.color(4), dashes=dashstyle3,
             label='C (%g\%% pv)' % (C * 100), linewidth=1.5, markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(4))
    plt.xlabel('Pressure (GPa)')
    plt.ylabel(
        'Shear velocity $V_{\mathlarger{\mathlarger{\mathlarger{s}}}}$ (km/s)')
    plt.xlim([30, 130])
    plt.legend(loc='lower right', prop=prop)
  #  plt.tight_layout()
 #   plt.savefig("opt_pv_2.pdf",bbox_inches='tight')
#    plt.show()

    plt.subplot(2, 2, 4)

#    plt.figure(dpi=100,figsize=figsize)
    plt.plot(seis_p / 1.e9, seis_vphi / 1.e3, color='k', linestyle='-',
             linewidth=2.0, markersize=6, markerfacecolor='None', label='PREM', mew=1.5)
    plt.plot(A_p / 1.e9, A_vphi / 1.e3, color=colors.color(3), linestyle='-',
             markevery=5, marker='s', markersize=5, markeredgecolor=colors.color(3), markerfacecolor='None', mew=1.5, label='A (%g\%% pv)' % (A * 100), linewidth=1.5)
    plt.plot(
        B_p / 1.e9, B_vphi / 1.e3, color=colors.color(1), dashes=dashstyle2,
        markevery=5, marker='s', markersize=5, markeredgecolor=colors.color(1), markerfacecolor='None', mew=1.5, label='B (%g\%% pv)' % (B * 100), linewidth=1.5)
    plt.plot(
        C_p / 1.e9, C_vphi / 1.e3, color=colors.color(4), dashes=dashstyle3,
        markevery=5, marker='s', markersize=5, markeredgecolor=colors.color(4), markerfacecolor='None', mew=1.5, label='C (%g\%% pv)' % (C * 100), linewidth=1.5)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel(
        "Bulk sound velocity $V_{\mathlarger{\mathlarger{\mathlarger{\phi}}}}$ (km/s)")
    plt.xlim([30, 130])
    plt.legend(loc='lower right', prop=prop)
    # plt.tight_layout()
    # plt.savefig("opt_pv_3.pdf",bbox_inches='tight')
    # plt.show()

    # plot percent differences
#    plt.figure(dpi=100,figsize=figsize)
    plt.subplot(2, 2, 2)
    plt.plot(seis_p / 1.e9, seis_vs * 0.0,
             color='k', linestyle='-', linewidth=2.0)
    plt.plot(seis_p / 1.e9, (A_vs - seis_vs) / seis_vs * 100.0, color=colors.color(3), label='$V_s$: A (%g\%% pv)' %
             (A * 100), linewidth=1.5, linestyle='-', markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(3))
    plt.plot(seis_p / 1.e9, (B_vs - seis_vs) / seis_vs * 100.0, color=colors.color(1), label='$V_s$: B (%g\%% pv)' %
             (B * 100), linewidth=1.5, dashes=dashstyle2, markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(1))
    plt.plot(seis_p / 1.e9, (C_vs - seis_vs) / seis_vs * 100.0, color=colors.color(4), label='$V_s$: C (%g\%% pv)' %
             (C * 100), linewidth=1.5, dashes=dashstyle3, markevery=5, marker='v', markerfacecolor='None', mew=1.5, markeredgecolor=colors.color(4))
    plt.plot(
        seis_p / 1.e9, (A_vphi - seis_vphi) / seis_vphi * 100.0, color=colors.color(3), markevery=5, marker='s', markersize=5,
        markeredgecolor=colors.color(3), markerfacecolor='None', mew=1.5, label='$V_\phi$: A', linewidth=1.5, linestyle='-')
    plt.plot(
        seis_p / 1.e9, (B_vphi - seis_vphi) / seis_vphi * 100.0, color=colors.color(1), markevery=5, marker='s', markersize=5,
        markeredgecolor=colors.color(1), markerfacecolor='None', mew=1.5, label='$V_\phi$: B', linewidth=1.5, dashes=dashstyle2)
    plt.plot(
        seis_p / 1.e9, (C_vphi - seis_vphi) / seis_vphi * 100.0, color=colors.color(4), markevery=5, marker='s', markersize=5,
        markeredgecolor=colors.color(4), markerfacecolor='None', mew=1.5, label='$V_\phi$: C', linewidth=1.5, dashes=dashstyle3)
    plt.xlabel('Pressure (GPa)')
    plt.ylabel('Difference from PREM (\%)')
    plt.ylim([-5, 4])
    plt.xlim([30, 130])
    plt.legend(loc='lower center', ncol=2, prop=prop)
#    plt.tight_layout()
#    plt.savefig("opt_pv_4.pdf",bbox_inches='tight')
#    plt.show()

    plt.tight_layout()
    if "RUNNING_TESTS" not in globals():
        plt.savefig("paper_opt_pv.pdf", bbox_inches='tight')
    plt.show()
