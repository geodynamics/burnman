from __future__ import absolute_import
from __future__ import print_function
# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
# hack to allow scripts to be placed in subdirectories next to burnman:
if not os.path.exists('burnman') and os.path.exists('../burnman'):
    sys.path.insert(1, os.path.abspath('..'))

import burnman
from burnman import minerals


import timeit

if True:
    def test_grueneisen():
        for aa in range(0, 500):
            a = burnman.eos.SLB3()
            m = minerals.SLB_2011.fayalite()
            x = 0
            for i in range(0, 10000):
                x += a.grueneisen_parameter(1e9, 1e4, 1.1, m.params)
        return x

    test_grueneisen()
    # r=timeit.timeit('__main__.test_grueneisen()', setup="import __main__", number=3)
    # print("test_grueneisen", r)

if True:
    seismic_model = burnman.seismic.PREM()
                                         # pick from .prem() .slow() .fast()
                                         # (see code/seismic.py)
    number_of_points = 10  # set on how many depth slices the computations should be done
    depths = np.linspace(1000e3, 2500e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    temperature = burnman.geotherm.brown_shankland(seis_p)

    def calc_velocities(a, b, c):
        amount_perovskite = a
        pv = minerals.SLB_2011.mg_fe_perovskite([b, 1.0 - b, 0.0])
        fp = minerals.SLB_2011.ferropericlase([c, 1.0 - c])
        rock = burnman.Composite(
            [pv, fp], [amount_perovskite, 1.0 - amount_perovskite])

        mat_rho, mat_vp, mat_vs = rock.evaluate(
            ['density', 'v_phi', 'v_s'], seis_p, temperature)
        return mat_vp, mat_vs, mat_rho

    def error(a, b, c):
        mat_vp, mat_vs, mat_rho = calc_velocities(a, b, c)

        vs_err = burnman.l2(depths, mat_vs, seis_vs) / 1e9
        vp_err = burnman.l2(depths, mat_vp, seis_vp) / 1e9
        den_err = burnman.l2(depths, mat_rho, seis_rho) / 1e9

        return vs_err + vp_err + den_err

    def run():
        n = 5
        m = 1e10
        for a in np.linspace(0.0, 1.0, n):
            for b in np.linspace(0.0, 1.0, n):
                for c in np.linspace(0.0, 1.0, n):
                    m = min(m, error(a, b, c))
        print(m)
        return m

    run()
    # r=timeit.timeit('__main__.run()', setup="import __main__", number=3)
    # print("test_seismic_velocities", r)
