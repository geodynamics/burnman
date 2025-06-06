# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2025 by the BurnMan team, released under the GNU
# GPL v2 or later.

import numpy as np

import burnman
from burnman import minerals
from burnman.utils.math import l2_norm_profile

if True:

    def test_set_state(runs=10000):
        g = burnman.minerals.SLB_2011.garnet()
        g.set_composition([0.1, 0.2, 0.4, 0.2, 0.1])
        x = 0.0
        for x in range(0, runs):
            g.set_state(10.0e9 + x, 1500.0)
            x += g.activities
        print(x)

    test_set_state(1)
    if False:
        import cProfile, pstats, StringIO

        pr = cProfile.Profile()
        pr.enable()
        test_set_state()
        pr.disable()
        pr.dump_stats("test_set_state.stats")
        s = StringIO.StringIO()
        sortby = "cumulative"
        ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
        ps.print_stats()
        print(s.getvalue())


if True:

    def test_grueneisen():
        for aa in range(0, 500):
            a = burnman.eos.SLB3()
            m = minerals.SLB_2011.fayalite()
            m.set_state(1.0e9, 1.0e3)
            x = 0
            for i in range(0, 10000):
                x += m.grueneisen_parameter
        return x

    test_grueneisen()
    # r=timeit.timeit('__main__.test_grueneisen()', setup="import __main__", number=3)
    # print("test_grueneisen", r)

if True:
    seismic_model = burnman.seismic.PREM()
    number_of_points = (
        10  # set on how many depth slices the computations should be done
    )
    depths = np.linspace(1000e3, 2500e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ["pressure", "density", "v_p", "v_s", "v_phi"], depths
    )

    temperature = burnman.geotherm.brown_shankland(depths)

    def calc_velocities(a, b, c):
        amount_perovskite = a
        pv = minerals.SLB_2011.mg_fe_perovskite([b, 1.0 - b, 0.0])
        fp = minerals.SLB_2011.ferropericlase([c, 1.0 - c])
        rock = burnman.Composite([pv, fp], [amount_perovskite, 1.0 - amount_perovskite])

        mat_rho, mat_vp, mat_vs = rock.evaluate(
            ["density", "v_phi", "v_s"], seis_p, temperature
        )
        return mat_vp, mat_vs, mat_rho

    def error(a, b, c):
        mat_vp, mat_vs, mat_rho = calc_velocities(a, b, c)

        vs_err = l2_norm_profile(depths, mat_vs, seis_vs) / 1.0e3
        vp_err = l2_norm_profile(depths, mat_vp, seis_vp) / 1.0e3
        den_err = l2_norm_profile(depths, mat_rho, seis_rho) / 1.0e3

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
