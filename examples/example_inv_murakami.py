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

import pymc
import cProfile
import scipy.stats as sp
import matplotlib.mlab as mlab


if __name__ == "__main__":
    seismic_model = burnman.seismic.PREM()
                                         # pick from .prem() .slow() .fast()
                                         # (see code/seismic.py)
    number_of_points = 10  # set on how many depth slices the computations should be done
    depths = np.linspace(1000e3, 2500e3, number_of_points)
    seis_p, seis_rho, seis_vp, seis_vs, seis_vphi = seismic_model.evaluate(
        ['pressure', 'density', 'v_p', 'v_s', 'v_phi'], depths)

    temperature = burnman.geotherm.brown_shankland(seis_p)

    print("preparations done")

    def calc_velocities(amount_pv, iron_pv, iron_fp):
        pv = minerals.SLB_2011.mg_fe_perovskite([1.0 - iron_pv, iron_pv, 0.0])
        fp = minerals.SLB_2011.ferropericlase([1.0 - iron_fp, iron_fp])
        rock = burnman.Composite([pv, fp], [amount_pv, 1.0 - amount_pv])

        mat_rho, mat_vp, mat_vs = rock.evaluate(
            ['density', 'v_p', 'v_s'], seis_p, temperature)
        return mat_vp, mat_vs, mat_rho

    def error(amount_pv, iron_pv, iron_fp):
        mat_vp, mat_vs, mat_rho = calc_velocities(amount_pv, iron_pv, iron_fp)

        vs_err = burnman.l2(depths, mat_vs, seis_vs) / 1e9
        vp_err = burnman.l2(depths, mat_vp, seis_vp) / 1e9
        den_err = burnman.l2(depths, mat_rho, seis_rho) / 1e9

        # print vs_err, vp_err, den_err

        return vs_err + vp_err + den_err

    # Priors on unknown parameters:
    amount_pv = pymc.Uniform('amount_pv', lower=0.0, upper=1.0, value=0.5)
    iron_pv = pymc.Uniform('iron_pv', lower=0.0, upper=1.0, value=0.5)
    iron_fp = pymc.Uniform('iron_fp', lower=0.0, upper=1.0, value=0.5)

    minerr = 1e100

    @pymc.deterministic
    def theta(amount_pv=amount_pv, iron_pv=iron_pv, iron_fp=iron_fp):
        global minerr

        try:
            e = error(amount_pv, iron_pv, iron_fp)
            if (e < minerr):
                minerr = e
                print("best fit", e, "values:", amount_pv, iron_pv, iron_fp)
            return e
        except ValueError:
            return 1e20  # float("inf")

    sig = 10.0
    misfit = pymc.Normal('d', mu=theta, tau=1.0 / (
        sig * sig), value=0, observed=True, trace=True)
    model = [amount_pv, iron_pv, iron_fp, misfit]
    things = ['amount_pv', 'iron_pv', 'iron_fp', 'misfit']

    whattodo = ""

    if len(sys.argv) < 3:
        print("options:")
        print("run <dbname>")
        print("continue <dbname>")
        print("plot <dbname1> <dbname2> ...")
        print("show amount_pv iron_pv iron_fp")
    else:
        whattodo = sys.argv[1]
        dbname = sys.argv[2]

    if whattodo == "run":
        S = pymc.MCMC(model, db='pickle', dbname=dbname)
        S.sample(iter=400, burn=200, thin=1)
        S.db.close()

    if whattodo == "continue":
        n_runs = 50
        for l in range(0, n_runs):
            db = pymc.database.pickle.load(dbname)
            print("*** run=%d/%d, # samples: %d" %
                  (l, n_runs, db.trace('amount_pv').stats()['n']))
            S = pymc.MCMC(model, db=db)
            S.sample(iter=500, burn=0, thin=1)
            S.db.close()

    if whattodo == "plot":
        files = sys.argv[2:]
        print("files:", files)

        toburn = 0
        plot_idx = 1

        for t in things:
            if t == 'misfit':
                continue
            trace = []
            print("trace:", t)
            for filename in files:
                db = pymc.database.pickle.load(filename)
                dir(db)
                newtrace = db.trace(t, chain=None).gettrace(
                    burn=toburn, chain=None)
                if (trace != []):
                    trace = np.append(trace, newtrace)
                else:
                    trace = newtrace
                print("   adding ", newtrace.size, "burn = ", toburn)
            print("  total size ", trace.size)
            print("mean = ", trace.mean())
            for bin in [10, 20, 50, 100]:
                hist, bin_edges = np.histogram(trace, bins=bin)
                a = np.argmax(hist)
                print("maxlike = ", bin_edges[a], bin_edges[
                      a + 1], (bin_edges[a] + bin_edges[a + 1]) / 2.0)

            plt.subplot(2, len(things) / 2, plot_idx)
            if plot_idx == 2:
                n, bins, patches = plt.hist(
                    np.array(trace), 50,  normed=1, facecolor='green', alpha=0.75)

                X = sp.gumbel_l.fit(np.array(trace))
                print(X)
                dist = sp.gumbel_l(X[0], X[1])
                x = np.array(bins)
                y = dist.pdf(x)
                print(y)
                plt.plot(x, y, 'k--', linewidth=2)

                X = sp.norm.fit(np.array(trace))
                print(X)
                dist = sp.norm(X[0], X[1])
                x = np.array(bins)
                y = dist.pdf(x)
                plt.plot(x, y, 'r--', linewidth=2)

                X = sp.genextreme.fit(np.array(trace))
                print(X)
                dist = sp.genextreme(X[0], X[1], X[2])
                x = np.array(bins)
                y = dist.pdf(x)
                plt.plot(x, y, 'b--', linewidth=2)
                plt.title("%s" % (t), fontsize='small')

            elif plot_idx == 3:
                n, bins, patches = plt.hist(
                    np.array(trace), 50,  normed=1, facecolor='green', alpha=0.75)

                # copied = np.append(np.array(trace), -np.array(trace))
                #(mu, sigma) = sp.norm.fit(copied)
                # y = mlab.normpdf( bins, mu, sigma)
                # l = plt.plot(bins, y, 'r--', linewidth=2)

                X = sp.expon.fit(np.array(trace), floc=0)

                print(X)

                # X = sp.burr.fit(np.array(trace))
                dist = sp.expon(X[0], X[1])
                # print X
                print(bins)
                print(dist.pdf(np.array(bins)))
                plt.plot(bins, dist.pdf(np.array(bins)), 'r--', linewidth=2)
                plt.title("%s, mu: %.3e, sigma: %.3e" %
                          (t, mu, sigma), fontsize='small')

            else:
                (mu, sigma) = sp.norm.fit(np.array(trace))
                print("mu, sigma: %e %e" % (mu, sigma))
                n, bins, patches = plt.hist(
                    np.array(trace), 50,  normed=1, facecolor='green', alpha=0.75)
                y = mlab.normpdf(bins, mu, sigma)
                l = plt.plot(bins, y, 'r--', linewidth=2)
                plt.title("%s, mean: %.3e, std dev.: %.3e" %
                          (t, mu, sigma), fontsize='small')

            plot_idx += 1

        plt.savefig("output_figures/example_inv_murakami.png")
        plt.show()

    if whattodo == "misfittest":

        for a in np.linspace(0.4, 0.8, 10):
            for b in np.linspace(0.05, 0.2, 5):
                for c in np.linspace(0, 0.2, 5):
                    mat_vp, mat_vs, mat_rho = calc_velocities(a, b, c)
                    misfit = error(a, b, c)
                    print("misfit: %s " % misfit)
                    if misfit < 25:
                        plt.subplot(2, 2, 1)
                        plt.plot(
                            seis_p / 1.e9, mat_vs / 1.e3, color='r', linestyle='-',
                            marker='x', markerfacecolor='r', markersize=4)
                        plt.subplot(2, 2, 2)
                        plt.plot(
                            seis_p / 1.e9, mat_vp / 1.e3, color='r', linestyle='-',
                            marker='x', markerfacecolor='r', markersize=4)
                        plt.subplot(2, 2, 3)
                        plt.plot(
                            seis_p / 1.e9, mat_rho / 1.e3, color='r', linestyle='-',
                            marker='x', markerfacecolor='r', markersize=4, label='model 1')

        plt.subplot(2, 2, 1)
        plt.plot(seis_p / 1.e9, seis_vs / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='None', markersize=6)
        plt.ylim([4, 8])
        plt.title("Vs (km/s)")

        # plot Vp
        plt.subplot(2, 2, 2)
        plt.plot(seis_p / 1.e9, seis_vp / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='k', markersize=4)
        plt.ylim([10, 14])
        plt.title("Vp (km/s)")

        # plot density
        plt.subplot(2, 2, 3)
        plt.plot(seis_p / 1.e9, seis_rho / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='k', markersize=4, label='ref')
        plt.title("density (kg/m^3)")
        plt.legend(loc='upper left')
        plt.ylim([4, 8])
        plt.savefig("output_figures/example_inv_murakami_2.png")
        plt.show()

    if whattodo == "test":
        db = pymc.database.pickle.load(dbname)
        S = pymc.MCMC(model, db=db)

        for t in things:
            print(db.trace(t).stats())

        print("means:")
        for t in things:
            print(t, "\t", db.trace(t).stats()['mean'])

        print("#samples: %d" % db.trace('mg_pv_K').stats()['n'])

        pymc.raftery_lewis(S, q=0.025, r=0.01)

        b = 1
        t = 1

        scores = pymc.geweke(S, intervals=20)

        pymc.Matplot.trace(db.trace('deviance', chain=None).gettrace(
            burn=1000, thin=t, chain=None), 'deviance', rows=2, columns=9, num=1)

        pymc.Matplot.trace(db.trace('mg_pv_K', chain=None).gettrace(
            thin=t, chain=None), 'mg_pv_K', rows=2, columns=9, num=2)
        pymc.Matplot.histogram(
            np.array(db.trace('mg_pv_K', chain=None).gettrace(burn=b, chain=None)), 'mg_pv_K', rows=2, columns=9, num=11)

        pymc.Matplot.trace(db.trace('mg_pv_K_prime', chain=None).gettrace(
            thin=t, chain=None), 'mg_pv_K_prime', rows=2, columns=9, num=3)
        pymc.Matplot.histogram(np.array(db.trace('mg_pv_K_prime', chain=None).gettrace(
            burn=b, chain=None)), 'mg_pv_K_prime', rows=2, columns=9, num=12)

        pymc.Matplot.trace(db.trace('mg_pv_G', chain=None).gettrace(
            thin=t, chain=None), 'mg_pv_G', rows=2, columns=9, num=4)
        pymc.Matplot.histogram(
            np.array(db.trace('mg_pv_G', chain=None).gettrace(burn=b, chain=None)), 'mg_pv_G', rows=2, columns=9, num=13)

        pymc.Matplot.trace(db.trace('mg_pv_G_prime', chain=None).gettrace(
            thin=t, chain=None), 'mg_pv_G_prime', rows=2, columns=9, num=5)
        pymc.Matplot.histogram(np.array(db.trace('mg_pv_G_prime', chain=None).gettrace(
            burn=b, chain=None)), 'mg_pv_G_prime', rows=2, columns=9, num=14)

        pymc.Matplot.trace(db.trace('fe_pv_K', chain=None).gettrace(
            thin=t, chain=None), 'fe_pv_K', rows=2, columns=9, num=6)
        pymc.Matplot.histogram(
            np.array(db.trace('fe_pv_K', chain=None).gettrace(burn=b, chain=None)), 'fe_pv_K', rows=2, columns=9, num=15)

        pymc.Matplot.trace(db.trace('fe_pv_K_prime', chain=None).gettrace(
            thin=t, chain=None), 'fe_pv_K_prime', rows=2, columns=9, num=7)
        pymc.Matplot.histogram(np.array(db.trace('fe_pv_K_prime', chain=None).gettrace(
            burn=b, chain=None)), 'fe_pv_K_prime', rows=2, columns=9, num=16)

        pymc.Matplot.trace(db.trace('fe_pv_G', chain=None).gettrace(
            thin=t, chain=None), 'fe_pv_G', rows=2, columns=9, num=8)
        pymc.Matplot.histogram(
            np.array(db.trace('fe_pv_G', chain=None).gettrace(burn=b, chain=None)), 'fe_pv_G', rows=2, columns=9, num=17)

        pymc.Matplot.trace(db.trace('fe_pv_G_prime', chain=None).gettrace(
            thin=t, chain=None), 'fe_pv_G_prime', rows=2, columns=9, num=9)
        pymc.Matplot.histogram(np.array(db.trace('fe_pv_G_prime', chain=None).gettrace(
            burn=b, chain=None)), 'fe_pv_G_prime', rows=2, columns=9, num=18)

        plt.show()

    if whattodo == "show":
        values = [float(i) for i in sys.argv[2:]]
        mat_vp, mat_vs, mat_rho = calc_velocities(
            values[0], values[1], values[2])

        misfit = error(values[0], values[1], values[2])
        print("misfit: %s " % misfit)

        plt.suptitle('misfit %.3e, amount_pv=%.4f, iron_pv=%.4f, iron_fp=%.4f' %
                     (misfit, values[0], values[1], values[2]))

        plt.subplot(2, 2, 1)
        plt.plot(seis_p / 1.e9, mat_vs / 1.e3, color='r', linestyle='-',
                 marker='x', markerfacecolor='r', markersize=4)
        plt.plot(seis_p / 1.e9, seis_vs / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='None', markersize=6)
        plt.ylim([4, 8])
        plt.title("Vs (km/s)")

        # plot Vp
        plt.subplot(2, 2, 2)
        plt.plot(seis_p / 1.e9, mat_vp / 1.e3, color='r', linestyle='-',
                 marker='x', markerfacecolor='r', markersize=4)
        plt.plot(seis_p / 1.e9, seis_vp / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='k', markersize=4)
        plt.ylim([10, 14])
        plt.title("Vp (km/s)")

        # plot density
        plt.subplot(2, 2, 3)
        plt.plot(seis_p / 1.e9, mat_rho / 1.e3, color='r', linestyle='-',
                 marker='x', markerfacecolor='r', markersize=4, label='model 1')
        plt.plot(seis_p / 1.e9, seis_rho / 1.e3, color='k', linestyle='-',
                 marker='o', markerfacecolor='k', markersize=4, label='ref')
        plt.title("density (kg/m^3)")
        plt.legend(loc='upper left')
        plt.ylim([4, 8])
        plt.savefig("output_figures/example_inv_murakami_2.png")
        plt.show()

    if whattodo == "profile2":
        # run with:
        # python -m cProfile -o output.pstats example_inv_big_pv.py profile2 1
        # gprof2dot.py -f pstats output.pstats | dot -Tpng -o output.png
        [error(235.654790318e9, 3.87724833477, 165.45907725e9, 1.61618366689, 273.888690109e9, 4.38543140228, 306.635371217e9, 1.42610871557)
         for i in range(0, 10)]

    if whattodo == "profile":
        # just run normally
        cProfile.run(
            "[error(235.654790318e9, 3.87724833477, 165.45907725e9, 1.61618366689, 273.888690109e9, 4.38543140228, 306.635371217e9, 1.42610871557) for i in range(0,10)]")
