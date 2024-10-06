from __future__ import absolute_import
from __future__ import print_function

import numpy as np
from scipy.optimize import minimize, differential_evolution, leastsq

from stishovite_data import (
    common_data,
    P_for_CN,
    V_for_CN,
    T_for_CN,
    KNR_GPa,
    KNR_err_GPa,
    CN_GPa,
    CN_err_GPa,
)
from scipy.optimize import root_scalar
from stishovite_model import make_models, make_scalar_model, modify_Zhang_elasticity
from copy import deepcopy
from burnman import RelaxedSolution
import matplotlib.pyplot as plt
from stishovite_parameters import scalar_bounds, cell_bounds, elastic_bounds
from stishovite_parameters import scalar_args, cell_args, elastic_args
from stishovite_parameters import scalar_and_cell_args, all_args

data = common_data()

delta_ijk = np.einsum("ij, lk, l->ijk", np.eye(6), np.eye(6), np.ones(6))

min_misfit_scalar = [1.0e7]
min_misfit_cell = [1.0e7]
min_misfit_elastic = [1.0e7]
min_misfit_combined = [1.0e7]

plot = False


def transition_pressure(solution, T):
    def diff_mu(P_GPa):
        solution.set_state(P_GPa * 1.0e9, T)
        solution.set_composition([0.499, 0.501])
        mu = solution.partial_gibbs
        return mu[1] - mu[0]

    sol = root_scalar(diff_mu, x0=80.0, x1=60.0)
    assert sol.converged
    return sol.root * 1.0e9


# We want to fit the Zhang room temperature data (a, b, c)
# and also the symmetry breaking pressure from Fischer at 3000 K
# (approximately 88 GPa)


def misfit_scalar(args):
    (
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        dDebye_0,
        P_tr_GPa,
        fP_Zhang,
        fP_Andrault,
        fP_Wang,
        fP_Nishihara,
        fP2_Zhang,
        fP2_Andrault,
        fP2_Wang,
        fP2_Nishihara,
    ) = args

    scalar_prms = make_scalar_model(
        dVQ0, dKQ0, dKpQ0, dgrQ0, dqQ0, V0Q1overV0Q0, dDebye_0, P_tr_GPa
    )
    stishovite_Q0, stishovite_Q1, scalar_stv, ESV_interactions = scalar_prms

    if plot:
        fig = plt.figure()
        ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

    relaxed_stv = RelaxedSolution(
        scalar_stv,
        relaxation_vectors=np.array([[1.0, -1.0]]),
        unrelaxed_vectors=np.array([[0.5, 0.5]]),
    )
    relaxed_stv.set_composition([1.0])

    # Fit the HT transition pressure from Fischer
    P_tr_3000K_model = transition_pressure(scalar_stv, 3000.0) / 1.0e9
    P_tr_3000K_obs = 90.0
    P_tr_3000K_err = 1.0
    misfits = [(P_tr_3000K_model - P_tr_3000K_obs) / P_tr_3000K_err]

    # Fit volumes for stishovite from Zhang, Andrault
    for phase, publication in [
        ("stv", "Nishihara_2005"),
        ("stv", "Wang_2012"),
        ("stv", "Zhang_2021"),
        ("poststv", "Zhang_2021"),
        ("stv", "Andrault_2003"),
        ("poststv", "Andrault_2003"),
    ]:
        if publication == "Zhang_2021":
            fP = fP_Zhang
            fP2 = fP2_Zhang
        elif publication == "Andrault_2003":
            fP = fP_Andrault
            fP2 = fP2_Andrault
        elif publication == "Wang_2012":
            fP = fP_Wang
            fP2 = fP2_Wang
        elif publication == "Nishihara_2005":
            fP = fP_Nishihara
            fP2 = fP2_Nishihara
        else:
            fP = 1.0
            fP2 = 0.0

        PTV = deepcopy(data["PTV"][phase][publication])
        PTV_err = deepcopy(data["PTV_err"][phase][publication])
        PTV_err[:, 0] = np.max([1.0e8 + 0.01 * PTV[:, 0], PTV_err[:, 0]])
        P_model = relaxed_stv.evaluate_with_volumes(["pressure"], PTV[:, 2], PTV[:, 1])[
            0
        ]
        P_actual = PTV[:, 0] * (fP + fP2 * PTV[:, 0])
        misfits.extend(list((P_model - P_actual) / PTV_err[:, 0]))

        if plot:
            sort = np.argsort(P_model)
            ax[0].plot(P_model[sort], PTV[:, 2][sort])
            ax[0].scatter(P_actual[sort], PTV[:, 2][sort])

    # Fit relative volumes from Ito
    PTV = deepcopy(data["PTV"]["stv"]["Ito_1974"])
    PTV_err = deepcopy(data["PTV_err"]["stv"]["Ito_1974"])
    scalar_stv.set_composition([0.5, 0.5])
    scalar_stv.set_state(1.0e5, 298.15)
    PTV[:, 2] = PTV[:, 2] / PTV[0, 2] * scalar_stv.V
    V_model = relaxed_stv.evaluate(["V"], PTV[:, 0], PTV[:, 1])[0]
    misfits.extend(list((V_model - PTV[:, 2]) / PTV_err[:, 2]))

    if plot:
        sort = np.argsort(PTV[:, 1])
        ax[1].plot(PTV[:, 1][sort], V_model[sort])
        ax[1].scatter(PTV[:, 1][sort], PTV[:, 2][sort])

    # Fit isentropic bulk moduli from Zhang
    KNR_obs = relaxed_stv.evaluate_with_volumes(
        ["isentropic_bulk_modulus_reuss"], V_for_CN, T_for_CN
    )[0]
    misfits.extend(list((KNR_GPa - KNR_obs / 1.0e9) / KNR_err_GPa))

    if plot:
        sort = np.argsort(V_for_CN)
        ax[2].plot(V_for_CN[sort], KNR_obs[sort] / 1.0e9)
        ax[2].scatter(V_for_CN[sort], KNR_GPa[sort])
        plt.show()

    misfits = np.array(misfits)
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_scalar[0]:
        min_misfit_scalar[0] = chisqr
        print(repr(args))
    print(chisqr)
    return misfits


def misfit_cell(args, scalar_args):
    (
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        dDebye_0,
        P_tr_GPa,
        fP_Zhang,
        fP_Andrault,
        fP_Wang,
        fP_Nishihara,
        fP2_Zhang,
        fP2_Andrault,
        fP2_Wang,
        fP2_Nishihara,
    ) = scalar_args
    (a0Q1, b0Q1, PsiI_33_a, PsiI_33_b, PsiI_33_c, PsiI_33_b2, PsiI_33_c2) = args
    (
        a11,
        a22,
        a33,
        a44,
        a55,
        a66,
        b11,
        b22,
        b33,
        b44,
        b55,
        b66,
        c44,
        c55,
        c66,
        b112,
        b222,
        b332,
    ) = np.zeros(18)

    models = make_models(
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        a0Q1,
        b0Q1,
        dDebye_0,
        P_tr_GPa,
        PsiI_33_a,
        PsiI_33_b,
        PsiI_33_c,
        PsiI_33_b2,
        PsiI_33_c2,
        a11,
        a22,
        a33,
        a44,
        a55,
        a66,
        b11,
        b22,
        b33,
        b44,
        b55,
        b66,
        c44,
        c55,
        c66,
        b112,
        b222,
        b332,
    )
    _, _, _, stishovite_relaxed = models
    stishovite_relaxed.set_composition([1.0])

    misfits = []

    if plot:
        fig = plt.figure()
        ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    # fit a, b and c for stishovite from Zhang and Andrault
    for phase, publication in [
        ("stv", "Zhang_2021"),
        ("poststv", "Zhang_2021"),
        # ("stv", "Andrault_2003"),
        ("poststv", "Andrault_2003"),
    ]:
        if publication == "Zhang_2021":
            fP = fP_Zhang
            fP2 = fP2_Zhang
        elif publication == "Andrault_2003":
            fP = fP_Andrault
            fP2 = fP2_Andrault
        elif publication == "Wang_2012":
            fP = fP_Wang
            fP2 = fP2_Wang
        elif publication == "Nishihara_2005":
            fP = fP_Nishihara
            fP2 = fP2_Nishihara
        else:
            fP = 1.0
            fP2 = 0.0

        PTV = data["PTV"][phase][publication]
        abc_obs = data["abc"][phase][publication]
        abc_err_obs = np.max(
            [
                data["abc"][phase][publication] * 0.0005,
                data["abc_err"][phase][publication],
            ],
            axis=0,
        )

        P_actual = PTV[:, 0] * (fP + fP2 * PTV[:, 0])
        cell_parameters_model = stishovite_relaxed.evaluate(
            ["cell_parameters"], P_actual, PTV[:, 1]
        )[0]
        abc_model = cell_parameters_model[:, :3]

        misfits.extend(list(np.ravel((abc_obs - abc_model) / abc_err_obs)))
        if plot:
            for i in range(3):
                j = 1
                if i < 2:
                    j = 0
                sort = np.argsort(P_actual)
                ax[j].errorbar(
                    P_actual[sort],
                    abc_obs[sort, i],
                    yerr=abc_err_obs[sort, i],
                    ls="none",
                )

    if plot:
        pressures = np.linspace(1.0e5, 120.0e9, 101)
        temperatures = pressures * 0.0 + 298.15
        cell_parameters_model = stishovite_relaxed.evaluate(
            ["cell_parameters"], pressures, temperatures
        )[0]
        abc_model = cell_parameters_model[:, :3]
        for i in range(3):
            j = 1
            if i < 2:
                j = 0
            ax[j].plot(pressures, abc_model[:, i])
        plt.show()

    misfits = np.array(misfits)
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_cell[0]:
        min_misfit_cell[0] = chisqr
        print(repr(args))
    print(chisqr)
    return misfits


def misfit_elastic(args, scalar_args, cell_args):
    (
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        dDebye_0,
        P_tr_GPa,
        fP_Zhang,
        fP_Andrault,
        fP_Wang,
        fP_Nishihara,
        fP2_Zhang,
        fP2_Andrault,
        fP2_Wang,
        fP2_Nishihara,
    ) = scalar_args
    (a0Q1, b0Q1, PsiI_33_a, PsiI_33_b, PsiI_33_c, PsiI_33_b2, PsiI_33_c2) = cell_args
    (a11, a22, a33, a44, a55, a66, b11i, b33i, b44i, b66i, c44, c66, b112i, b332i) = (
        args
    )

    b55i = b44i
    b22i = b11i
    b222i = b112i
    c55 = c44

    frel = -0.14
    b11 = b11i / (PsiI_33_c * np.exp(PsiI_33_c * frel) - 1.0)
    b22 = b22i / (PsiI_33_c * np.exp(PsiI_33_c * frel) - 1.0)
    b33 = b33i / (PsiI_33_c * np.exp(PsiI_33_c * frel) - 1.0)
    b112 = b112i / (PsiI_33_c2 * np.exp(PsiI_33_c2 * frel) - 1.0)
    b222 = b222i / (PsiI_33_c2 * np.exp(PsiI_33_c2 * frel) - 1.0)
    b332 = b332i / (PsiI_33_c2 * np.exp(PsiI_33_c2 * frel) - 1.0)
    b44 = b44i / (c44 * np.exp(c44 * frel) - 1.0)
    b55 = b55i / (c55 * np.exp(c55 * frel) - 1.0)
    b66 = b66i / (c66 * np.exp(c66 * frel) - 1.0)

    models = make_models(
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        a0Q1,
        b0Q1,
        dDebye_0,
        P_tr_GPa,
        PsiI_33_a,
        PsiI_33_b,
        PsiI_33_c,
        PsiI_33_b2,
        PsiI_33_c2,
        a11,
        a22,
        a33,
        a44,
        a55,
        a66,
        b11,
        b22,
        b33,
        b44,
        b55,
        b66,
        c44,
        c55,
        c66,
        b112,
        b222,
        b332,
    )
    _, _, _, stishovite_relaxed = models
    stishovite_relaxed.set_composition([1.0])

    misfits = []

    if plot:
        fig = plt.figure()
        ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

    # Fit compliance data
    fP = fP_Zhang
    fP2 = fP2_Zhang
    P_actual = P_for_CN * (fP + fP2 * P_for_CN)
    CN_model = stishovite_relaxed.evaluate(
        ["isentropic_stiffness_tensor"], P_actual, T_for_CN
    )[0]
    sort = np.argsort(P_actual)
    for i, j in [
        (0, 0),
        (1, 1),
        (2, 2),
        (0, 1),
        (0, 2),
        (1, 2),
        (3, 3),
        (4, 4),
        (5, 5),
    ]:
        chis = (CN_GPa[:, i, j] - CN_model[:, i, j] / 1.0e9) / CN_err_GPa[:, i, j]
        misfits.extend(list(chis))

        if plot:
            axi = 2
            if i < 2 and j < 2:
                axi = 0
            elif i < 3:
                axi = 1

            lc = ax[axi].plot(
                P_actual[sort], CN_model[sort, i, j] / 1.0e9, label=f"{i+1}{j+1}"
            )
            ax[axi].errorbar(
                P_actual[sort],
                CN_GPa[sort, i, j],
                yerr=CN_err_GPa[sort, i, j],
                ls="none",
                c=lc[0].get_color(),
            )

    if plot:
        for axi in range(3):
            ax[axi].legend()
        plt.show()

    misfits = np.array(misfits)
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_elastic[0]:
        min_misfit_elastic[0] = chisqr
        print(repr(args))
    print(chisqr)
    return misfits


def misfit_scalar_and_cell(args):
    scalar_args = args[:16]
    cell_args = args[16:]
    misfits = np.concatenate(
        (misfit_scalar(scalar_args), misfit_cell(cell_args, scalar_args))
    )
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_combined[0]:
        min_misfit_combined[0] = chisqr
        print(repr(args))
    print(chisqr)
    return misfits


def misfit_cell_and_elastic(args, scalar_args):
    cell_args = args[:7]
    elastic_args = args[7:]
    misfits_1 = misfit_cell(cell_args, scalar_args)
    modify_Zhang_elasticity(scalar_args, cell_args, elastic_args)
    misfits_2 = misfit_elastic(elastic_args, scalar_args, cell_args)

    misfits = np.concatenate((misfits_1, misfits_2))
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_combined[0]:
        min_misfit_combined[0] = chisqr
        print(repr(args))
    print(chisqr)
    return misfits


def misfit_all(args):
    scalar_args = args[:16]
    cell_args = args[16:23]
    elastic_args = args[23:]
    misfits_1 = misfit_scalar(scalar_args)
    misfits_2 = misfit_cell(cell_args, scalar_args)
    modify_Zhang_elasticity(scalar_args, cell_args, elastic_args)
    misfits_3 = misfit_elastic(elastic_args, scalar_args, cell_args)

    misfits = np.concatenate((misfits_1, misfits_2, misfits_3))
    chisqr = np.sum(np.power(misfits, 2.0))

    if chisqr < min_misfit_combined[0]:
        min_misfit_combined[0] = chisqr
        print(repr(args[:16]))
        print(repr(args[16:23]))
        print(repr(args[23:]))
    print(chisqr)
    return misfits


if __name__ == "__main__":
    if False:
        print("Starting misfit_scalar")
        sol = leastsq(misfit_scalar, scalar_args, full_output=True)
        print(sol)

    if False:

        def misfit_de(cell_args):
            return misfit_cell(cell_args, scalar_args)

        sol = differential_evolution(func=misfit_de, bounds=cell_bounds, x0=cell_args)

        sol = minimize(
            misfit_cell,
            cell_args,
            bounds=cell_bounds,
            args=(scalar_args),
            method="Nelder-Mead",
        )

    if False:
        print("Starting misfit_scalar_and_cell")
        sol = leastsq(
            misfit_scalar_and_cell, scalar_and_cell_args, factor=0.01, full_output=True
        )  # factor = 0.01 definitely helps
        print(sol)

    if True:
        print("Starting misfit_all")
        sol = leastsq(
            misfit_all, all_args, factor=0.01, full_output=True
        )  # factor = 0.01 definitely helps

        # Estimating covariance. See here:
        # https://stackoverflow.com/q/14854339/6272561
        x = sol[0]
        residuals = misfit_all(x)
        n_data = len(residuals)
        n_params = len(x)
        residual_variance = np.sum(np.power(residuals, 2.0)) / (n_data - n_params)

        sig = np.sqrt(np.diag(sol[1] * residual_variance))
        np.savetxt("model_output/covariance_matrix.dat", sol[1] * residual_variance)
        print(repr(sig))
        for i in range(len(x)):
            print(f"{x[i]:.4e} +/- {sig[i]:.4e}")

    if False:

        def lsq(args):
            return np.sum(np.power(misfit_all(args), 2.0))

        sol = minimize(
            misfit_all,
            all_args,
            bounds=scalar_bounds + cell_bounds + elastic_bounds,
            method="Nelder-Mead",
        )

    if False:
        cell_and_elastic_args = np.concatenate((cell_args, elastic_args))
        modify_Zhang_elasticity(scalar_args, cell_args, elastic_args)
        sol = minimize(
            misfit_cell_and_elastic,
            cell_and_elastic_args,
            bounds=cell_bounds + elastic_bounds,
            args=(scalar_args),
            method="Nelder-Mead",
        )

    if False:
        sol = minimize(
            misfit_elastic,
            elastic_args,
            bounds=elastic_bounds,
            args=(scalar_args, cell_args),
            method="Powell",
        )

    if False:
        sol = minimize(
            misfit_cell,
            cell_args,
            bounds=cell_bounds,
            args=(scalar_args),
            method="Nelder-Mead",
        )
