import numpy as np
from plagioclase_model import make_scalar_model, make_anisotropic_model
from plagioclase_data import get_C1_data
from scipy.optimize import minimize
from plagioclase_parameters import scalar_args, cell_args, elastic_args
from scipy.optimize import differential_evolution

data = get_C1_data()

min_chisqr = [1.0e18]
iteration = [0]


def misfit_scalar(args):
    # For scalar EoS, fit isothermal compressibilities and volumes
    ss = make_scalar_model(args)

    molar_fractions = np.array([data["cell"]["p_an"], 1.0 - data["cell"]["p_an"]]).T
    pressures = data["cell"]["P"]
    temperatures = data["cell"]["T"]

    V_obs = ss.evaluate(["V"], pressures, temperatures, molar_fractions)[0]

    chisqr = np.sum(
        np.power((V_obs - data["cell"]["V"]) / (data["cell"]["V_err"]), 2.0)
    )

    molar_fractions = np.array([data["beta"]["p_an"], 1.0 - data["beta"]["p_an"]]).T
    pressures = data["beta"]["P"]
    temperatures = data["beta"]["T"]

    beta_model = ss.evaluate(
        ["isothermal_compressibility_reuss"], pressures, temperatures, molar_fractions
    )[0]
    chisqr += np.sum(
        np.power((beta_model - data["beta"]["bTR"]) / (data["beta"]["bTR_err"]), 2.0)
    )

    # add a chisqr for extreme nonlinear mixing behaviour in K_T
    # chisqr += np.power(args[5]/4., 2.)

    if chisqr < min_chisqr[0]:
        min_chisqr[0] = chisqr
        print(repr(args))
        print(chisqr)
    print()
    return chisqr


def misfit_cell(args, scalar_args, elastic_args):
    # For cell EoS, fit cell parameters and isothermal compressibilities

    cell_args = args
    ss = make_anisotropic_model(scalar_args, cell_args, elastic_args)

    # cell parameters
    pressures = data["cell"]["P"]
    temperatures = data["cell"]["T"]
    molar_fractions = np.array([data["cell"]["p_an"], 1.0 - data["cell"]["p_an"]]).T
    cell_model = ss.evaluate(
        ["cell_parameters"], pressures, temperatures, molar_fractions
    )[0]

    chisqr = 0.0
    f_cell = 0.1  # cell and compressibility uncertainties are poorly scaled
    for i, prm in enumerate(["a", "b", "c", "alpha", "beta", "gamma"]):
        chi = (cell_model[:, i] - data["cell"][prm]) / data["cell"][prm + "_err"]
        chisqr += f_cell * np.sum(np.power(chi, 2.0))

    # isothermal compressibilities
    if False:
        pressures = data["beta"]["P"]
        temperatures = data["beta"]["T"]
        molar_fractions = np.array([data["beta"]["p_an"], 1.0 - data["beta"]["p_an"]]).T
        beta_model = ss.evaluate(
            ["isothermal_compressibility_tensor"],
            pressures,
            temperatures,
            molar_fractions,
        )[0]

        f_beta = 10.0  # cell and compressibility uncertainties are poorly scaled
        for idx, (i, j) in enumerate([(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]):
            beta_model_i = beta_model[:, i, j]
            chi = (beta_model_i - data["beta"]["b"][:, idx]) / data["beta"]["b_err"][
                :, idx
            ]
            chisqr += f_beta * np.sum(np.power(chi, 2.0))

    if False:
        iteration[0] = iteration[0] + 1
        if chisqr < min_chisqr[0]:
            min_chisqr[0] = chisqr
            print(repr(args))
            print(chisqr)
        else:
            print(f"{iteration[0]}: {chisqr} \r")

    return chisqr


def misfit_elastic(args, scalar_args, cell_args):
    elastic_args = args
    ss = make_anisotropic_model(scalar_args, cell_args, elastic_args)

    pressures = data["CN"]["P"]
    temperatures = data["CN"]["T"]
    molar_fractions = np.array([data["CN"]["p_an"], 1.0 - data["CN"]["p_an"]]).T
    CN_model = ss.evaluate(
        ["isentropic_stiffness_tensor"], pressures, temperatures, molar_fractions
    )[0]

    chi = (CN_model - data["CN"]["CN"]) / data["CN"]["CN_err"]
    chisqr = np.sum(np.tri(6, 6) * np.power(chi, 2.0))

    if False:
        iteration[0] = iteration[0] + 1
        if chisqr < min_chisqr[0]:
            min_chisqr[0] = chisqr
            print(repr(args))
            print(chisqr)
        else:
            print(f"\r{iteration[0]}: {chisqr}", end="", flush=True)

    return chisqr


def misfit_cell_and_elastic(args, scalar_args):
    cell_args = args[:20]
    elastic_args = args[20:]
    chisqr = misfit_cell(cell_args, scalar_args, elastic_args)
    chisqr += misfit_elastic(elastic_args, scalar_args, cell_args)

    iteration[0] = iteration[0] + 1
    if chisqr < min_chisqr[0]:
        min_chisqr[0] = chisqr
        print(repr(args[:20]))
        print(repr(args[20:]))
        print(chisqr)
    else:
        print(f"\r{iteration[0]}: {chisqr}", end="", flush=True)

    return chisqr


def misfit_elastic_2(args, scalar_args, cell_args):
    elastic_args = args
    ss = make_anisotropic_model(scalar_args, cell_args, elastic_args)

    pressures = data["CN"]["P"]
    temperatures = data["CN"]["T"]
    molar_fractions = np.array([data["CN"]["p_an"], 1.0 - data["CN"]["p_an"]]).T
    SN_model, KNR_model = ss.evaluate(
        ["isentropic_compliance_tensor", "isentropic_bulk_modulus_reuss"],
        pressures,
        temperatures,
        molar_fractions,
    )
    phi_model = np.einsum("ijk, i->ijk", SN_model, KNR_model)

    chi = phi_model - data["CN"]["phiN"]
    chisqr = np.sum(np.power(chi, 2.0))

    iteration[0] = iteration[0] + 1
    if chisqr < min_chisqr[0]:
        min_chisqr[0] = chisqr
        print(repr(args))
        print(chisqr)
    else:
        print(f"\r{iteration[0]}: {chisqr}", end="", flush=True)

    return chisqr


if False:
    args = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    sol = minimize(misfit_scalar, args)
    exit()

if False:
    sol = minimize(
        misfit_cell,
        cell_args,
        args=(scalar_args, elastic_args),
    )

if False:
    bounds = [(-20, 20) for i in range(30)]
    sol = minimize(
        misfit_elastic,
        elastic_args,
        args=(scalar_args, cell_args),
    )  # bounds=bounds)

if False:
    bounds = [(-10, 10) for i in range(30)]
    sol = differential_evolution(
        misfit_elastic, bounds=bounds, x0=elastic_args, args=(scalar_args, cell_args)
    )

if True:
    ce_args = np.concatenate((cell_args, elastic_args))
    sol = minimize(misfit_cell_and_elastic, ce_args, args=(scalar_args))

    sol = minimize(
        misfit_cell_and_elastic, sol.x, args=(scalar_args), method="Nelder-Mead"
    )

    sol = minimize(
        misfit_cell_and_elastic, sol.x, args=(scalar_args), method="Nelder-Mead"
    )

    sol = minimize(misfit_cell_and_elastic, sol.x, args=(scalar_args))
