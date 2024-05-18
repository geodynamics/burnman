from burnman import Solution
from burnman.minerals import SLB_2022
from burnman.classes.solutionmodel import PolynomialSolution
from burnman import AnisotropicMineral, AnisotropicSolution
import numpy as np
from copy import deepcopy
from burnman.utils.unitcell import cell_parameters_to_vectors
from tabulate import tabulate


def make_scalar_model(args):
    dV_ab, dV_an, dV_aban = args[:3] / 1.0e6
    dK_ab, dK_an, dK_aban = args[3:6] * 1.0e9

    ab_scalar = SLB_2022.albite()
    an_scalar = SLB_2022.anorthite()

    ab_scalar.params["V_0"] = ab_scalar.params["V_0"] + dV_ab
    an_scalar.params["V_0"] = an_scalar.params["V_0"] + dV_an
    ab_scalar.params["K_0"] = ab_scalar.params["K_0"] + dK_ab
    an_scalar.params["K_0"] = an_scalar.params["K_0"] + dK_an

    aban_linear = SLB_2022.anorthite()
    aban_linear.params["V_0"] = 0.5 * (
        ab_scalar.params["V_0"] + an_scalar.params["V_0"]
    )
    aban_linear.params["K_0"] = 0.5 * (
        ab_scalar.params["K_0"] + an_scalar.params["K_0"]
    )
    aban_linear.params["Kprime_0"] = 0.5 * (
        ab_scalar.params["Kprime_0"] + an_scalar.params["Kprime_0"]
    )

    aban_scalar = deepcopy(aban_linear)

    aban_scalar.params["V_0"] = aban_scalar.params["V_0"] + dV_aban
    aban_scalar.params["K_0"] = aban_scalar.params["K_0"] + dK_aban

    table = [
        ["", "ab", "an", "an$_{50}$ (1)", "an$_{50}$ (2)"],
        [
            "V_0 (cm$^3$/mol)",
            ab_scalar.params["V_0"] * 1.0e6,
            an_scalar.params["V_0"] * 1.0e6,
            aban_linear.params["V_0"] * 1.0e6,
            aban_scalar.params["V_0"] * 1.0e6,
        ],
        [
            "K_0 (GPa)",
            ab_scalar.params["K_0"] / 1.0e9,
            an_scalar.params["K_0"] / 1.0e9,
            aban_linear.params["K_0"] / 1.0e9,
            aban_scalar.params["K_0"] / 1.0e9,
        ],
    ]
    print(tabulate(table, headers="firstrow", tablefmt="latex_raw", floatfmt=".6e"))

    class plagioclase_scalar(Solution):
        def __init__(self, molar_fractions=None):
            self.name = "plagioclase (NCAS)"
            self.solution_model = PolynomialSolution(
                endmembers=[
                    [an_scalar, "[Ca][Al]2Si2O8"],
                    [ab_scalar, "[Na][Al1/2Si1/2]2Si2O8"],
                ],
                interaction_endmembers=[aban_scalar, aban_linear],
                endmember_coefficients_and_interactions=[[4.0, -4.0, 0, 1, 1, 1]],
            )
            Solution.__init__(self, molar_fractions=molar_fractions)

    return plagioclase_scalar()


def make_anisotropic_model(scalar_args, cell_args, elastic_args):

    ss = make_scalar_model(scalar_args)
    an_scalar = ss.endmembers[0][0]
    ab_scalar = ss.endmembers[1][0]
    aban_scalar = ss.solution_model.interaction_endmembers[0]
    aban_linear = ss.solution_model.interaction_endmembers[1]

    b_over_a_ab, c_over_a_ab = cell_args[:2]
    alpha_ab, beta_ab, gamma_ab = cell_args[2:5]
    ab_cell_parameters = np.array(
        [1.0, b_over_a_ab, c_over_a_ab, alpha_ab, beta_ab, gamma_ab]
    )

    frame_convention = [1, 2, 0]

    M = cell_parameters_to_vectors(ab_cell_parameters, frame_convention)
    f = np.cbrt(ab_scalar.params["V_0"] / np.linalg.det(M))
    ab_cell_parameters[:3] = ab_cell_parameters[:3] * f

    b_over_a_an, c_over_a_an = cell_args[5:7]
    alpha_an, beta_an, gamma_an = cell_args[7:10]
    an_cell_parameters = np.array(
        [1.0, b_over_a_an, c_over_a_an, alpha_an, beta_an, gamma_an]
    )
    M = cell_parameters_to_vectors(an_cell_parameters, frame_convention)
    f = np.cbrt(an_scalar.params["V_0"] / np.linalg.det(M))
    an_cell_parameters[:3] = an_cell_parameters[:3] * f

    table = [
        ["", "$a$", "$b$", "$c$", "$\\alpha$", "$\\beta$", "$\\gamma$"],
        ["ab"],
        ["an"],
    ]
    table[1].extend(ab_cell_parameters)
    table[2].extend(an_cell_parameters)

    print(tabulate(table, headers="firstrow", tablefmt="latex_raw", floatfmt=".6e"))

    an_a = np.zeros((6, 6))
    ab_a = np.zeros((6, 6))

    # the following are all the elastic arguments
    # 15x2 = 30 in total
    # corresponding to all the block off-diagonals
    idx = 0
    # 3*2 = 6; 44, 55, 66
    for i in range(3, 6):
        an_a[i, i] = elastic_args[idx]
        ab_a[i, i] = elastic_args[idx + 1]
        idx = idx + 2
    # 6*2 = 12; 12, 13, 23, 15, 16, 26
    for j_shift in [0, 3]:
        for i in range(2):
            for j in range(i + 1, 3):
                an_a[i, j + j_shift] = elastic_args[idx]
                ab_a[i, j + j_shift] = elastic_args[idx + 1]
                an_a[j + j_shift, i] = elastic_args[idx]
                ab_a[j + j_shift, i] = elastic_args[idx + 1]
                idx = idx + 2

    # 6*2 = 12; 24, 34, 35, 45, 46, 56
    for i_shift in [0, 3]:
        for i in range(1, 3):
            for j in range(3, i + 3):
                an_a[i + i_shift, j] = elastic_args[idx]
                ab_a[i + i_shift, j] = elastic_args[idx + 1]
                an_a[j, i + i_shift] = elastic_args[idx]
                ab_a[j, i + i_shift] = elastic_args[idx + 1]
                idx = idx + 2

    # the following is another 5 x 2 = 10 cell arguments
    # corresponding to all the block diagonals
    # bringing the total to idx + 10 = 20
    an_b = np.sum(an_a[:3, :3], axis=0)
    an_c = np.sum(an_a[:3, 3:], axis=0)
    ab_b = np.sum(ab_a[:3, :3], axis=0)
    ab_c = np.sum(ab_a[:3, 3:], axis=0)

    idx = 10
    # 2x2 = 4; 22, 33
    for i in range(1, 3):
        j = i
        an_a[i, i] = cell_args[idx] - an_b[i]
        ab_a[i, i] = cell_args[idx + 1] - ab_b[i]
        idx = idx + 2
    # 3x2 = 6; 14, 25, 36
    for i in range(3):
        j = i + 3
        an_a[i, j] = cell_args[idx] - an_c[i]
        ab_a[i, j] = cell_args[idx + 1] - ab_c[i]
        an_a[j, i] = cell_args[idx] - an_c[i]
        ab_a[j, i] = cell_args[idx + 1] - ab_c[i]
        idx = idx + 2
    # Finally, 11
    an_a[0, 0] = 1.0 - np.sum(an_a[:3, :3])
    ab_a[0, 0] = 1.0 - np.sum(ab_a[:3, :3])

    table = ab_a
    print(tabulate(table, tablefmt="latex_raw", floatfmt=".5f"))

    table = an_a
    print(tabulate(table, tablefmt="latex_raw", floatfmt=".5f"))

    # print(an_a)
    # print(ab_a)
    # exit()

    an_anisotropic_parameters = {"a": an_a}
    ab_anisotropic_parameters = {"a": ab_a}

    def psi_func(f, Pth, params):
        dPsidf = params["a"]
        Psi = params["a"] * f
        dPsidPth = np.zeros((6, 6))
        return (Psi, dPsidf, dPsidPth)

    an_aniso = AnisotropicMineral(
        an_scalar,
        an_cell_parameters,
        an_anisotropic_parameters,
        psi_function=psi_func,
        frame_convention=frame_convention,
        orthotropic=False,
    )

    ab_aniso = AnisotropicMineral(
        ab_scalar,
        ab_cell_parameters,
        ab_anisotropic_parameters,
        psi_function=psi_func,
        frame_convention=frame_convention,
        orthotropic=False,
    )

    # An ideal mixing model with two endmembers
    def ideal_fn(lnV, Pth, X, params):
        z = np.zeros((6, 6))
        f = np.zeros((3, 3, 2))
        return (z, z, z, f)

    prm = {}

    plagioclase_anisotropic = AnisotropicSolution(
        name="plagioclase (C1)",
        solution_model=PolynomialSolution(
            endmembers=[
                [an_aniso, "[Ca][Al]2Si2O8"],
                [ab_aniso, "[Na][Al1/2Si1/2]2Si2O8"],
            ],
            interaction_endmembers=[aban_scalar, aban_linear],
            endmember_coefficients_and_interactions=[[4.0, -4.0, 0, 1, 1, 1]],
        ),
        psi_excess_function=ideal_fn,
        anisotropic_parameters=prm,
    )

    return plagioclase_anisotropic
