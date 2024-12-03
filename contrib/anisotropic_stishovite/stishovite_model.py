import burnman
from burnman.classes.solutionmodel import PolynomialSolution
from burnman import AnisotropicMineral, AnisotropicSolution
from burnman import RelaxedAnisotropicSolution
import numpy as np
from copy import copy, deepcopy
from tabulate import tabulate

stv_SLB = burnman.minerals.SLB_2022.st()
stv_SLB.property_modifiers = []


def make_scalar_model(dVQ0, dKQ0, dKpQ0, dgrQ0, dqQ0, V0Q1overV0Q0, dDebye_0, P_tr_GPa):
    # Make Q0
    stishovite_Q0 = deepcopy(stv_SLB)
    for prm, dv in [
        ("V_0", dVQ0 * 1.0e-6),  # cm^3
        ("K_0", dKQ0 * 1.0e9),  # GPa
        ("Kprime_0", dKpQ0),
        ("grueneisen_0", dgrQ0),
        ("q_0", dqQ0),
    ]:
        stishovite_Q0.params[prm] = stv_SLB.params[prm] + dv

    # Make Q1
    stishovite_Q1 = deepcopy(stv_SLB)
    Q_scale_factor = 6.25e10  # user-chosen variable (scales Q)
    V0Q1 = V0Q1overV0Q0 * stishovite_Q0.params["V_0"]
    dV = V0Q1 - stishovite_Q0.params["V_0"]
    dF = -dV * Q_scale_factor
    for prm, dv in [
        ("V_0", dV),
        ("K_0", 0.0),
        ("Kprime_0", 0.0),
        ("grueneisen_0", 0.0),
        ("q_0", 0.0),
        ("F_0", dF),
        ("Debye_0", dDebye_0),
    ]:
        stishovite_Q1.params[prm] = stishovite_Q0.params[prm] + dv

    """
    # (1^6 - Q^6) = 12*(pa^5*pb + pb^5*pa) + 40*(pa^3*pb^3)
    # (1^4 - Q^4) = 8*(pa^3*pb + pb^3*pa)
    # (1^2 - Q^2) = 4*(pa*pb)

    # G_xs = (a-b)*(1 - Q^2) + b*(1 - Q^n)

    # At Q = 0
    # G_xs = a
    # Let a = GQ0 - GQ1

    # At Q = 1
    # G_xs = 0

    # At Q = 0
    # d2G/dQdQ = -2*(a - b) - c*b*Q^(n-2)
    # a = b

    # Thus at Q0, the symmetry breaking occurs when b = (GQ0-GQ1)
    """

    T_tr = 298.15  # user-chosen variable
    P_tr = P_tr_GPa * 1.0e9
    stishovite_Q0.set_state(P_tr, T_tr)
    stishovite_Q1.set_state(P_tr, T_tr)
    b_xs = stishovite_Q0.gibbs - stishovite_Q1.gibbs

    # Q^4
    ESV_interactions = [
        [-4.0 * b_xs, 0.0, 0.0, 0, 1, 1, 1],
        [8.0 * b_xs, 0.0, 0.0, 0, 1, 1, 3],
        [8.0 * b_xs, 0.0, 0.0, 0, 3, 1, 1],
    ]

    """
    # Q^6
    ESV_interactions = [[-4.*b_xs, 0., 0., 0, 1, 1, 1],
                        [12.*b_xs, 0., 0., 0, 5, 1, 1],
                        [40.*b_xs, 0., 0., 0, 3, 1, 3],
                        [12.*b_xs, 0., 0., 0, 1, 1, 5]]
    """

    class stv_poststv(burnman.Solution):
        def __init__(self, molar_fractions=None):
            self.name = "stishovite-post-stishovite"
            self.solution_model = PolynomialSolution(
                endmembers=[[stishovite_Q1, "[Si]O2"], [stishovite_Q1, "[Si]O2"]],
                ESV_interactions=ESV_interactions,
                interaction_endmembers=[stishovite_Q0, stishovite_Q1],
                endmember_coefficients_and_interactions=[[4.0, -4.0, 0, 1, 1, 1]],
            )
            self.b_xs = b_xs  # for easy access later
            self.endmember_coefficients_and_interactions = [[4.0, -4.0, 0, 1, 1, 1]]
            burnman.Solution.__init__(self, molar_fractions=molar_fractions)

    return stishovite_Q0, stishovite_Q1, stv_poststv(), ESV_interactions


def make_models(
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
    PsiI_33_d,
    f_PsiI_22,
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
    d44,
    d55,
    d66,
):

    scalar_prms = make_scalar_model(
        dVQ0, dKQ0, dKpQ0, dgrQ0, dqQ0, V0Q1overV0Q0, dDebye_0, P_tr_GPa
    )
    stishovite_Q0, stishovite_Q1, scalar_stv, ESV_interactions = scalar_prms

    c0Q1 = stishovite_Q1.params["V_0"] / (a0Q1 * b0Q1)

    """
    ANISOTROPIC MODELS
    """

    def psi_func(f, Pth, params):
        a = params["a"]
        b = params["b"]
        c = params["c"]
        d = params["d"]

        dPsidf = a + b * np.tanh(c * f + d)
        Psi = a * f + b / c * np.log(np.cosh(c * f + d) / np.cosh(d))
        dPsidPth = np.zeros((6, 6))
        return (Psi, dPsidf, dPsidPth)

    # Make the ordered model
    # Note that c is not dependent on Q, so take the value of
    # c from the standard state volume of the Q1 structure
    # in normal Q0 stishovite (at room temperature):
    # c_Q1(V0, T0) =  c_Q0(V0(Q1), T0)

    # Note also that we assume that the Q-related strain is
    # not a function of V or T (note not in agreement with
    # Fischer data), so we can use the splitting
    # for the Q1-structure at any V or T where it is stable
    # to calculate the relative splitting of a and b at any
    # other state where Q is known.

    cell_parameters = np.array([a0Q1, b0Q1, c0Q1, 90.0, 90.0, 90.0])  # m, degrees
    anisotropic_parameters = {
        "a": np.zeros((6, 6)),
        "b": np.zeros((6, 6)),
        "c": np.ones((6, 6)) * PsiI_33_c,
        "d": np.ones((6, 6)) * PsiI_33_d,
    }

    anisotropic_parameters["c"][3][3] = c44
    anisotropic_parameters["c"][4][4] = c55
    anisotropic_parameters["c"][5][5] = c66
    anisotropic_parameters["d"][3][3] = d44
    anisotropic_parameters["d"][4][4] = d55
    anisotropic_parameters["d"][5][5] = d66
    """
    beta_T/beta_RT = d(Psi*I)/df
    F_ij = exp(Psi I)

    M=FM0

    for orthotropic systems:
    ln(M_ii) = ln(F_ii) + ln(M0_ii)
    ln(M_ii) - ln(M0_ii) = ln(F_ii) = (Psi I)_ii
    ln(a) - ln(a0) = (Psi I)_11
    ln(b) - ln(b0) = (Psi I)_22
    ln(c) - ln(c0) = (Psi I)_33
    ln(c/c0Q1) = aQ1 * fQ1 + bQ1 * (np.exp(cQ1 * fQ1) - 1.)

    note that bQ1 needs to be scaled relative to b fit using Q0,
    because of the change in volume:
    Let f(Q1) be the constant ln(V0Q1/V0Q0)
    fQ1 is ln(VQ1/V0Q1)
    ln(c) - ln(c0Q1) = a * (f - f(Q1)) + b * (np.exp(c * f) - np.exp(c * f(Q1)))
    ln(c) - ln(c0Q1) = a * fQ1 + b * (np.exp(c)*(np.exp(f) - np.exp(f(Q1))))
    ln(c) - ln(c0Q1) = a * fQ1 + b*np.exp(f(Q1)) * (np.exp(c)*(np.exp(f - f(Q1)) - 1))
    ln(c) - ln(c0Q1) = a * fQ1 + b*np.exp(f(Q1)) * (np.exp(c * fQ1) - 1)

    Thus
    aQ1 = a
    bQ1 = b*np.exp(f(Q1)) = b*V0Q1/V0Q0
    cQ1 = c

    # lna_Q1 = 0.5*(lnV - lnc_Q1) - lna_minus_lnb/2.
    # ln(a/a0)_Q1 = 0.5*(f - ln(c/c0)_Q1)
    """

    # the following give the a, b1, c1
    PsiI_33 = {
        "a": PsiI_33_a,
        "b": PsiI_33_b,
    }
    # The following assumes that the compression of the a and
    # b axes potentially differs
    PsiI_11 = {
        "a": (1.0 - f_PsiI_22) * (1.0 - PsiI_33_a),
        "b": -0.5 * PsiI_33_b,
    }
    PsiI_22 = {
        "a": f_PsiI_22 * (1.0 - PsiI_33_a),
        "b": -0.5 * PsiI_33_b,
    }

    anisotropic_parameters["a"][0, 0] = a11
    anisotropic_parameters["b"][0, 0] = b11
    anisotropic_parameters["a"][1, 1] = a22
    anisotropic_parameters["b"][1, 1] = b22
    anisotropic_parameters["a"][2, 2] = a33
    anisotropic_parameters["b"][2, 2] = b33
    anisotropic_parameters["a"][3, 3] = a44
    anisotropic_parameters["b"][3, 3] = b44
    anisotropic_parameters["a"][4, 4] = a55
    anisotropic_parameters["b"][4, 4] = b55
    anisotropic_parameters["a"][5, 5] = a66
    anisotropic_parameters["b"][5, 5] = b66

    # Fill the rest
    for p in ["a", "b"]:

        anisotropic_parameters[p][0, 1] = 0.5 * (
            (PsiI_11[p] + PsiI_22[p] - PsiI_33[p])
            - (
                anisotropic_parameters[p][0, 0]
                + anisotropic_parameters[p][1, 1]
                - anisotropic_parameters[p][2, 2]
            )
        )
        anisotropic_parameters[p][0, 2] = 0.5 * (
            (PsiI_11[p] + PsiI_33[p] - PsiI_22[p])
            - (
                anisotropic_parameters[p][0, 0]
                + anisotropic_parameters[p][2, 2]
                - anisotropic_parameters[p][1, 1]
            )
        )
        anisotropic_parameters[p][1, 2] = 0.5 * (
            (PsiI_22[p] + PsiI_33[p] - PsiI_11[p])
            - (
                anisotropic_parameters[p][1, 1]
                + anisotropic_parameters[p][2, 2]
                - anisotropic_parameters[p][0, 0]
            )
        )

        anisotropic_parameters[p][1, 0] = anisotropic_parameters[p][0, 1]
        anisotropic_parameters[p][2, 0] = anisotropic_parameters[p][0, 2]
        anisotropic_parameters[p][2, 1] = anisotropic_parameters[p][1, 2]

    # Make the antiordered model
    # Remember, Voigt form to standard notation is:
    # 11->1, 22->2, 33->3, 23->4, 13->5, 12->6
    # So ordered to antiordered transition through
    # a tetragonal phase involves the
    # following change of rows and columns
    # 1 <-> 2 (i.e. indices 0 and 1)
    # 4 <-> 5 (i.e. indices 3 and 4)

    anti_cell_parameters = deepcopy(cell_parameters)
    anti_cell_parameters[0] = cell_parameters[1]
    anti_cell_parameters[1] = cell_parameters[0]
    anti_anisotropic_parameters = deepcopy(anisotropic_parameters)
    for p in ["a", "b"]:
        for i, j in [(0, 1), (3, 4)]:
            l1 = copy(anti_anisotropic_parameters[p][j, :])
            l2 = copy(anti_anisotropic_parameters[p][i, :])
            anti_anisotropic_parameters[p][i, :] = l1
            anti_anisotropic_parameters[p][j, :] = l2
            l1 = copy(anti_anisotropic_parameters[p][:, j])
            l2 = copy(anti_anisotropic_parameters[p][:, i])
            anti_anisotropic_parameters[p][:, i] = l1
            anti_anisotropic_parameters[p][:, j] = l2

    anisotropic_stv_Q1 = AnisotropicMineral(
        stishovite_Q1,
        cell_parameters,
        anisotropic_parameters,
        psi_function=psi_func,
        orthotropic=True,
    )
    anisotropic_anti_stv_Q1 = AnisotropicMineral(
        stishovite_Q1,
        anti_cell_parameters,
        anti_anisotropic_parameters,
        psi_function=psi_func,
        orthotropic=True,
    )

    # An ideal mixing model
    def fn(lnV, Pth, X, params):
        z = np.zeros((6, 6))
        f = np.zeros((3, 3, 2))
        return (z, z, z, f)

    prm = {}

    stishovite_anisotropic = AnisotropicSolution(
        name="stishovite",
        solution_model=PolynomialSolution(
            endmembers=[
                [anisotropic_stv_Q1, "[Si]O2"],
                [anisotropic_anti_stv_Q1, "[Si]O2"],
            ],
            ESV_interactions=ESV_interactions,
            interaction_endmembers=[stishovite_Q0, stishovite_Q1],
            endmember_coefficients_and_interactions=[[4.0, -4.0, 0, 1, 1, 1]],
        ),
        psi_excess_function=fn,
        anisotropic_parameters=prm,
    )

    stishovite_relaxed = RelaxedAnisotropicSolution(
        stishovite_anisotropic,
        relaxation_vectors=np.array([[1.0, -1.0]]),
        unrelaxed_vectors=np.array([[0.5, 0.5]]),
    )

    return stishovite_Q1, scalar_stv, stishovite_anisotropic, stishovite_relaxed


def get_models(all_args):

    scalar_args = all_args[:8]
    cell_args = all_args[8:15]
    elastic_args = all_args[15:]

    (
        dVQ0,
        dKQ0,
        dKpQ0,
        dgrQ0,
        dqQ0,
        V0Q1overV0Q0,
        dDebye_0,
        P_tr_GPa,
    ) = scalar_args
    (a0Q1, b0Q1, PsiI_33_a, PsiI_33_b, PsiI_33_c, PsiI_33_d, f_PsiI_22) = cell_args
    (a11, a22, a33, a44, a55, a66, b11, b33, b44, b66, c44, c66, d44, d66) = (
        elastic_args
    )

    b22 = b11
    b55 = b44
    c55 = c44
    d55 = d44

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
        PsiI_33_d,
        f_PsiI_22,
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
        d44,
        d55,
        d66,
    )

    # The following assumes that the compression of the a and
    # b axes potentially differs
    PsiI_11_a = (1.0 - f_PsiI_22) * (1.0 - PsiI_33_a)
    PsiI_11_b = -0.5 * PsiI_33_b
    PsiI_22_a = f_PsiI_22 * (1.0 - PsiI_33_a)
    PsiI_22_b = -0.5 * PsiI_33_b

    sm = models[1]
    Q0, Q1 = sm.solution_model.interaction_endmembers
    Q1.params["F_0"] -= Q0.params["F_0"]
    Q0.params["F_0"] = 0
    properties = ["F_0", "V_0", "K_0", "Kprime_0", "Debye_0", "grueneisen_0", "q_0"]
    pps = [
        "$\\mathcal{F}_0$",
        "$V_0$",
        "$K_0$",
        "$K'_0$",
        "Debye $T_0$",
        "$\\gamma_0$",
        "$q_0$",
    ]
    table = [["", "$Q = 0$", "$Q = 1$"]]
    for i, p in enumerate(properties):
        if p in ["K_0", "Kprime_0", "grueneisen_0", "q_0"]:
            table.append([pps[i], f"{Q0.params[p]:.3e}", f"({Q1.params[p]:.3e})"])
        elif p == "Debye_0":
            table.append(
                [pps[i], f"({Q0.params[p]:.3e}; SLB22)", f"{Q1.params[p]:.3e}"]
            )
        else:
            table.append([pps[i], f"{Q0.params[p]:.3e}", f"{Q1.params[p]:.3e}"])

    table[1][1] = "-"
    table.append(["$a_0$", "-", a0Q1])
    table.append(["$b_0$", "-", b0Q1])

    print(tabulate(table, headers="firstrow", tablefmt="latex_raw", floatfmt=".3e"))

    print("b_xs", sm.b_xs)

    headers = [""]
    headers.extend([f"$\\Psi_{{{i+1}}}$" for i in range(3)])
    headers.extend([f"$\\Psi_{{{i+1}{i+1}}}$" for i in range(6)])
    table = [headers]
    table.extend(
        [
            [
                "$a$",
                PsiI_11_a,
                f"({PsiI_22_a:.3e})",
                PsiI_33_a,
                a11,
                a22,
                a33,
                a44,
                a55,
                a66,
            ],
            [
                "$b$",
                PsiI_11_b,
                f"({PsiI_22_b:.3e})",
                f"({PsiI_33_b:.3e})",
                b11,
                f"({b22:.3e})",
                b33,
                b44,
                f"({b55:.3e})",
                b66,
            ],
            [
                "$c$",
                PsiI_33_c,
                f"({PsiI_33_c:.3e})",
                f"({PsiI_33_c:.3e})",
                f"({PsiI_33_c:.3e})",
                f"({PsiI_33_c:.3e})",
                f"({PsiI_33_c:.3e})",
                c44,
                f"({c55:.3e})",
                c66,
            ],
            [
                "$d$",
                PsiI_33_d,
                f"({PsiI_33_d:.3e})",
                f"({PsiI_33_d:.3e})",
                f"({PsiI_33_d:.3e})",
                f"({PsiI_33_d:.3e})",
                f"({PsiI_33_d:.3e})",
                d44,
                f"({d55:.3e})",
                d66,
            ],
        ]
    )

    table = [
        [f"{item:.3e}" if type(item) is np.float64 else item for item in row]
        for row in table
    ]
    print(
        tabulate(
            list(map(list, zip(*table))),
            headers="firstrow",
            tablefmt="latex_raw",
            floatfmt=".3e",
        )
    )

    return models
