import burnman
from burnman.classes.solutionmodel import PolynomialSolution
from burnman import AnisotropicMineral, AnisotropicSolution
from burnman import RelaxedAnisotropicSolution
import numpy as np
from copy import copy, deepcopy
from stishovite_parameters import all_args
from stishovite_data import P_for_CN, T_for_CN, SN_invGPa, CN_GPa
from tabulate import tabulate

stv_SLB = burnman.minerals.SLB_2022.st()
stv_SLB.property_modifiers = []


def modify_Zhang_elasticity(scalar_args, cell_args, elastic_args):
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
        fP2_Zhang,
        fP2_Andrault,
    ) = scalar_args
    (a0Q1, b0Q1, PsiI_33_a, PsiI_33_b, PsiI_33_c, PsiI_33_b2, PsiI_33_c2) = cell_args
    (a11, a22, a33, a44, a55, a66, b11i, b33i, b44i, b66i, c44, c66, b112i, b332i) = (
        elastic_args
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

    # Fit compliance data
    fP = fP_Zhang
    fP2 = fP2_Zhang
    P_actual = P_for_CN * (fP + fP2 * P_for_CN)

    beta_T_model = stishovite_relaxed.evaluate(
        ["isentropic_compressibility_tensor"], P_actual, T_for_CN
    )[0]
    beta_ii_model = np.einsum("ijj->ij", beta_T_model)
    SN = deepcopy(SN_invGPa) / 1.0e9
    beta_ii_obs = np.einsum("ijk->ij", SN[:, :3, :3])
    beta_RT_obs = np.einsum("ijk->i", SN[:, :3, :3])

    if False:
        # Modify in an "L" shape (C13, C23, C33, C32, C31)
        f = beta_ii_model[:, 2] / beta_ii_obs[:, 2]
        g_obs = 2.0 * np.sum(SN[:, :3, 2], axis=1) - SN[:, 2, 2]
        SN[:, 2, :] = np.einsum("ij, i->ij", SN[:, 2, :], f)
        SN[:, :, 2] = np.einsum("ij, i->ij", SN[:, :, 2], f)
        SN[:, 2, 2] = SN[:, 2, 2] / f

        f2 = (beta_RT_obs - g_obs * f) / (beta_RT_obs - g_obs)
        SN[:, :2, :2] = np.einsum("ijk, i->ijk", SN[:, :2, :2], f2)
    elif False:
        # Modify only the diagonals
        for i in range(3):
            SN[:, i, i] = SN[:, i, i] + beta_ii_model[:, i] - beta_ii_obs[:, i]
    else:
        # Modify only the off-diagonals
        dbeta = beta_ii_model - beta_ii_obs
        SN[:, 0, 1] = SN[:, 0, 1] + (dbeta[:, 0] + dbeta[:, 1] - dbeta[:, 2]) / 2.0
        SN[:, 0, 2] = SN[:, 0, 2] + (dbeta[:, 0] + dbeta[:, 2] - dbeta[:, 1]) / 2.0
        SN[:, 1, 2] = SN[:, 1, 2] + (dbeta[:, 1] + dbeta[:, 2] - dbeta[:, 0]) / 2.0
        SN[:, 1, 0] = SN[:, 0, 1]
        SN[:, 2, 0] = SN[:, 0, 2]
        SN[:, 2, 1] = SN[:, 1, 2]

    # modify in place
    # SN is used to calculate, so overwriting repeatedly is ok
    CN_GPa[:, :, :] = np.linalg.inv(SN) / 1.0e9


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
        b1 = params["b1"]
        b2 = params["b2"]
        dPsidf = (
            params["a"]
            + b1 * params["c1"] * (np.exp(params["c1"] * f) - 1.0)
            + b2 * params["c2"] * (np.exp(params["c2"] * f) - 1.0)
        )
        Psi = (
            0.0
            + (params["a"] - b1 * params["c1"] - b2 * params["c2"]) * f
            + b1 * (np.exp(params["c1"] * f) - 1.0)
            + b2 * (np.exp(params["c2"] * f) - 1.0)
        )
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
        "b1": np.zeros((6, 6)),
        "c1": np.zeros((6, 6)),
        "b2": np.zeros((6, 6)),
        "c2": np.zeros((6, 6)),
    }

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
        "b1": PsiI_33_b,
        "c1": PsiI_33_c,
        "b2": PsiI_33_b2,
        "c2": PsiI_33_c2,
    }
    # The following assumes that the compression of the a and
    # b axes is the same
    PsiI_11 = {
        "a": 0.5 * (1.0 - PsiI_33_a),
        "b1": -0.5 * PsiI_33_b,
        "c1": PsiI_33_c,
        "b2": -0.5 * PsiI_33_b2,
        "c2": PsiI_33_c2,
    }
    PsiI_22 = {
        "a": 0.5 * (1.0 - PsiI_33_a),
        "b1": -0.5 * PsiI_33_b,
        "c1": PsiI_33_c,
        "b2": -0.5 * PsiI_33_b2,
        "c2": PsiI_33_c2,
    }

    anisotropic_parameters["a"][0, 0] = a11
    anisotropic_parameters["b1"][0, 0] = b11
    anisotropic_parameters["c1"][0, 0] = PsiI_33_c
    anisotropic_parameters["b2"][0, 0] = b112
    anisotropic_parameters["c2"][0, 0] = PsiI_33_c2
    anisotropic_parameters["a"][1, 1] = a22
    anisotropic_parameters["b1"][1, 1] = b22
    anisotropic_parameters["c1"][1, 1] = PsiI_33_c
    anisotropic_parameters["b2"][1, 1] = b222
    anisotropic_parameters["c2"][1, 1] = PsiI_33_c2
    anisotropic_parameters["a"][2, 2] = a33
    anisotropic_parameters["b1"][2, 2] = b33
    anisotropic_parameters["c1"][2, 2] = PsiI_33_c
    anisotropic_parameters["b2"][2, 2] = b332
    anisotropic_parameters["c2"][2, 2] = PsiI_33_c2
    anisotropic_parameters["a"][3, 3] = a44
    anisotropic_parameters["b1"][3, 3] = b44
    anisotropic_parameters["c1"][3, 3] = c44
    anisotropic_parameters["a"][4, 4] = a55
    anisotropic_parameters["b1"][4, 4] = b55
    anisotropic_parameters["c1"][4, 4] = c55
    anisotropic_parameters["a"][5, 5] = a66
    anisotropic_parameters["b1"][5, 5] = b66
    anisotropic_parameters["c1"][5, 5] = c66

    # Fill the rest
    for p in ["a", "b1", "b2"]:

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

    for i, j in [[0, 1], [0, 2], [1, 2]]:
        anisotropic_parameters["c1"][i, j] = PsiI_33_c
        anisotropic_parameters["c1"][j, i] = PsiI_33_c
        anisotropic_parameters["c2"][i, j] = PsiI_33_c2
        anisotropic_parameters["c2"][j, i] = PsiI_33_c2

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
    for p in ["a", "b1", "c1", "b2", "c2"]:
        for i, j in [(0, 1), (3, 4)]:
            l1 = copy(anti_anisotropic_parameters[p][j, :])
            l2 = copy(anti_anisotropic_parameters[p][i, :])
            anti_anisotropic_parameters[p][i, :] = l1
            anti_anisotropic_parameters[p][j, :] = l2
            l1 = copy(anti_anisotropic_parameters[p][:, j])
            l2 = copy(anti_anisotropic_parameters[p][:, i])
            anti_anisotropic_parameters[p][:, i] = l1
            anti_anisotropic_parameters[p][:, j] = l2

    #  print(np.sum(anti_anisotropic_parameters["a"][:3], axis=0)[:3])
    #  print(np.sum(anti_anisotropic_parameters["b1"][:3], axis=0)[:3])
    #  print(np.sum(anti_anisotropic_parameters["c1"][:3], axis=0)[:3])
    #  exit()
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


def get_models():

    scalar_args = all_args[:12]
    cell_args = all_args[12:19]
    elastic_args = all_args[19:]

    modify_Zhang_elasticity(scalar_args, cell_args, elastic_args)

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
        fP2_Zhang,
        fP2_Andrault,
    ) = scalar_args
    (a0Q1, b0Q1, PsiI_33_a, PsiI_33_b, PsiI_33_c, PsiI_33_b2, PsiI_33_c2) = cell_args
    (a11, a22, a33, a44, a55, a66, b11i, b33i, b44i, b66i, c44, c66, b112i, b332i) = (
        elastic_args
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

    # The following assumes that the compression of the a and
    # b axes is the same
    PsiI_11_a = 0.5 * (1.0 - PsiI_33_a)
    PsiI_11_b = -0.5 * PsiI_33_b
    PsiI_11_b2 = -0.5 * PsiI_33_b2

    sm = models[1]
    Q0, Q1 = sm.solution_model.interaction_endmembers
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
        table.append([pps[i], f"{Q0.params[p]:.6e}", f"{Q1.params[p]:.6e}"])

    table.append(["$a_0$", "-", a0Q1])
    table.append(["$b_0$", "-", b0Q1])

    print(tabulate(table, headers="firstrow", tablefmt="latex_raw", floatfmt=".6e"))

    print("b_xs", sm.b_xs)

    headers = [""]
    headers.extend([f"$\\Psi_{{{i+1}}}$" for i in range(3)])
    headers.extend([f"$\\Psi_{{{i+1}{i+1}}}$" for i in range(6)])
    table = [headers]
    table.extend(
        [
            ["$a$", PsiI_11_a, "as above", PsiI_33_a, a11, a22, a33, a44, a55, a66],
            [
                "$b_1$",
                PsiI_11_b,
                "as above",
                PsiI_33_b,
                b11,
                "as above",
                b33,
                b44,
                "as above",
                b66,
            ],
            [
                "$c_1$",
                PsiI_33_c,
                "as above",
                "as above",
                "as above",
                "as above",
                "as above",
                c44,
                "as above",
                c66,
            ],
            [
                "$b_2$",
                PsiI_11_b2,
                "as above",
                PsiI_33_b2,
                b112,
                "as above",
                b332,
                "-",
                "-",
                "-",
            ],
            [
                "$c_2$",
                PsiI_33_c2,
                "as above",
                "as above",
                "as above",
                "as above",
                "as above",
                "-",
                "-",
                "-",
            ],
        ]
    )

    table = [
        [f"{item:.6e}" if type(item) is np.float64 else item for item in row]
        for row in table
    ]
    print(
        tabulate(
            list(map(list, zip(*table))),
            headers="firstrow",
            tablefmt="latex_raw",
            floatfmt=".6e",
        )
    )

    print("Corrections to pressure:")
    print(f"    Zhang: {fP_Zhang}, {fP2_Zhang}")
    print(f"    Andrault: {fP_Andrault}, {fP2_Andrault}")
    return models, fP_Zhang, fP_Andrault, fP2_Zhang, fP2_Andrault
