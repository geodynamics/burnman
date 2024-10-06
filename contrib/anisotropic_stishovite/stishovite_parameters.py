import numpy as np

scalar_param_names = [
    "$V_0 (Q=0)$",
    "$K_0 (Q=0)$",
    "$K'_0 (Q=0)$",
    "$\\gamma_0 (Q=0)$",
    "$q_0 (Q=0)$",
    "$V_0 (Q=1) / V_0 (Q=0)$",
    "$\\Theta_0 (Q=1) - \\Theta_0 (Q=0)$",
    "$P_{{tr}} (T = T_{{ref}})$",
    "$f_{{P}}$ (Zhang)",
    "$f_{{P}}$ (Andrault)",
    "$f_{{P}}$ (Wang)",
    "$f_{{P}}$ (Nishihara)",
    "$f_{{P2}}$ (Zhang)",
    "$f_{{P2}}$ (Andrault)",
    "$f_{{P2}}$ (Wang)",
    "$f_{{P2}}$ (Nishihara)",
]

# 76
scalar_args = [
    -1.15714715e-02,
    -2.49554064e00,
    -7.69118482e-04,
    -2.08731542e-01,
    -6.08850674e-01,
    9.93834483e-01,
    1.76227286e01,
    4.99610469e01,
    8.91475407e-01,
    9.89762819e-01,
    9.87845187e-01,
    9.95720410e-01,
    2.92917650e-13,
    -9.08667999e-13,
    -2.20268875e-13,
    1.55926921e-12,
]


cell_param_names = [
    "$a_0 (Q=1)$",
    "$b_0 (Q=1)$",
    "$\\Psi_{{3a}}$",
    "$\\Psi_{{3b}}$",
    "$\\Psi_{{3c}}$",
    "$\\Psi_{{3b2}}$",
    "$\\Psi_{{3c2}}$",
]

# 1117
cell_args = [
    2.72779892e-02,
    2.85698440e-02,
    2.10385858e-01,
    -6.02813854e-01,
    1.10903361e00,
    8.86818639e-05,
    -1.09619314e01,
]


elastic_param_names = [
    "$\\Psi_{{11a}}$",
    "$\\Psi_{{22a}}$",
    "$\\Psi_{{33a}}$",
    "$\\Psi_{{44a}}$",
    "$\\Psi_{{55a}}$",
    "$\\Psi_{{66a}}$",
    "$\\Psi_{{11b}}$",
    "$\\Psi_{{33b}}$",
    "$\\Psi_{{44b}}$",
    "$\\Psi_{{66b}}$",
    "$\\Psi_{{44c}}$",
    "$\\Psi_{{66c}}$",
    "$\\Psi_{{11b2}}$",
    "$\\Psi_{{33b2}}$",
]

# let the b parameters be equal to
# frel = -0.14
# b1_input = b1 * c1 * exp(c1 * frel) - 1)
# 197
elastic_args = [
    3.24807928e-01,
    8.37663793e-01,
    4.93000185e-01,
    1.11959128e00,
    1.20914946e00,
    9.49237642e-01,
    -4.29973046e-02,
    7.10161257e-04,
    3.85475084e-01,
    1.13312302e-01,
    9.74598492e-01,
    9.97390953e-01,
    5.85741559e-01,
    2.13037695e-01,
]

scalar_and_cell_args = np.concatenate((scalar_args, cell_args))
all_args = np.concatenate((scalar_args, cell_args, elastic_args))
all_param_names = scalar_param_names + cell_param_names + elastic_param_names

scalar_bounds = (
    (-8.0e-3, -4.0e-3),  # dVQ0
    (-10.0, 0.0),  # dKQ0
    (-0.5, 2.0),  # dKpQ0  4.0292,
    (-0.6, 1.0),  # dgrQ0  1.55674,
    (-2.0, 3.0),  # dqQ0  2.2096,
    (0.99, 0.999),  # V0Q1overV0Q0
    (0.0, 80.0),  # dDebye_0
    (46, 53),  # P_tr_GPa,
    (0.9, 1.1),  # s
    (0.9, 1.1),  # s
    (0.9, 1.1),  # s
    (0.9, 1.1),  # s
    (-0.1, 0.1),  # s
    (-0.1, 0.1),
    (-0.1, 0.1),  # s
    (-0.1, 0.1),
)  # s

cell_bounds = (
    (2.6e-2, 2.8e-2),  # a0Q1
    (2.7e-2, 3.0e-2),  # b0Q1
    (0.0, 1.0),  # PsiI_33_a
    (-1.0, -0.1),  # PsiI_33_b
    (0.5, 2),  # PsiI_33_c
    (-1.0e-2, 1.0e-2),  # PsiI_33_b2
    (-20, -2),
)  # PsiI_33_c2

elastic_bounds = tuple((None, None) for i in range(14))
