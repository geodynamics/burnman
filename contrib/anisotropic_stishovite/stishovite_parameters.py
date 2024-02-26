import numpy as np

# dVQ0, dKQ0, dKpQ0, dgrQ0, dqQ0,
# V0Q1overV0Q0, dDebye_0, P_tr_GPa,
# fP_Zhang, fP_Andrault, fP2_Zhang, fP2_Andrault
# 74
scalar_args = [
    -5.28259962e-03,
    -2.78232366e00,
    1.86241455e-03,
    -2.16473596e-01,
    -8.76004122e-02,
    9.93918688e-01,
    1.76286306e01,
    4.97878625e01,
    9.30217075e-01,
    9.93494742e-01,
    -4.54240420e-13,
    -9.20410207e-13,
]


# a0Q1, b0Q1,
# PsiI_33_a, PsiI_33_b, PsiI_33_c,
# PsiI_33_b2, PsiI_33_c2
# 1151
cell_args = [
    2.72829954e-02,
    2.85720276e-02,
    2.10496560e-01,
    -6.02250364e-01,
    1.11223359e00,
    8.20161603e-05,
    -1.04689736e01,
]


# a11, a22, a33, a44, a55, a66
# b11, b33, b44, b66
# c44, c66
# b112, b332

# let the b parameters be equal to
# frel = -0.14
# b1_input = b1 * c1 * exp(c1 * frel) - 1)
# 178
elastic_args = [
    3.21306666e-01,
    8.34212845e-01,
    4.93022578e-01,
    1.12011984e00,
    1.20957700e00,
    9.49364288e-01,
    -4.47565581e-02,
    7.06681926e-04,
    3.85553457e-01,
    1.13867050e-01,
    9.74629124e-01,
    9.98028989e-01,
    5.64360733e-01,
    2.12770903e-01,
]

scalar_and_cell_args = np.concatenate((scalar_args, cell_args))
all_args = np.concatenate((scalar_args, cell_args, elastic_args))

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
