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
]

cell_param_names = [
    "$a_0 (Q=1)$",
    "$b_0 (Q=1)$",
    "$\\Psi_{{3a}}$",
    "$\\Psi_{{3b}}$",
    "$\\Psi_{{3c}}$",
    "$\\Psi_{{3d}}$",
    "$f\\Psi_{{2}}$",
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
    "$\\Psi_{{44d}}$",
    "$\\Psi_{{66d}}$",
]

# let the b parameters be equal to
# frel = -0.14
# b1_input = b1 * c1 * exp(c1 * frel) - 1)
# 197

scalar_args = np.array(
    [
        6.97812599e-03,
        -3.18199894e00,
        2.02568387e-01,
        2.58439405e-01,
        -4.82534573e-05,
        9.97227979e-01,
        7.24798635e00,
        5.05759030e01,
    ]
)
cell_args = np.array(
    [
        0.02754879,
        0.02839605,
        0.19781129,
        -0.13650963,
        7.93215475,
        0.03176523,
        0.46676305,
    ]
)
elastic_args = np.array(
    [
        0.56655959,
        0.71555477,
        0.49784435,
        1.33005083,
        1.45154562,
        0.97764314,
        0.17631275,
        -0.0516391,
        -0.70121757,
        -0.18867047,
        3.6433892,
        4.63773363,
        0.31665864,
        0.09046939,
    ]
)


"""
# Removing the Fischer transition constraint
scalar_args = np.array([ 3.99584627e-03, -1.87172492e-01,  1.17470694e-01, -1.31671294e-01,
       -1.99994808e+00,  9.97305028e-01,  1.69419908e+01,  4.97481580e+01])
cell_args = np.array([ 0.02755226,  0.0284011 ,  0.1912243 , -0.14349027,  7.76739196,
        0.02422077,  0.46590064])
elastic_args = np.array([ 0.54560767,  0.66594729,  0.47720876,  1.33951091,  1.4542009 ,
        0.98069674,  0.2092898 , -0.07160769, -0.72256214, -0.16372882,
        3.56540327,  5.38951401,  0.32085463,  0.12083754])
"""


scalar_and_cell_args = np.concatenate((scalar_args, cell_args))
all_args = np.concatenate((scalar_args, cell_args, elastic_args))
all_param_names = scalar_param_names + cell_param_names + elastic_param_names

scalar_bounds = (
    (-8.0e-3, 7.0e-3),  # dVQ0
    (-10.0, 10.0),  # dKQ0
    (0.0, 2.0),  # dKpQ0  4.0292,
    (-0.6, 0.3),  # dgrQ0  1.55674,
    (-2.0, 0.0),  # dqQ0  2.2096,
    (0.99, 0.999),  # V0Q1overV0Q0
    (0.0, 100.0),  # dDebye_0
    (46, 53),  # P_tr_GPa,
)

cell_bounds = (
    (2.6e-2, 2.8e-2),  # a0Q1
    (2.7e-2, 3.0e-2),  # b0Q1
    (None, None),  # PsiI_33_a
    (None, None),  # PsiI_33_b
    (None, None),  # PsiI_33_c
    (None, None),  # PsiI_33_d
    (None, None),  # f
)

elastic_bounds = tuple((None, None) for i in range(14))


print(f"Total number of arguments: {len(all_args)}")
print(f"Scalar arguments: {len(scalar_args)}")
print(f"Cell arguments: {len(cell_args)}")
print(f"Elastic arguments: {len(elastic_args)}")
