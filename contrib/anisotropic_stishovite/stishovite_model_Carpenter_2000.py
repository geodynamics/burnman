import numpy as np
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt
from burnman.minerals.HP_2011_ds62 import stv
from string import ascii_lowercase

# Expressions for the bare elastic constant values
# are taken from Table 2
# All units are GPa


def C11(P_GPa):
    return 578 + 5.38 * P_GPa


def C33(P_GPa):
    return 776 + 4.94 * P_GPa


def C12(P_GPa):
    return 86 + 5.38 * P_GPa


def C13(P_GPa):
    return 191 + 2.72 * P_GPa


def C44(P_GPa):
    return 252 + 1.88 * P_GPa


def C66(P_GPa):
    return 323 + 3.10 * P_GPa


# The following parameters are also taken from Table 2
Pastc = 49
Pc = 49.0 + 50.7

l1 = -8.0
l2 = 24.62
l3 = 17.0
a = -0.04856
b = 10.94

# l4 and l6 are taken from the text
l4 = 20.0
l6 = 20.0


def gibbs_over_tetragonal_volume_1(P_GPa):
    """
    The Gibbs energy divided by the
    volume of the (bare) tetragonal phase,
    as defined in Equation 1.
    """
    tet_mask = np.array([P < Pastc for P in P_GPa])
    gibbs_values = np.empty_like(P_GPa)
    gibbs_values[tet_mask] = 0.0

    P = P_GPa[~tet_mask]
    Q2 = Qsqr(P)
    Q = np.sqrt(Q2)
    Q4 = Q2 * Q2

    e1se2 = -l2 / (0.5 * (C11(P) - C12(P))) * Q
    e1pe2 = -l1 / (0.5 * (C11(P) + C12(P))) * Q2
    e1 = (e1se2 + e1pe2) / 2.0
    e2 = (-e1se2 + e1pe2) / 2.0
    e3 = -l3 / C33(P) * Q2

    gibbs_values[~tet_mask] = (
        0.5 * a * (P - Pc) * Q2
        + 0.25 * b * Q4
        + l1 * (e1 + e2) * Q2
        + l2 * (e1 - e2) * Q
        + l3 * e3 * Q2
        + 0.25 * (C11(P) + C12(P)) * e1pe2 * e1pe2
        + 0.25 * (C11(P) - C12(P)) * e1se2 * e1se2
        + C13(P) * e1pe2 * e3
        + 0.5 * C33(P) * e3 * e3
    )
    return gibbs_values


def gibbs_over_tetragonal_volume_4(P_GPa):
    """
    The Gibbs energy divided by the
    volume of the (bare) tetragonal phase,
    as defined in Equation 4.
    """
    tet_mask = np.array([P < Pastc for P in P_GPa])
    P = P_GPa[~tet_mask]
    gibbs_values = np.empty_like(P_GPa)
    gibbs_values[tet_mask] = 0.0
    gibbs_values[~tet_mask] = 0.5 * a * (P - Pastc) * Qsqr(P) + 0.25 * bast(P) * Qsqr(
        P
    ) * Qsqr(P)
    return gibbs_values


def Qsqr(P_GPa):
    """
    The value of Q squared, as defined in Equation 3
    and used in both Equations 1 and 4.
    """
    tet_mask = np.array([P < Pastc for P in P_GPa])
    Q2 = np.empty_like(P_GPa)
    Q2[tet_mask] = 0
    P = P_GPa[~tet_mask]
    Q2[~tet_mask] = a / bast(P) * (Pastc - P)
    return Q2


def bast(P_GPa):
    """
    The asterisked value of b, as given in
    Equation 5 and used in Equation 4.
    """
    return b - 2.0 * (
        (
            l3 * l3 * (C11(P_GPa) + C12(P_GPa))
            + 2.0 * l1 * l1 * C33(P_GPa)
            - 4.0 * l1 * l3 * C13(P_GPa)
        )
        / ((C11(P_GPa) + C12(P_GPa)) * C33(P_GPa) - 2.0 * C13(P_GPa) * C13(P_GPa))
    )


def chi(P_GPa):
    """
    The value of chi, as defined in
    Equations 9 and 10 and used to define the Cijs.
    """
    tet_mask = np.array([P < Pastc for P in P_GPa])
    chi_values = np.empty_like(P_GPa)
    chi_values[tet_mask] = 1.0 / (a * (P_GPa[tet_mask] - Pc))
    chi_values[~tet_mask] = 1.0 / (
        2.0 * a * b / bast(P_GPa[~tet_mask]) * (Pastc - P_GPa[~tet_mask])
        + a * (Pastc - Pc)
    )
    return chi_values


def CT(P_GPa, mode="eqm"):
    # Expressions in Table 1
    C = np.zeros((6, 6, len(P_GPa)))

    if mode == "eqm":
        X = chi(P_GPa)
        Q = np.sqrt(Qsqr(P_GPa))
    else:
        X = np.ones_like(P_GPa)
        Q = np.ones_like(P_GPa)

    tet_mask = np.array([(P < Pastc and mode == "eqm") or mode == "tet" for P in P_GPa])
    C[0, 0, tet_mask] = C11(P_GPa[tet_mask]) - l2 * l2 * X[tet_mask]
    C[1, 1, tet_mask] = C[0, 0, tet_mask]
    C[2, 2, tet_mask] = C33(P_GPa[tet_mask])
    C[0, 1, tet_mask] = C12(P_GPa[tet_mask]) + l2 * l2 * X[tet_mask]
    C[0, 2, tet_mask] = C13(P_GPa[tet_mask])
    C[1, 2, tet_mask] = C[0, 2, tet_mask]
    C[3, 3, tet_mask] = C44(P_GPa[tet_mask])
    C[4, 4, tet_mask] = C44(P_GPa[tet_mask])
    C[5, 5, tet_mask] = C66(P_GPa[tet_mask])

    Q = Q[~tet_mask]
    X = X[~tet_mask]
    C[0, 0, ~tet_mask] = (
        C11(P_GPa[~tet_mask])
        - (4.0 * l1 * l1 * Q * Q + l2 * l2 + 4.0 * l1 * l2 * Q) * X
    )
    C[1, 1, ~tet_mask] = (
        C11(P_GPa[~tet_mask])
        - (4.0 * l1 * l1 * Q * Q + l2 * l2 - 4.0 * l1 * l2 * Q) * X
    )
    C[2, 2, ~tet_mask] = C33(P_GPa[~tet_mask]) - 4.0 * l3 * l3 * Q * Q * X
    C[0, 1, ~tet_mask] = C12(P_GPa[~tet_mask]) - (4.0 * l1 * l1 * Q * Q - l2 * l2) * X
    C[0, 2, ~tet_mask] = (
        C13(P_GPa[~tet_mask]) - (4 * l1 * l3 * Q * Q + 2.0 * l2 * l3 * Q) * X
    )
    C[1, 2, ~tet_mask] = (
        C13(P_GPa[~tet_mask]) - (4 * l1 * l3 * Q * Q - 2.0 * l2 * l3 * Q) * X
    )
    C[3, 3, ~tet_mask] = C44(P_GPa[~tet_mask]) + 2.0 * l4 * Q
    C[4, 4, ~tet_mask] = C44(P_GPa[~tet_mask]) - 2.0 * l4 * Q
    C[5, 5, ~tet_mask] = C66(P_GPa[~tet_mask]) + 2.0 * l6 * Q * Q

    C[1, 0] = C[0, 1]
    C[2, 0] = C[0, 2]
    C[2, 1] = C[1, 2]

    return np.moveaxis(C, 2, 0)


def strains(P_GPa):
    tet_mask = np.array([P < Pastc for P in P_GPa])
    epsilons = np.zeros((3, len(P_GPa)))
    P = P_GPa[~tet_mask]
    Q2 = Qsqr(P)
    Q = np.sqrt(Q2)

    e1se2 = -l2 / (0.5 * (C11(P) - C12(P))) * Q
    e1pe2 = -l1 / (0.5 * (C11(P) + C12(P))) * Q2
    epsilons[0, ~tet_mask] = (e1se2 + e1pe2) / 2.0
    epsilons[1, ~tet_mask] = (-e1se2 + e1pe2) / 2.0
    epsilons[2, ~tet_mask] = -l3 / C33(P) * Q2
    return epsilons


pressures = np.linspace(0.0e9, 120.0e9, 1001)
temperatures = 300.0 + pressures * 0.0

Ceqm = CT(pressures / 1.0e9, mode="eqm") * 1.0e9
Ctet = CT(pressures / 1.0e9, mode="tet") * 1.0e9
Cort = CT(pressures / 1.0e9, mode="ort") * 1.0e9
ST = np.linalg.inv(Ceqm)
KTR = 1.0 / np.sum(ST[:, :3, :3], axis=(1, 2))
ST_tet = np.linalg.inv(Ctet)
KTR_tet = 1.0 / np.sum(ST_tet[:, :3, :3], axis=(1, 2))

eps = strains(pressures / 1.0e9)
plt.plot(pressures / 1.0e9, np.power(eps[0] - eps[1], 2.0))
plt.show()

for i, j in [(1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3), (2, 3)]:
    (ln,) = plt.plot(
        pressures / 1.0e9, Ceqm[:, i - 1, j - 1] / 1.0e9, label=f"$C_{{{i}{j}}}$"
    )
    # plt.plot(pressures, Cort[:, i, j],
    #         label=f'$C_{{{i}{j}}}$', linestyle=':', c=ln.get_color())

plt.plot(pressures / 1.0e9, KTR / 1.0e9, label="$K_{{TR}}$")
plt.legend()
plt.show()

V_0 = 1.401e-05
# V_0 = 1.401e-05 / 6.0223e23 * 1.e30 * 2.
lnV = cumulative_trapezoid(-1.0 / KTR, pressures, initial=0) + np.log(V_0)
V_orig = np.exp(lnV)

fig = plt.figure(figsize=(12, 4))
ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

G = gibbs_over_tetragonal_volume_1(pressures / 1.0e9) * V_orig * 1.0e9
ax[0].plot(pressures / 1.0e9, G, label="$\\mathcal{G}_{{eq 1}} \\cdot V_{{tet}}$")
G = gibbs_over_tetragonal_volume_4(pressures / 1.0e9) * V_orig * 1.0e9
ax[0].plot(pressures / 1.0e9, G, label="$\\mathcal{G}_{{eq 4}} \\cdot V_{{tet}}$")

ax[0].set_ylabel("Excess Gibbs energy ($\\mathcal{{G}}_{{xs}}$; J/mol)")

ax[1].plot(
    pressures / 1.0e9,
    V_orig * 1.0e6,
    linestyle=":",
    linewidth=3.0,
    label="tetragonal: from $V_0$ and $C_{{Tij0}}$s",
)

V_diff = np.gradient(G, pressures, edge_order=2)
V = V_orig + V_diff

ax[1].plot(
    pressures / 1.0e9,
    V * 1.0e6,
    linestyle="--",
    label="eqm: from $V_0 + C_{{Tij0}}$s + $\\partial \\mathcal{{G}}_{{xs}}/\\partial P$",
)

V_HP = stv().evaluate(["V"], pressures, temperatures)[0]
# ax[1].plot(pressures / 1.0e9, V_HP * 1.0e6, linestyle=":", label="eqm: V(HP)")

ax[1].set_ylabel("Volume (cm$^3$/mol)")

KTR2 = -np.gradient(pressures, np.log(V))

ax[2].plot(
    pressures / 1.0e9,
    KTR_tet / 1.0e9,
    linestyle=":",
    linewidth=3.0,
    label="tetragonal: from $C_{{Tij0}}$s",
)
ax[2].plot(
    pressures / 1.0e9, KTR / 1.0e9, label="eqm: from eqm $C_{{Tij}}$s", c="black"
)
ax[2].plot(
    pressures / 1.0e9,
    KTR2 / 1.0e9,
    linestyle="--",
    label="eqm: from $-V (\\partial P / \\partial V)_T$",
)

ax[2].set_ylabel("Isothermal Reuss bulk modulus (GPa)")

for i in range(3):
    ax[i].set_title(f"({ascii_lowercase[i]})")
    ax[i].set_xlabel("Pressure (GPa)")
    ax[i].legend()

fig.set_layout_engine("tight")
fig.savefig("figures/Carpenter_2000_model.pdf")
plt.show()
