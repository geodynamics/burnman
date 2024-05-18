import numpy as np
import matplotlib.pyplot as plt
from plagioclase_parameters import scalar_args, cell_args, elastic_args
from plagioclase_model import make_anisotropic_model
from plagioclase_data import get_data
from burnman.minerals.SLB_2022 import plagioclase

ss_SLB = plagioclase()

ss = make_anisotropic_model(scalar_args, cell_args, elastic_args)

# Model data
p_ans = np.linspace(0.0, 1.0, 101)
molar_fractions = np.array([p_ans, 1.0 - p_ans]).T
pressures = p_ans * 0.0 + 1.0e5
temperatures = p_ans * 0.0 + 300.0
prps = ss.evaluate(
    [
        "molar_volume",
        "isothermal_bulk_modulus_reuss",
        "isothermal_stiffness_tensor",
        "isentropic_stiffness_tensor",
        "isothermal_compliance_tensor",
        "isothermal_compressibility_tensor",
        "cell_parameters",
    ],
    pressures,
    temperatures,
    molar_fractions,
)

prps_SLB = ss_SLB.evaluate(
    ["molar_volume", "isothermal_bulk_modulus_reuss"],
    pressures,
    temperatures,
    molar_fractions,
)


labels = ["V", "KTR", "CT", "CN", "ST", "betaT", "cell"]
prps = {labels[i]: prps[i] for i in range(7)}
prps["psi"] = np.einsum("ijk, i->ijk", prps["ST"], prps["KTR"])

prps_SLB = {labels[i]: prps_SLB[i] for i in range(2)}

d = get_data()


# Scalar properties
fig = plt.figure(figsize=(8, 4))
ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]


mask2 = p_ans < 0.5

ln = ax[0].plot(p_ans[mask2], prps["V"][mask2] * 1.0e6, label="C$\\bar{1}$, this study")
ln = ax[0].plot(
    p_ans,
    prps_SLB["V"] * 1.0e6,
    linestyle="--",
    color=ln[0].get_color(),
    label="SLB2022",
)

ax[0].plot(
    p_ans[~mask2], prps["V"][~mask2] * 1.0e6, color=ln[0].get_color(), linestyle=":"
)

mask = d["cell"]["p_an"] < 0.5
ax[0].scatter(
    d["cell"]["p_an"][mask], d["cell"]["V"][mask] * 1.0e6, color=ln[0].get_color()
)
ax[0].errorbar(
    d["cell"]["p_an"][mask],
    d["cell"]["V"][mask] * 1.0e6,
    d["cell"]["V_err"][mask] * 1.0e6,
    color=ln[0].get_color(),
    ls="none",
)
ax[0].scatter(
    d["cell"]["p_an"][~mask],
    d["cell"]["V"][~mask] * 1.0e6,
    color=ln[0].get_color(),
    alpha=0.5,
)
ax[0].errorbar(
    d["cell"]["p_an"][~mask],
    d["cell"]["V"][~mask] * 1.0e6,
    d["cell"]["V_err"][~mask] * 1.0e6,
    color=ln[0].get_color(),
    alpha=0.5,
    ls="none",
)


ax[1].plot(p_ans[mask2], prps["KTR"][mask2] / 1.0e9, color=ln[0].get_color())
ax[1].plot(p_ans, prps_SLB["KTR"] / 1.0e9, linestyle="--", color=ln[0].get_color())
ax[1].plot(
    p_ans[~mask2], prps["KTR"][~mask2] / 1.0e9, color=ln[0].get_color(), linestyle=":"
)

mask = d["beta"]["p_an"] < 0.47
ax[1].scatter(
    d["beta"]["p_an"][mask],
    1.0 / d["beta"]["bTR"][mask] / 1.0e9,
    color=ln[0].get_color(),
)
ax[1].errorbar(
    d["beta"]["p_an"][mask],
    1.0 / d["beta"]["bTR"][mask] / 1.0e9,
    (d["beta"]["bTR_err"] / d["beta"]["bTR"] / d["beta"]["bTR"])[mask] / 1.0e9,
    color=ln[0].get_color(),
    ls="none",
)
ax[1].scatter(
    d["beta"]["p_an"][~mask],
    1.0 / d["beta"]["bTR"][~mask] / 1.0e9,
    color=ln[0].get_color(),
    alpha=0.5,
)
ax[1].errorbar(
    d["beta"]["p_an"][~mask],
    1.0 / d["beta"]["bTR"][~mask] / 1.0e9,
    (d["beta"]["bTR_err"] / d["beta"]["bTR"] / d["beta"]["bTR"])[~mask] / 1.0e9,
    color=ln[0].get_color(),
    alpha=0.5,
    ls="none",
)

for i in range(2):
    ax[i].set_xlabel("$p_{an}$")

ax[0].set_ylabel("$V$ (cm$^3$/mol)")
ax[1].set_ylabel("$K_{\\text{TR}}$ (GPa)")

ax[0].legend()

fig.set_tight_layout(True)
fig.savefig("plag_V_KT.pdf")

# Elastic stiffness data
ijs = np.array(
    [
        [1, 1],
        [2, 2],
        [3, 3],
        [4, 4],
        [5, 5],
        [6, 6],
        [1, 2],
        [1, 3],
        [2, 3],
        [1, 5],
        [2, 5],
        [3, 5],
        [4, 6],
        [1, 4],
        [1, 6],
        [2, 4],
        [2, 6],
        [3, 4],
        [3, 6],
        [4, 5],
        [5, 6],
    ]
)

inds = ijs - 1
mask = d["CN"]["p_an"] < 0.5
mask2 = p_ans < 0.5


fig = plt.figure(figsize=(10, 8))
ax = [fig.add_subplot(3, 3, i) for i in range(1, 8)]
for irow, (i, j) in enumerate(inds):
    name = f"$C_{{{i+1}{j+1}}}$"
    axi = int((irow - (irow) % 3) / 3)
    ln = ax[axi].plot(p_ans[mask2], prps["CN"][mask2, i, j] / 1.0e9)
    ax[axi].plot(
        p_ans[~mask2],
        prps["CN"][~mask2, i, j] / 1.0e9,
        linestyle=":",
        color=ln[0].get_color(),
    )
    ax[axi].scatter(
        d["CN"]["p_an"][mask],
        d["CN"]["CN"][mask, i, j] / 1.0e9,
        label=name,
        color=ln[0].get_color(),
    )
    ax[axi].errorbar(
        d["CN"]["p_an"][mask],
        d["CN"]["CN"][mask, i, j] / 1.0e9,
        d["CN"]["CN_err"][mask, i, j] / 1.0e9,
        ls="none",
        color=ln[0].get_color(),
    )
    ax[axi].scatter(
        d["CN"]["p_an"][~mask],
        d["CN"]["CN"][~mask, i, j] / 1.0e9,
        color=ln[0].get_color(),
        alpha=0.5,
    )
    ax[axi].errorbar(
        d["CN"]["p_an"][~mask],
        d["CN"]["CN"][~mask, i, j] / 1.0e9,
        d["CN"]["CN_err"][~mask, i, j] / 1.0e9,
        ls="none",
        color=ln[0].get_color(),
        alpha=0.5,
    )

for i in range(7):
    ax[i].legend()
    ax[i].set_xlabel("$p_{an}$")
    ax[i].set_ylabel("$C_{Nij}$ (GPa)")

fig.set_tight_layout(True)
fig.savefig("plag_stiffnesses.pdf")

# psi
fig = plt.figure(figsize=(10, 8))
ax = [fig.add_subplot(3, 3, i) for i in range(1, 8)]
for irow, (i, j) in enumerate(inds):
    name = f"$S_{{N{i+1}{j+1}}} / \\beta_{{NR}}$"
    axi = int((irow - (irow) % 3) / 3)
    ln = ax[axi].plot(p_ans[mask2], prps["psi"][mask2, i, j])
    ax[axi].plot(
        p_ans[~mask2], prps["psi"][~mask2, i, j], linestyle=":", color=ln[0].get_color()
    )
    ax[axi].scatter(
        d["CN"]["p_an"][mask],
        d["CN"]["psiN"][mask, i, j],
        label=name,
        color=ln[0].get_color(),
    )
    ax[axi].scatter(
        d["CN"]["p_an"][~mask],
        d["CN"]["psiN"][~mask, i, j],
        color=ln[0].get_color(),
        alpha=0.5,
    )

for i in range(7):
    ax[i].legend()
    ax[i].set_xlabel("$p_{an}$")
    ax[i].set_ylabel("$S_{Nij} / \\beta_{NR}$")

fig.set_tight_layout(True)
fig.savefig("plag_psi.pdf")

# Cell parameters
fig = plt.figure(figsize=(10, 6))
ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]

labels = ["a", "b", "c", "alpha", "beta", "gamma"]

inds = ijs - 1
mask = d["cell"]["p_an"] < 0.5
mask2 = p_ans < 0.5
for i, label in enumerate(labels):

    ln = ax[i].plot(p_ans[mask2], prps["cell"][mask2, i])
    ax[i].plot(
        p_ans[~mask2], prps["cell"][~mask2, i], linestyle=":", color=ln[0].get_color()
    )
    ax[i].scatter(
        d["cell"]["p_an"][mask], d["cell"][label][mask], color=ln[0].get_color()
    )
    ax[i].errorbar(
        d["cell"]["p_an"][mask],
        d["cell"][label][mask],
        d["cell"][label + "_err"][mask],
        ls="none",
        color=ln[0].get_color(),
    )
    ax[i].scatter(
        d["cell"]["p_an"][~mask],
        d["cell"][label][~mask],
        color=ln[0].get_color(),
        alpha=0.5,
    )
    ax[i].errorbar(
        d["cell"]["p_an"][~mask],
        d["cell"][label][~mask],
        d["cell"][label + "_err"][~mask],
        ls="none",
        color=ln[0].get_color(),
        alpha=0.5,
    )

    if i < 3:
        ax[i].set_ylabel(f"${label}$ (m)")
    else:
        ax[i].set_ylabel(f"$\\{label}$ ($^{{\\circ}}$)")
    ax[i].set_xlabel("$p_{an}$")

fig.set_tight_layout(True)
fig.savefig("plag_cell_parameters.pdf")

# Compressibilities
fig = plt.figure(figsize=(10, 6))
ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]

ijs = np.array(
    [
        [1, 1],
        [2, 2],
        [3, 3],
        [2, 3],
        [1, 3],
        [1, 2],
    ]
)

inds = ijs - 1

mask = d["beta"]["p_an"] < 0.5
mask2 = p_ans < 0.5
for axi, (i, j) in enumerate(inds):

    ln = ax[axi].plot(p_ans[mask2], prps["betaT"][mask2, i, j] * 1.0e9)
    ax[axi].plot(
        p_ans[~mask2],
        prps["betaT"][~mask2, i, j] * 1.0e9,
        linestyle=":",
        color=ln[0].get_color(),
    )
    ax[axi].scatter(
        d["beta"]["p_an"][mask],
        d["beta"]["b"][mask, axi] * 1.0e9,
        color=ln[0].get_color(),
    )
    ax[axi].errorbar(
        d["beta"]["p_an"][mask],
        d["beta"]["b"][mask, axi] * 1.0e9,
        d["beta"]["b_err"][mask, axi] * 1.0e9,
        ls="none",
        color=ln[0].get_color(),
    )
    ax[axi].scatter(
        d["beta"]["p_an"][~mask],
        d["beta"]["b"][~mask, axi] * 1.0e9,
        color=ln[0].get_color(),
        alpha=0.5,
    )
    ax[axi].errorbar(
        d["beta"]["p_an"][~mask],
        d["beta"]["b"][~mask, axi] * 1.0e9,
        d["beta"]["b_err"][~mask, axi] * 1.0e9,
        ls="none",
        color=ln[0].get_color(),
        alpha=0.5,
    )

    ax[axi].set_ylabel(f"$\\beta_{{T{axi+1}}}$ (GPa$^{{-1}}$)")
    ax[axi].set_xlabel("$p_{an}$")

fig.set_tight_layout(True)
fig.savefig("plag_compressibilities.pdf")
plt.show()
