import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from stishovite_parameters import (
    all_args,
    scalar_param_names,
    cell_param_names,
    all_param_names,
)
import burnman
from tabulate import tabulate

popt = all_args
pcov = np.loadtxt("model_output/covariance_matrix.dat")

stv_SLB = burnman.minerals.SLB_2022.st()
stv_SLB.property_modifiers = []

popt[0] += stv_SLB.params["V_0"] * 1.0e6
popt[1] += stv_SLB.params["K_0"] / 1.0e9
popt[2] += stv_SLB.params["Kprime_0"]
popt[3] += stv_SLB.params["grueneisen_0"]
popt[4] += stv_SLB.params["q_0"]

sigma = np.sqrt(np.diag(pcov))
rel_error = [
    "{:g}".format(float("{:.1g}".format(i))) for i in sigma / np.abs(popt) * 100.0
]

table = [all_param_names, popt, sigma, rel_error]
table = list(zip(*table))
print(
    tabulate(
        table,
        headers=["Parameter", "Value", "Sigma", "Relative error (\\%, 1 s.f.)"],
        tablefmt="latex_raw",
        floatfmt=["s", ".6e", ".6e"],
    )
)


def correlation_from_covariance(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation


corr = np.around(correlation_from_covariance(pcov) * 100.0, 0)

d = pd.DataFrame(data=corr, columns=all_param_names, index=all_param_names)


sns.set_theme(style="white")

# Set up the matplotlib figure
fig, ax = plt.subplots(figsize=(11, 9))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(
    d,
    cmap=cmap,
    center=0,
    vmin=-100,
    annot=True,
    fmt="g",
    square=True,
    linewidths=0.5,
    annot_kws={"size": 6},
    cbar_kws={"shrink": 0.5, "label": "Correlation coefficient (%)"},
)

ax.hlines(
    [
        len(scalar_param_names) - 8,
        len(scalar_param_names),
        len(scalar_param_names + cell_param_names),
    ],
    *ax.get_xlim()
)
ax.vlines(
    [
        len(scalar_param_names) - 8,
        len(scalar_param_names),
        len(scalar_param_names + cell_param_names),
    ],
    *ax.get_ylim()
)


fig.set_tight_layout(True)
fig.savefig("figures/correlation_matrix.pdf")
plt.show()
