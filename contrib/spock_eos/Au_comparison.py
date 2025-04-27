from scipy.special import gamma, gammainc, exp1
import numpy as np
import matplotlib.pyplot as plt
from burnman.eos.reciprocal_kprime import RKprime
from burnman.eos.birch_murnaghan import BM3
from burnman.eos.birch_murnaghan_4th import BM4
from burnman.eos.vinet import Vinet
from burnman.eos.modified_tait import MT
from burnman.eos.spock import SPOCK
from burnman.eos.macaw import MACAW, make_params
from burnman import Mineral
from burnman.tools.eos import check_eos_consistency
import warnings

"""
Au_comparison
-------------
"""

# Colorblind friendly colors
# https://www.nature.com/articles/nmeth.1618.pdf
colors = {
    "SPOCK": "#000000",
    "BM3": "#E69F00",
    "none": "#56B4E9",
    "MT": "#009E73",
    "Vinet": "#F0E442",
    "MACAW": "#0072B2",
    "BM4": "#D55E00",
    "RK": "#CC79A7",
}


# Define the equations of state we wish to compare
eoses = {
    "Vinet": Vinet(),
    "BM3": BM3(),
    "BM4": BM4(),
    "MT": MT(),
    "RK": RKprime(),
    "MACAW": MACAW(),
    "SPOCK": SPOCK(),
}

# Create the parameter dictionary with the required values for
# all the equations of state.
params = {
    "E_0": 0.0,
    "P_0": 1.0e5,
    "T_0": 0.0,
    "V_0": 10.1959e-6,
    "K_0": 156.6e9,
    "Kprime_0": 6.76,
    "Kprime_prime_0": -15.6 / 156.6e9,
    "Kdprime_0": -15.6 / 156.6e9,
    "Kprime_inf": 2.97,
    "molar_mass": 0.06008,
    "equation_of_state": None,
}


# Define the range of volumes over which to
# compare the equations of state
volumes = np.linspace(1.0e-6, 1.4, 1001) * params["V_0"]

# Initialise the figure and subplots
fig = plt.figure(figsize=(10, 7))
ax = [fig.add_subplot(2, 2, i) for i in [4, 2, 1, 3]]

for name, eos in eoses.items():

    # Set the equation of state and initialise the mineral object
    params["equation_of_state"] = eos
    m = Mineral(params)

    # Check that the SPOCK and MACAW equations of state are internally consistent
    if name == "SPOCK" or name == "MACAW":
        assert check_eos_consistency(m, 1.0e10, 300.0, including_shear_properties=False)

    # Define some plotting values
    if name == "SPOCK":
        linestyle = "-"
        linewidth = 2.0
        alpha = 1.0
    else:
        linestyle = "-"
        linewidth = 1.5
        alpha = 1.0

    # Get the values of the properties to plot from the equations of state
    pressures = np.empty_like(volumes)
    KTs = np.empty_like(pressures)
    Vs = np.empty_like(pressures)
    Es = np.empty_like(pressures)
    for i, volume in enumerate(volumes):
        with warnings.catch_warnings(record=True) as w:
            try:
                Vs[i] = volumes[i]
                pressures[i] = eos.pressure(0.0, volumes[i], params)
                KTs[i] = eos.isothermal_bulk_modulus_reuss(
                    pressures[i], 0.0, volumes[i], params
                )
                Es[i] = (
                    eos.molar_internal_energy(pressures[i], 0.0, volumes[i], params)
                    / 1.0e9
                )
            except ValueError:
                pressures[i] = np.nan
                KTs[i] = np.nan
                Es[i] = np.nan

    with warnings.catch_warnings(record=True) as w:
        Kprimes = -np.gradient(np.log(KTs), np.log(Vs), edge_order=2)

    # Plot the properties
    ax[0].plot(
        Vs / params["V_0"],
        Kprimes,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=alpha,
        label=name,
        c=colors[name],
    )
    ax[1].plot(
        Vs / params["V_0"],
        KTs / 1.0e9,
        alpha=alpha,
        linestyle=linestyle,
        linewidth=linewidth,
        label=name,
        c=colors[name],
    )
    ax[2].plot(
        Vs / params["V_0"],
        pressures / 1.0e9,
        alpha=alpha,
        linestyle=linestyle,
        linewidth=linewidth,
        label=name,
        c=colors[name],
    )
    ax[3].plot(
        Vs / params["V_0"],
        Es * 1.0e9,
        alpha=alpha,
        linestyle=linestyle,
        linewidth=linewidth,
        label=name,
        c=colors[name],
    )


def Kdprime_0K_0_MACAW(params):
    A, B, C = make_params(params["K_0"], params["Kprime_0"], params["Kprime_inf"])
    return (
        C
        * (
            -3.0 * B**4
            - 12.0 * B**3 * C
            + 3.0 * B**3
            - 18.0 * B**2 * C**2
            + 9.0 * B**2 * C
            + 3.75 * B**2
            - 12.0 * B * C**3
            + 9.0 * B * C**2
            + 21.0 * B * C
            - 2.25 * B
            - 3.0 * C**4
            + 3.0 * C**3
            - 3.0 * C**2
        )
        / (
            2.0 * B**4
            + 8.0 * B**3 * C
            + 4.0 * B**3
            + 12.0 * B**2 * C**2
            + 6.0 * B**2 * C
            + 2.0 * B**2
            + 8.0 * B * C**3
            - 2.0 * B * C
            + 2.0 * C**4
            - 2.0 * C**3
            + 0.5 * C**2
        )
    )


print(f"MACAW value of K''_0 * K_0: {Kdprime_0K_0_MACAW(params):.2f}")

# Assign useful ranges and labels to the plots
ax[0].set_ylim(0.0, 15.0)
ax[1].set_ylim(1.0, 1.0e6)
ax[2].set_ylim(-50.0, 1000.0)
ax[3].set_ylim(1.0e3, 1.0e9)
ax[1].set_yscale("log")
ax[3].set_yscale("log")

for i, param in enumerate(["Kprime_0", "K_0", "P_0", "E_0"]):
    p = params[param]
    if i == 1 or i == 2:
        p = p / 1.0e9
    ax[i].set_xlim(0.0, 1.4)
    ax[i].plot(ax[i].get_xlim(), [p, p], c="grey", linewidth=0.5)
    ax[i].set_xlabel("$V/V_0$")

for i, param in enumerate(["Kprime_inf"]):
    p = params[param]
    if i == 1 or i == 2:
        p = p / 1.0e9
    ax[i].plot(ax[i].get_xlim(), [p, p], c="grey", linewidth=1.0)

for i in range(4):
    ax[i].set_ylim(ax[i].get_ylim())
    ax[i].plot([1.0, 1.0], ax[i].get_ylim(), c="grey", linewidth=0.5)

ax[0].set_ylabel("$K'$")
ax[1].set_ylabel("$K$ (GPa)")
ax[2].set_ylabel("Pressure (GPa)")
ax[3].set_ylabel("Internal energy (J/mol)")

bbox = dict(facecolor="white", edgecolor="none", alpha=1.0)
for i, label in enumerate(["d", "b", "a", "c"]):
    ax[i].text(
        0.04,
        0.93,
        label,
        fontsize=12,
        ha="center",
        va="center",
        transform=ax[i].transAxes,
        bbox=bbox,
    )

ax[2].legend()
fig.set_layout_engine("tight")

# Save and show the figure
fig.savefig("figures/Au_comparison.pdf")
plt.show()
