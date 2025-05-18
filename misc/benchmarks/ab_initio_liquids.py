import burnman
from burnman.minerals import DKS_2013_liquids
from burnman import constants
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from burnman.tools.eos import check_eos_consistency
from burnman.tools.chemistry import hugoniot


SiO2_liq = DKS_2013_liquids.SiO2_liquid()
MgO_liq = DKS_2013_liquids.MgO_liquid()

MgO_temperatures = np.array([2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0, 10000.0])
MgO_volumes = np.linspace(7e-6, 18e-6, 101)

SiO2_temperatures = np.array([2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0])
SiO2_volumes = np.linspace(9e-6, 30e-6, 101)

plots = [
    [
        ["figures/SiO2_liquid_PVT.png", [9, 30, -10, 220]],
        ["figures/SiO2_liquid_SelVT.png", [9, 30, -0.03, 0.75]],
        ["figures/SiO2_liquid_EVT.png", [9, 30, -2400, -1200]],
    ],
    [
        ["figures/MgO_liquid_PVT.png", [7, 18, -6, 240]],
        ["figures/MgO_liquid_SelVT.png", [6, 18, -0.04, 0.84]],
        ["figures/MgO_liquid_EVT.png", [7, 18, -1200, -200]],
    ],
]


for i, (phase, n_atoms, temperatures, volumes) in enumerate(
    [
        (SiO2_liq, 3.0, SiO2_temperatures, SiO2_volumes),
        (MgO_liq, 2.0, MgO_temperatures, MgO_volumes),
    ]
):
    print(
        "EoS consistent for {0} model: {1}".format(
            phase.name,
            check_eos_consistency(phase, tol=1.0e-4, including_shear_properties=False),
        )
    )
    fig = plt.figure()
    ax_P = fig.add_subplot(1, 3, 1)
    ax_S = fig.add_subplot(1, 3, 2)
    ax_E = fig.add_subplot(1, 3, 3)
    ax_P.set_xlabel("Volume (cm^3/mol)")
    ax_S.set_xlabel("Volume (cm^3/mol)")
    ax_E.set_xlabel("Volume (cm^3/mol)")

    ax_P.set_ylabel("Pressure (GPa)")
    ax_S.set_ylabel("Entropy/nR")
    ax_E.set_ylabel("Internal Energy (kJ/mol)")

    ax_P.imshow(mpimg.imread(plots[i][0][0]), extent=plots[i][0][1], aspect="auto")
    ax_S.imshow(mpimg.imread(plots[i][1][0]), extent=plots[i][1][1], aspect="auto")
    ax_E.imshow(mpimg.imread(plots[i][2][0]), extent=plots[i][2][1], aspect="auto")

    for temperature in temperatures:
        pressures = np.empty_like(volumes) + np.nan
        entropies = np.empty_like(volumes) + np.nan
        energies = np.empty_like(volumes) + np.nan

        for j, volume in enumerate(volumes):
            try:
                phase.set_state_with_volume(volume, temperature)
                pressures[j] = phase.pressure
                entropies[j] = phase.method._S_el(temperature, volume, phase.params)
                energies[j] = phase.method._molar_internal_energy(
                    phase.pressure, temperature, volume, phase.params
                )
                energies[j] = phase.molar_internal_energy
            except Exception as e:
                print(e.message)

        ax_P.plot(
            volumes * 1.0e6,
            pressures / 1.0e9,
            linewidth=2,
            label="{0:.0f} K".format(temperature),
        )
        ax_S.plot(
            volumes * 1.0e6,
            entropies / n_atoms / constants.gas_constant,
            linewidth=2,
            label="{0:.0f} K".format(temperature),
        )
        ax_E.plot(
            volumes * 1.0e6,
            energies / 1.0e3,
            linewidth=2,
            label="{0:.0f} K".format(temperature),
        )

        ax_E.legend(loc="upper right")

    fig.set_layout_engine("tight")
    plt.show()


# Fe2SiO4
# Plot EoS and hugoniot

fa_liq = burnman.minerals.RS_2014_liquids.Fe2SiO4_liquid()

print(
    "EoS consistent for {0} model: {1}".format(
        fa_liq.name,
        check_eos_consistency(fa_liq, tol=1.0e-4, including_shear_properties=False),
    )
)
fig = plt.figure()
ax_P = fig.add_subplot(1, 2, 1)
ax_hugoniot = fig.add_subplot(1, 2, 2)


ax_P.imshow(
    mpimg.imread("figures/Fe2SiO4_liquid_PVT.png"),
    extent=[3.5, 7.5, 0.0, 200],
    aspect="auto",
)
ax_hugoniot.imshow(
    mpimg.imread("figures/Fe2SiO4_liquid_hugoniot.png"),
    extent=[3.5, 7.5, 0, 200],
    aspect="auto",
)

pressures = np.linspace(1.0e5, 200.0e9, 51)
for T in [3000.0, 4000.0, 6000.0]:
    rhos = fa_liq.evaluate(["rho"], pressures, [T] * len(pressures))[0]
    ax_P.plot(rhos / 1.0e3, pressures / 1.0e9, label="{0:.0f} K".format(T))

ax_P.set_xlabel("Densities (Mg/m^3)")
ax_P.set_ylabel("Pressures (GPa)")
ax_hugoniot.set_title("PVT")
ax_P.legend(loc="upper left")


temperatures, volumes = hugoniot(fa_liq, 1.0e5, 1573.0, pressures)
ax_hugoniot.plot(fa_liq.molar_mass / volumes / 1.0e3, pressures / 1.0e9)

ax_hugoniot.set_xlabel("Densities (Mg/m^3)")
ax_hugoniot.set_ylabel("Pressures (GPa)")
ax_hugoniot.set_title("hugoniot")

fig.set_layout_engine("tight")
plt.show()
