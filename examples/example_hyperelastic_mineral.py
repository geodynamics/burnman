# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_hyperelastic_mineral
---------------------------

This example asks the question:
Does a crystal change thickness when sheared under a constant normal load?

In the purely incompressible case, the answer is no:
the plate separation is fixed. It is easy to see why this is so;
the area of a parallelogram is the product of its base and height.
If the area and base are fixed, the height must be fixed as well.

In a compressible model, the answer is not as straightforward,
because the volume is free to change. The change in thickness is
a kind of Poynting effect, named after John Poynting, who described
it in 1909 (https://doi.org/10.1098/rspa.1909.0059).
The sign and magnitude of the change depends on the constitutive law.

In this example, we use a hyperelastic model of an orthorhombic mineral
and shear it between parallel plates at fixed temperature. We then solve
for the plate separation that maintains the normal traction. The example
compares sliding along a and c, with plate normals along b.

In this model, shearing makes the crystal thicken. Since the area
of the base of the parallelogram is fixed, its volume increases by the
same fraction as the plate separation. The hydrostatic EoS pressure
at that larger volume and the same temperature falls, as expected.
However, the pressure (minus one third of the trace of the stress tensor)
of the deformed crystal actually rises. This is because
the imposed shear induces compressive normal stresses. There is no conflict
with the requirement that the material must have positive compressibility,
because the loading path changes the unit cell shape as well as its volume.

*Uses:*

* :class:`burnman.AnisotropicMineral`
* :class:`burnman.HyperelasticMineral`

*Demonstrates:*

* Solving a mixed displacement/traction boundary condition at finite strain.
* Orientation-dependent shear response and normal strain.
* Checking mechanical work against a thermodynamic potential.
* Plotting the loading paths and their departure from linear elasticity.
* Visualizing the cell shape at successive shear displacements.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.integrate import simpson
from scipy.optimize import brentq

from burnman import AnisotropicMineral, HyperelasticMineral
from burnman.minerals import SLB_2011


def make_mineral():
    """
    Make a hyperelastic model based on the SLB_2011 forsterite EoS.
    """
    fo = SLB_2011.forsterite()
    lengths = np.array([4.7646, 10.2296, 5.9942])
    lengths *= np.cbrt(fo.params["V_0"] / np.prod(lengths))
    cell_parameters = np.concatenate((lengths, [90.0, 90.0, 90.0]))

    # Leading compliance coefficients in the expansion of Psi in ln(V/V0).
    constants = np.zeros((6, 6, 2, 1))
    constants[:, :, 1, 0] = np.array(
        [
            [0.44, -0.12, -0.10, 0.0, 0.0, 0.0],
            [-0.12, 0.78, -0.22, 0.0, 0.0, 0.0],
            [-0.10, -0.22, 0.66, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.97, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.61, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.55],
        ]
    )
    return HyperelasticMineral(AnisotropicMineral(fo, cell_parameters, constants))


def analytical_normal_strain_coefficient(hydrostatic_mineral, i, j):
    """
    Predict the normal strain (lambda) at a given shear strain (gamma).
    The governing equation can be derived by minimizing the appropriate
    potential (A + sigma_n * V), where sigma_n is the fixed normal stress
    and V is the volume.

    :return: The coefficient of the normal_strain_coefficient
        (lambda - 1) / gamma**2 for shear along orthorhombic crystal axes.
    :rtype: float
    """
    pressure, temperature = (
        hydrostatic_mineral.pressure,
        hydrostatic_mineral.temperature,
    )
    C = hydrostatic_mineral.full_isothermal_stiffness_tensor.copy()
    K = hydrostatic_mineral.isothermal_bulk_modulus_reuss
    dP = 1.0e6  # Pa; a small hydrostatic pressure increment
    try:
        hydrostatic_mineral.set_state(pressure + dP, temperature)
        G_plus = hydrostatic_mineral.full_isothermal_stiffness_tensor[i, j, i, j]
        hydrostatic_mineral.set_state(pressure - dP, temperature)
        G_minus = hydrostatic_mineral.full_isothermal_stiffness_tensor[i, j, i, j]
    finally:
        hydrostatic_mineral.set_state(pressure, temperature)
    Gprime = (G_plus - G_minus) / (2.0 * dP)
    return (2.0 * C[i, j, i, j] + 2.0 * K * Gprime - C[j, j, j, j] + C[j, j, i, i]) / (
        4.0 * C[j, j, j, j]
    )


def equilibrate_gap(m, cell_vectors_initial, temperature, normal_stress, gamma, i, j):
    """
    Solve for the normal stretch at prescribed shear displacement and load.
    """

    def residual(log_stretch):
        F = np.eye(3)
        F[i, j] = gamma
        F[j, j] = np.exp(log_stretch)
        # Cell vectors are rows; deformation gradients act on column vectors.
        m.set_state_with_molar_cell_vectors((F @ cell_vectors_initial.T).T, temperature)
        return (m.cauchy_stress[j, j] + normal_stress) / 1.0e9

    # The modest shears below have equilibria within 10% of the initial gap.
    root = brentq(residual, np.log(0.9), np.log(1.1), xtol=1.0e-11)
    residual(root)  # Leave the mineral at the converged state.
    return np.exp(root)


def plot_cell_shapes(cell_paths, gammas):
    """
    Show the cell parallelepipeds in fixed axes without magnifying the strain.
    """
    corners = np.array(list(np.ndindex(2, 2, 2)))
    edges = np.array(
        [
            (i, j)
            for i in range(8)
            for j in range(i + 1, 8)
            if np.count_nonzero(corners[i] - corners[j]) == 1
        ]
    )
    # Normalize all states by the initial b length,
    # preserving cell proportions.
    paths = [
        cell_vectors / np.linalg.norm(cell_vectors[0, 1])
        for cell_vectors in cell_paths.values()
    ]
    vertices = np.concatenate(
        [corners @ cell_vectors for path in paths for cell_vectors in path]
    )
    lower = vertices.min(axis=0) - 0.08
    upper = vertices.max(axis=0) + 0.08
    steps = [0, len(gammas) // 2, len(gammas) - 1]
    fig, axes = plt.subplots(
        2, 3, figsize=(12, 8), subplot_kw={"projection": "3d"}, layout="constrained"
    )
    fig.suptitle(
        "Cell shape during shear (plate normal along b)\n"
        "Dashed grey: initial cell; lengths scaled by initial b; no strain magnification"
    )
    for row, (label, path) in enumerate(zip(cell_paths, paths)):
        initial_vertices = corners @ path[0]
        for ax, step in zip(axes[row], steps):
            # Cell vectors are rows; corners contain fractional coordinates.
            current_vertices = corners @ path[step]
            ax.add_collection3d(
                Line3DCollection(
                    initial_vertices[edges], colors="0.6", linestyles="--", linewidths=1
                )
            )
            ax.add_collection3d(
                Line3DCollection(
                    current_vertices[edges], colors=f"C{row}", linewidths=2
                )
            )
            ax.set(
                xlim=(lower[0], upper[0]),
                ylim=(lower[1], upper[1]),
                zlim=(lower[2], upper[2]),
                xlabel=r"$x/b_0$",
                ylabel=r"$y/b_0$",
                zlabel=r"$z/b_0$",
                title=f"Slide along {label}: shear {100.0 * gammas[step]:g}%",
            )
            ax.set_box_aspect(upper - lower, zoom=0.85)
            ax.view_init(elev=20, azim=-55, vertical_axis="y")
    return fig


if __name__ == "__main__":
    pressure = 5.0e9  # Pa; positive for compression
    temperature = 1000.0  # K
    gammas = np.linspace(0.0, 0.10, 11)

    print("Shearing mineral between plates")
    print(f"Normal load: {pressure / 1.0e9:.1f} GPa; temperature: {temperature:.0f} K")
    print("In-plane stretches are fixed; the plate separation can relax.")
    print("Positive gap change means thickening.\n")

    fig_response, ax_response = plt.subplots(
        1, 3, figsize=(15, 4), layout="constrained"
    )
    fig_energy, ax_energy = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    conditions = f"{pressure / 1.0e9:g} GPa normal load, {temperature:g} K"
    fig_response.suptitle(f"Shear and thickening of mineral: {conditions}")
    fig_energy.suptitle("Work balance and normal stress redistribution")
    cell_paths = {}

    for i, j, label in [(0, 1, "a"), (2, 1, "c")]:
        m = make_mineral()
        m.hydrostatic.set_state(pressure, temperature)
        cell_vectors_initial = m.hydrostatic.cell_vectors.copy()
        volume_initial = m.hydrostatic.V
        shear_modulus = m.hydrostatic.full_isothermal_stiffness_tensor[i, j, i, j]
        expansion_coeff = analytical_normal_strain_coefficient(m.hydrostatic, i, j)
        m.set_state_with_molar_cell_vectors(cell_vectors_initial, temperature)
        potential_initial = m.helmholtz + pressure * m.V

        shear_stresses = []
        stretches = []
        relative_volumes = []
        potentials = []
        mean_pressures = []
        hydrostatic_pressures = []
        cell_vectors = []
        print(f"Sliding along {label}, plate normal along b")
        print(f"Initial shear modulus: {shear_modulus / 1.0e9:.3f} GPa")
        print(f"Prediction: lambda - 1 = {expansion_coeff:.3f} gamma^2")
        print(
            " gamma   gap change (%)   shear (GPa)        P (GPa)"
            "     P_hyd (GPa)   delta Phi (J/mol)"
        )
        for gamma in gammas:
            stretch = equilibrate_gap(
                m, cell_vectors_initial, temperature, pressure, gamma, i, j
            )
            stress = m.cauchy_stress.copy()
            delta_potential = m.helmholtz + pressure * m.V - potential_initial
            shear_stresses.append(stress[i, j])
            stretches.append(stretch)
            relative_volumes.append(m.V / volume_initial - 1.0)
            potentials.append(delta_potential)
            mean_pressures.append(-np.trace(stress) / 3.0)
            hydrostatic_pressures.append(m.hydrostatic.pressure)
            residual = np.abs(stress[j, j] + pressure)
            assert residual < 1.0, f"Residual normal stress {residual:.2f} Pa"
            cell_vectors.append(m.cell_vectors.copy())
            print(
                f"{gamma:6.2f}"
                f"   {100.0 * (stretch - 1.0):14.2f}"
                f"   {np.abs(stress[i, j]) / 1.0e9:11.2f}"
                f"   {mean_pressures[-1] / 1.0e9:12.2f}"
                f"   {hydrostatic_pressures[-1] / 1.0e9:13.2f}"
                f"   {delta_potential:17.0f}"
            )

        # The current plate area is unchanged, so V_initial * tau is the
        # generalized force conjugate to gamma (not V_current * tau).
        work = volume_initial * simpson(shear_stresses, x=gammas)
        relative_error = np.abs(work - delta_potential) / np.abs(delta_potential)
        print(f"Integrated shear work: {work:.0f} J/mol")
        print(f"Final delta Phi:       {delta_potential:.0f} J/mol")
        print(f"Relative work/energy discrepancy: {relative_error:.2f} J/mol")
        print(
            f"Final relative volume change: {100.0 * np.abs(relative_volumes[-1]):.3f} %"
        )
        print(
            f"Leading-order prediction: {100.0 * expansion_coeff * gammas[-1]**2:.3f} %"
        )
        assert relative_error < 1.0e-5
        cell_paths[label] = np.array(cell_vectors)

        color = "C0" if label == "a" else "C1"
        shear_percent = 100.0 * gammas
        ax_response[0].plot(
            shear_percent,
            np.array(shear_stresses) / 1.0e9,
            color=color,
            label=f"Slide along {label}",
        )
        ax_response[0].plot(
            shear_percent,
            shear_modulus * gammas / 1.0e9,
            "--",
            color=color,
            label=f"{label}: linear prediction",
        )
        ax_response[1].plot(
            shear_percent,
            100.0 * (np.array(stretches) - 1.0),
            color=color,
            label=f"Slide along {label}",
        )
        ax_response[2].plot(
            shear_percent,
            100.0 * np.array(relative_volumes),
            color=color,
            label=f"Slide along {label}",
        )
        for ax in ax_response[1:]:
            ax.plot(
                shear_percent,
                100.0 * expansion_coeff * gammas**2,
                "--",
                color=color,
                label=f"{label}: quadratic prediction",
            )
        ax_energy[0].plot(
            shear_percent,
            np.array(potentials) / 1.0e3,
            color=color,
            label=rf"{label}: $\Delta\Phi$",
        )
        # Use an even number of intervals for each composite Simpson integral.
        # The markers should lie on the independently calculated energy curves.
        work_indices = np.arange(0, len(gammas), 2)
        cumulative_work = [0.0] + [
            volume_initial * simpson(shear_stresses[: k + 1], x=gammas[: k + 1])
            for k in work_indices[1:]
        ]
        ax_energy[0].plot(
            shear_percent[work_indices],
            np.array(cumulative_work) / 1.0e3,
            "o",
            color=color,
            markerfacecolor="none",
            label=f"{label}: shear work",
        )
        ax_energy[1].plot(
            shear_percent,
            np.array(mean_pressures) / 1.0e9,
            color=color,
            label=f"{label}: total mean pressure",
        )
        ax_energy[1].plot(
            shear_percent,
            np.array(hydrostatic_pressures) / 1.0e9,
            "--",
            color=color,
            label=f"{label}: hydrostatic at same V, T",
        )

    print("The two shear directions give different gap changes and shear stresses.")

    ax_response[0].set_ylabel("Shear traction (GPa)")
    ax_response[1].set_ylabel("Plate gap change (%)")
    ax_response[2].set_ylabel(
        r"Relative molar volume change $\Delta V/V_{\rm initial}$ (%)"
    )
    ax_energy[0].set_ylabel("Energy change (kJ/mol)")
    ax_energy[1].set_ylabel("Pressure (GPa)")
    ax_energy[1].axhline(
        pressure / 1.0e9, color="0.4", linestyle=":", label="Fixed normal load"
    )
    for ax in [*ax_response, *ax_energy]:
        ax.set_xlabel(r"Imposed shear $100\gamma$ (%)")
        ax.set_xlim(0.0, shear_percent[-1])
        ax.grid(alpha=0.2)
        ax.legend()

    fig_cell_vectors = plot_cell_shapes(cell_paths, gammas)

    plt.show()
