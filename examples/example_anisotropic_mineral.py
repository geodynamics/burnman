# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2021 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""
example_anisotropic_mineral
---------------------------

This example illustrates how to create and interrogate an AnisotropicMineral
object.

*Specifically uses:*

* :class:`burnman.AnisotropicMineral`

*Demonstrates:*

* anisotropic functions

"""
from __future__ import absolute_import
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from burnman import AnisotropicMineral
from burnman.minerals import SLB_2011
from burnman.tools.eos import check_anisotropic_eos_consistency
from burnman.tools.plot import plot_projected_elastic_properties


if __name__ == "__main__":
    # Let's create a first approximation to an olivine mineral.
    # Olivine is orthorhombic.
    # BurnMan has an AnisotropicMineral class, which requires as input:
    # 1) An isotropic mineral (from which it gets the V-P-T relations)
    # 2) A standard state cell parameters vector
    # 3) A 4D numpy array containing constants describing the
    #    anisotropic tensor state function.

    fo = SLB_2011.forsterite()
    cell_lengths = np.array([4.7646, 10.2296, 5.9942])
    cell_lengths *= np.cbrt(fo.params["V_0"] / np.prod(cell_lengths))

    cell_parameters = np.array(
        [cell_lengths[0], cell_lengths[1], cell_lengths[2], 90, 90, 90]
    )

    # The constants function is given as an expansion in ln(V/V0) and
    # thermal pressure (see paper). Here we are only interested in a single
    # 6x6 block with indexing constants[:, :, 1, 0], which corresponds to
    # S_{Tpq} / beta_RT, where S_{Tpq} is the isothermal compliance tensor
    # in Voigt form (values in the off-diagonal blocks and lower diagonal block
    # are multiplied by factors of 2 and 4 respectively, and beta_RT is
    # the Reuss isothermal compressibility,
    # both measured at the reference state.
    # This block is the most important; the other blocks are pressure and
    # temperature corrections.
    constants = np.zeros((6, 6, 2, 1))
    constants[:, :, 1, 0] = np.array(
        [
            [0.44, -0.12, -0.1, 0.0, 0.0, 0.0],
            [-0.12, 0.78, -0.22, 0.0, 0.0, 0.0],
            [-0.1, -0.22, 0.66, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.97, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.61, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.55],
        ]
    )

    m = AnisotropicMineral(fo, cell_parameters, constants)

    # The following line checks that the mineral we just created is
    # internally consistent.
    assert check_anisotropic_eos_consistency(m)

    # Now we can set state and interrogate the
    # mineral for various anisotropic properties.
    # Here is a choice selection.

    P = 1.0e9
    T = 1600.0
    m.set_state(P, T)
    print(f"Model forsterite properties at {P/1.e9:.2f} GPa and {T:.2f} K:")
    np.set_printoptions(formatter={"all": lambda x: f"{x:0.4f}"})
    print("Cell vectors:")
    print(m.cell_vectors)
    print("Cell parameters:")
    print(m.cell_parameters)
    print("Thermal expansivity (1/MK):")
    print(m.thermal_expansivity_tensor * 1.0e6)
    print("Isothermal compressibility (1/GPa):")
    print(m.isothermal_compressibility_tensor * 1.0e9)
    print("Thermal stress tensor (MPa/K):")
    print(m.thermal_stress_tensor / 1.0e6)
    print("Isothermal stiffness_tensor (GPa):")
    print(m.isothermal_stiffness_tensor / 1.0e9)
    print("Isentropic stiffness_tensor (GPa):")
    print(m.isentropic_stiffness_tensor / 1.0e9)
    print("Grueneisen tensor:")
    print(m.grueneisen_tensor)

    # We can also obtain the scalar properties that are inherited from
    # the isotropic equation of state
    print(f"Volume: {m.V * 1.e6:.4f} cm^3/mol")
    print(f"Entropy: {m.S:.2f} J/K/mol")
    print(f"C_P: {m.molar_heat_capacity_p:.2f} J/K/mol")
    print(f"C_V: {m.molar_heat_capacity_v:.2f} J/K/mol")
    print(f"C_eps: {m.molar_isometric_heat_capacity:.2f} J/K/mol")

    # Plot thermal expansion figure
    fig = plt.figure(figsize=(8, 4))
    ax = [fig.add_subplot(1, 2, i) for i in range(1, 3)]

    temperatures = np.linspace(10.0, 1600.0, 101)
    alphas = np.empty((101, 4))
    extensions = np.empty((101, 3))
    vectors = np.empty((101, 4))

    labels = ["a", "b", "c", "V"]

    for i, T in enumerate(temperatures):
        m.set_state(1.0e5, T)
        alphas[i, :3] = np.diag(m.thermal_expansivity_tensor) * 1.0e5
        alphas[i, 3] = m.alpha * 1.0e5 / 3.0
        extensions[i] = (
            (np.diag(m.cell_vectors) / np.diag(m.cell_vectors_0)) - 1.0
        ) * 1.0e4
        vectors[i, :3] = np.diag(m.cell_vectors)

        vectors[i, 3] = m.V

    for i in range(4):
        label = f"$\\alpha_{{{labels[i]}}}$"
        if i == 3:
            ln = ax[0].plot(temperatures, alphas[:, i], label=label + "/3")
        else:
            ax[0].plot(temperatures, alphas[:, i], label=label)

    for i in range(3):
        ax[1].plot(temperatures, extensions[:, i], label=labels[i])

    Vthird_expansion = 1.0e4 * (
        np.power(np.prod(extensions * 1.0e-4 + 1, axis=1), 1.0 / 3.0) - 1.0
    )
    ln = ax[1].plot(temperatures, Vthird_expansion, label="$V^{1/3}$")

    ax[0].set_ylim(
        0.0,
    )

    for i in range(2):
        ax[i].set_xlim(0.0, 1600.0)
        ax[i].set_xlabel("Temperature (K)")
        ax[i].legend()

    ax[0].set_ylabel("Thermal expansivity (10$^{-5}$/K)")
    ax[1].set_ylabel("Relative length change ($10^{4} (x/x_0 - 1)$)")

    fig.set_tight_layout(True)
    # fig.savefig('example_anisotropic_mineral_Figure_1.png')
    plt.show()

    # Start plotting Cij figure
    fig = plt.figure(figsize=(12, 12))
    ax = [fig.add_subplot(3, 3, i) for i in range(1, 10)]

    pressures = np.linspace(1.0e7, 30.0e9, 101)
    G_iso = np.empty_like(pressures)
    G_aniso = np.empty_like(pressures)
    C = np.empty((len(pressures), 6, 6))

    f = np.empty_like(pressures)
    dXdf = np.empty_like(pressures)

    i_pq = ((1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6), (1, 2), (1, 3), (2, 3))

    temperatures = [500.0, 1000.0, 1500.0, 2000.0]
    for T in temperatures:
        for i, P in enumerate(pressures):
            m.set_state(P, T)
            C[i] = m.isentropic_stiffness_tensor

        for i, (p, q) in enumerate(i_pq):
            ln = ax[i].plot(
                pressures / 1.0e9, C[:, p - 1, q - 1] / 1.0e9, label=f"{T} K"
            )

    for i, (p, q) in enumerate(i_pq):
        ax[i].set_xlabel("Pressure (GPa)")
        ax[i].set_ylabel(f"$C_{{N {p}{q}}}$ (GPa)")
        ax[i].legend()

    fig.set_tight_layout(True)
    # fig.savefig('example_anisotropic_mineral_Figure_2.png')
    plt.show()

    # Finally, we make a pretty plot of various elastic/seismic properties
    # at a fixed pressure and temperature.
    fig = plt.figure(figsize=(12, 10.5))
    ax = [fig.add_subplot(3, 3, i, projection="polar") for i in range(1, 10)]

    P = 3.0e9
    T = 1600.0
    m.set_state(P, T)
    plot_types = [
        "vp",
        "vs1",
        "vp/vs1",
        "vp/vs2",
        "s anisotropy",
        "linear compressibility",
        "youngs modulus",
        "minimum poisson ratio",
        "maximum poisson ratio",
    ]

    contour_sets, ticks, lines = plot_projected_elastic_properties(m, plot_types, ax)
    for i in range(len(contour_sets)):
        cbar = fig.colorbar(contour_sets[i], ax=ax[i], ticks=ticks[i], pad=0.1)
        cbar.add_lines(lines[i])

    fig.set_tight_layout(True)
    plt.show()
