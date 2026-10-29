# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit for
# the Earth and Planetary Sciences
# Copyright (C) 2012 - 2026 by the BurnMan team, released under the GNU
# GPL v2 or later.

"""
example_hyperelastic_stress_rates
---------------------------------

In the Earth, tectonic and other geodynamic processes cause large-scale
movements of rocks. This movement can be seen as a combination of
translation, rotation, and deformation.

The most natural way to describe the deformation of a local rock element
is to track the deformation gradient F, which maps the initial configuration
of the rock element to its current configuration. This removes the
influence of rigid translation and includes both rotation and deformation.
In BurnMan, we can compute the thermodynamic and elastic response
to any applied elastic deformation gradient F. This includes the
internal energy, Cauchy stress, and elastic stiffnesses, in addition to
the more directly calculated quantities of volume change and finite
strain.

Geodynamic simulations usually process the physical response in a different
way. Instead of evaluating stress directly from the current state,
they evolve stress incrementally using the velocity gradient
and the stress retained from the previous time step. This approach fits
naturally into a calculation that includes viscous relaxation and
plastic flow. However, the retained stress must be transported with the
moving rock. In particular, rotation changes its components in fixed axes,
even if no further stretching occurs. An update based only on the rate of
stretching would miss this change.

The Truesdell, Green-Naghdi and Zaremba-Jaumann rates provide different
ways to express evolving stress in the same elastic material,
taking into account changes in configuration. The rates themselves are
exact definitions; the integration schemes used here are first-order
approximations. They differ in how they transport the stress and how they
process the constitutive response. Their integration over finite time steps
can give different answers, raising the questions: How large are the
accumulated errors, and how do they depend on the number of steps?

BurnMan's hyperelastic equation of state provides a reference for answering
these questions. It gives the exact nonlinear constitutive response of the
chosen material model, up to numerical differentiation and roundoff errors.
Stress can be evaluated directly at each prescribed deformation gradient F,
without accumulating errors through incremental updates.

In this example we run two experiments, comparing the hyperelastic
Cauchy stress with the stress derived from the three objective rates. The
first experiment rotates a prestretched crystal through 90 degrees.
Its principal stresses remain constant, while its stress components rotate.
All three objective updates reproduce the BurnMan result. Omitting rotation
leaves the stress pointing in the original direction.

The second experiment combines shear, dilation and rotation at constant
temperature. The plots compare the accumulated stresses with BurnMan's
direct calculation, then show how the errors decrease as the steps become
smaller. All three converge to the same hyperelastic response
in the limit of small time steps.

To implement the comparison, we define L = Fdot inv(F), D = sym(L),
W = skew(L), and Omega = Rdot R.T, where F = R U is the polar
decomposition of the total deformation gradient. The Cauchy-stress rates
are given below; sigma_dot follows the moving rock element::

    Zaremba-Jaumann stress rate: sigma_dot - W sigma + sigma W
    Green-Naghdi stress rate:    sigma_dot - Omega sigma + sigma Omega
    Truesdell stress rate:       sigma_dot - L sigma - sigma L.T + tr(D) sigma

BurnMan provides the current isothermal elastic stiffness B, which relates
the rate of stretching D to the Jaumann stress rate at constant temperature.
Their contraction is written B:D, with (B:D)_ij = sum_kl B_ijkl D_kl.
The full rate of change of the Cauchy stress in fixed axes is::

    sigma_dot = B:D + W sigma - sigma W

Moving the rotation terms to the left gives the Jaumann rate listed above.
Green-Naghdi subtracts rotation based on Omega instead of W. To keep the
same sigma_dot, the difference between these rotation terms must remain
on the right. Truesdell also subtracts terms involving stretching and
volume change, which must likewise be accounted for on the right.
Rearranging the same stress evolution equation therefore gives::

    Jaumann rate      = B:D
    Green-Naghdi rate = B:D + (W - Omega) sigma - sigma (W - Omega)
    Truesdell rate    = B:D - D sigma - sigma D + tr(D) sigma

References:

* Beuchert and Podladchikov (2010):
  Viscoelastic mantle convection and lithospheric stresses
  https://academic.oup.com/gji/article/183/1/35/594736
* Wallace (1967):
  Thermoelasticity of Stressed Materials and Comparison of
  Various Elastic Constants; see Eq. 2.35
  for the definition of BurnMan's stress--strain tangent B.
  https://doi.org/10.1103/PhysRev.162.776
* Objective stress rates:
  https://en.wikipedia.org/wiki/Objective_stress_rate
  https://mooseframework.inl.gov/source/materials/lagrangian/ComputeLagrangianObjectiveStress.html

*Uses:*

* :class:`burnman.HyperelasticMineral`

*Demonstrates:*

* Rotation of an initially non-hydrostatic stress.
* Conversion of the constitutive response between objective stress rates.
* Numerical convergence against a hyperelastic reference.
"""

import numpy as np
import matplotlib.pyplot as plt

if __package__:
    from .example_hyperelastic_mineral import make_mineral
else:
    from example_hyperelastic_mineral import make_mineral


RATES = ("Truesdell", "Green-Naghdi", "Zaremba-Jaumann")
SPIN_Z = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])


def rotation_z(angle):
    """
    Construct an active rotation about the z axis.

    :param angle: Rotation angle in radians, positive from x towards y.
    :return: A (3, 3) numpy array acting on column vectors.
    """
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def kinematics(t, shear, rotation, log_volume):
    """
    Evaluate shear, isotropic dilation and rotation along a prescribed path.

    All rates are derivatives with respect to the dimensionless load
    parameter t, rather than physical time.

    :param t: Load parameter in [0, 1].
    :param shear: Final dimensionless simple shear, with x displacement
        proportional to the initial y coordinate.
    :param rotation: Final superposed rotation about z, in radians.
    :param log_volume: Final logarithmic volume ratio, ln(V_final / V_initial).
    :return: Dictionary of (3, 3) numpy arrays: deformation gradient ``F``,
        stretching rate ``D``, spin ``W``, polar-rotation rate ``Omega``,
        polar rotation ``Green-Naghdi``, and integrated-spin rotation
        ``Zaremba-Jaumann``. Both rotation matrices start at the identity.
    """
    gamma = t * shear
    Q = rotation_z(t * rotation)
    F_shear = np.eye(3)
    F_shear[0, 1] = gamma
    L_shear = np.zeros((3, 3))
    L_shear[0, 1] = shear
    L = log_volume / 3.0 * np.eye(3) + rotation * SPIN_Z + Q @ L_shear @ Q.T

    # The polar angle of simple shear is -atan(gamma/2),
    # whereas integrating its spin gives -gamma/2.
    # Their difference starts at third order in gamma.
    return {
        "F": np.exp(t * log_volume / 3.0) * Q @ F_shear,
        "D": 0.5 * (L + L.T),
        "W": 0.5 * (L - L.T),
        "Omega": (rotation - 2.0 * shear / (4.0 + gamma**2)) * SPIN_Z,
        "Green-Naghdi": rotation_z(t * rotation - np.arctan(gamma / 2.0)),
        "Zaremba-Jaumann": rotation_z(t * rotation - gamma / 2.0),
    }


def sample_path(
    mineral, cell_vectors, temperature, n_steps, shear, rotation, log_volume
):
    """
    Sample the hyperelastic stress and loading (B:D) along a simple path
    in the space of deformation gradients.

    The mineral is left at the final sampled state.

    :param mineral: HyperelasticMineral with output_frame set to ``current``.
    :param cell_vectors: Initial molar cell vectors as the rows of a (3, 3)
        array, with determinant equal to the molar volume in m^3/mol.
    :param temperature: Constant temperature in K.
    :param n_steps: Positive integer number of intervals; both endpoints
        are included, giving n_steps + 1 samples.
    :param shear: Final dimensionless simple shear in the x direction
        with gradient in y.
    :param rotation: Final superposed rotation about z, in radians.
    :param log_volume: Final logarithmic volume ratio, ln(V_final / V_initial).
    :return: Dictionary with ``times`` (dimensionless load parameters,
        shape (n_steps + 1,)), ``states`` (kinematics dictionaries), and
        ``stress`` and ``loading`` arrays of shape (n_steps + 1, 3, 3).
        Stress is Cauchy stress in Pa; loading is B:D in Pa per unit load
        parameter. Both are expressed in the fixed spatial axes.
    """
    times = np.linspace(0.0, 2.0, n_steps + 1)
    states = [kinematics(t, shear, rotation, log_volume) for t in times]
    stresses, loading = [], []
    for state in states:
        mineral.set_state_with_molar_cell_vectors(
            (state["F"] @ cell_vectors.T).T, temperature
        )
        stresses.append(mineral.cauchy_stress.copy())
        if np.linalg.norm(state["D"]) == 0.0:
            loading.append(np.zeros((3, 3)))
        else:
            B = mineral.full_isothermal_stiffness_tensor
            loading.append(np.einsum("ijkl,kl->ij", B, state["D"]))
    return {
        "times": times,
        "states": states,
        "stress": np.array(stresses),
        "loading": np.array(loading),
    }


def integrate_rates(data, stride=1, consistent=True):
    """
    Integrate each stress-rate equation using a first-order update.

    Only the initial EoS stress is used to initialise the stress histories.
    Later stresses are accumulated using the rate equations, with B:D
    sampled from the EoS at the start of each update. The
    right-hand side of each stress-rate equation is evaluated
    at the start of the step.

    :param data: Dictionary returned by sample_path.
    :param stride: Positive integer number of sampled intervals per update.
        Must divide the sampled interval count; 1 uses every sample.
    :param consistent: If True, add the stress-dependent terms required
        when rewriting the Jaumann equation in Green-Naghdi or Truesdell
        form, so all three equations describe the same stress evolution.
        B:D alone is the right-hand side of the Jaumann equation.
        If False, deliberately omit these terms and set the right-hand
        side of all three equations to B:D.
    :return: Dictionary keyed by the three names in RATES and ``No rotation``.
        Each value is a Cauchy-stress history in Pa, with shape
        (n_updates + 1, 3, 3), where n_updates is the sampled interval count
        divided by stride. ``No rotation`` adds B:D without transporting stress.
    """
    n_steps = len(data["times"]) - 1
    if stride < 1 or n_steps % stride:
        raise ValueError("stride must be positive and divide the sampled step count")
    result = {name: [data["stress"][0].copy()] for name in (*RATES, "No rotation")}
    for i in range(0, n_steps, stride):
        start, end = data["states"][i], data["states"][i + stride]
        h = data["times"][i + stride] - data["times"][i]
        D, W, Omega = start["D"], start["W"], start["Omega"]
        for name, history in result.items():
            sigma = history[-1]
            source = data["loading"][i].copy()
            jacobian = 1.0
            if name == "Truesdell":
                if consistent:
                    source += -D @ sigma - sigma @ D + np.trace(D) * sigma
                transport = end["F"] @ np.linalg.inv(start["F"])
                jacobian = np.linalg.det(transport)
            elif name == "Green-Naghdi":
                if consistent:
                    spin_difference = W - Omega
                    source += spin_difference @ sigma - sigma @ spin_difference
                transport = end[name] @ start[name].T
            elif name == "Zaremba-Jaumann":
                transport = end[name] @ start[name].T
            else:
                transport = np.eye(3)
            history.append(transport @ (sigma + h * source) @ transport.T / jacobian)
    return {name: np.array(history) for name, history in result.items()}


def relative_error(predicted, reference):
    """
    Return the maximum Frobenius error divided by the largest reference norm.

    :param predicted: Predicted stress history, a numpy array of shape
        (n_samples, 3, 3).
    :param reference: Reference history at the same samples, with the same
        shape and units. Its maximum Frobenius norm must be nonzero.
    :return: Dimensionless scalar relative error, as a fraction, not a
        percentage. Each maximum is taken over all samples.
    """
    return np.max(np.linalg.norm(predicted - reference, axis=(1, 2))) / np.max(
        np.linalg.norm(reference, axis=(1, 2))
    )


def rotation_experiment(mineral, cell_vectors, temperature, n_steps):
    """
    Rotate a prestretched crystal and compare stress components in fixed axes.

    Apply a volume-preserving prestretch, then rotate through 90 degrees.
    Print the maximum stress errors and leave the mineral at the final state.

    :param mineral: HyperelasticMineral with output_frame set to ``current``.
    :param cell_vectors: Initial molar cell vectors as rows of a (3, 3) array,
        before prestretching; its determinant is the molar volume in m^3/mol.
    :param temperature: Constant temperature in K.
    :param n_steps: Positive integer number of rotation-update intervals.
    :return: A matplotlib.figure.Figure comparing the normal and shear
        stresses with the hyperelastic reference, including an update
        that omits rotation.
    """
    prestretch = np.diag(np.exp([0.04, -0.04, 0.0]))
    data = sample_path(
        mineral,
        (prestretch @ cell_vectors.T).T,
        temperature,
        n_steps,
        0.0,
        np.pi / 2.0,
        0.0,
    )
    predicted = integrate_rates(data)
    print("\nTest: Rigid rotation of a prestretched crystal (maximum relative error):")
    for name, stress in predicted.items():
        error = 100.0 * relative_error(stress, data["stress"])
        # Suppress digits below the chosen reporting precision.
        reported_error = "<0.01 %" if error < 0.01 else f"{error:.2f} %"
        print(f"  {name:20s} {reported_error}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    degrees = 90.0 * data["times"]
    for ax, component, label in zip(
        axes, [(0, 0), (0, 1)], [r"$\sigma_{xx}$", r"$\sigma_{xy}$"]
    ):
        i, j = component
        ax.plot(
            degrees,
            data["stress"][:, i, j] / 1.0e9,
            "ko",
            ms=3,
            label="Hyperelastic EoS",
        )
        # All three objective updates coincide in this rigid-rotation test.
        ax.plot(
            degrees,
            predicted["Truesdell"][:, i, j] / 1.0e9,
            label="All three objective rates",
        )
        ax.plot(
            degrees,
            predicted["No rotation"][:, i, j] / 1.0e9,
            "--",
            label="Rotation omitted",
        )
        ax.set(xlabel="Rigid rotation (degrees)", ylabel=f"{label} (GPa)")
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize="small")
    fig.suptitle(
        "A rotating crystal carries its stress with it, even when there is no deformation."
    )
    return fig


def deformation_experiment(mineral, cell_vectors, temperature, n_steps):
    """
    Conduct a convergence study of three objective stress rates.
    The study runs the same deformation path with different numbers
    of stress updates, then compares the accumulated stresses
    with the hyperelastic reference.

    The deformation path combines shear, dilation and rotation at
    constant temperature. Errors and convergence orders
    are calculated, printed to stdout and plotted.
    The mineral is left at the final state.

    :param mineral: HyperelasticMineral with output_frame set to ``current``.
    :param cell_vectors: Initial molar cell vectors as rows of a (3, 3) array,
        with determinant equal to the molar volume in m^3/mol.
    :param temperature: Constant temperature in K.
    :param n_steps: Positive integer number of updates in the stress plots.
        The convergence study uses n_steps, 2*n_steps, 4*n_steps and 8*n_steps.
    :return: A matplotlib.figure.Figure showing stress histories, polar
        rotation versus integrated spin, and stress-error convergence.
    """
    step_counts = n_steps * 2 ** np.arange(4)
    shear, rotation, log_volume = 1.0, np.pi / 2.0, 0.02
    print(
        "\nTest: Increase shear to gamma = 1, rotate 90 degrees, "
        f"and increase volume by {100*(np.exp(log_volume)-1):.2f} %."
    )
    print("Evaluating the hyperelastic stress and tangent along the path...")
    data = sample_path(
        mineral, cell_vectors, temperature, step_counts[-1], shear, rotation, log_volume
    )
    errors = {name: [] for name in RATES}
    wrong_errors = {name: [] for name in RATES[:2]}
    print("\nMaximum stress error (% of maximum EoS Frobenius norm):")
    print(f"{'Steps':>6s}" + "".join(f"{name:>20s}" for name in RATES))
    for n in step_counts:
        stride = step_counts[-1] // n
        predicted = integrate_rates(data, stride)
        wrong = integrate_rates(data, stride, consistent=False)
        reference = data["stress"][::stride]
        for name in RATES:
            errors[name].append(100.0 * relative_error(predicted[name], reference))
        for name in wrong_errors:
            wrong_errors[name].append(100.0 * relative_error(wrong[name], reference))
        print(f"{n:6d}" + "".join(f"{errors[name][-1]:20.2f}" for name in RATES))

    print("\nObserved convergence orders from the final step halving:")
    for name in RATES:
        order = np.log2(errors[name][-2] / errors[name][-1])
        print(f"  {name:20s} {order:.2f}")
    print(f"Using B:D unchanged with {step_counts[-1]} steps:")
    for name in wrong_errors:
        print(f"  {name:20s} {wrong_errors[name][-1]:.2f} % error")
    print("All consistent rates converge; the finite-step ranking is not universal.")

    stride = step_counts[-1] // n_steps
    predicted = integrate_rates(data, stride)
    times = data["times"]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    for ax, component, label in zip(
        axes[0], [(0, 1), (0, 0)], [r"$\sigma_{xy}$", r"$\sigma_{xx}$"]
    ):
        i, j = component
        ax.plot(
            times, data["stress"][:, i, j] / 1.0e9, "k", lw=2, label="Hyperelastic EoS"
        )
        for name in RATES:
            ax.plot(times[::stride], predicted[name][:, i, j] / 1.0e9, ".-", label=name)
        ax.plot(
            times[::stride],
            predicted["No rotation"][:, i, j] / 1.0e9,
            "--",
            color="0.5",
            label="Rotation omitted",
        )
        ax.set(
            xlabel="Load parameter t",
            ylabel=f"{label} (GPa)",
            title=f"{n_steps} stress-update steps",
        )
        ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize="small")

    ax = axes[1, 0]
    gamma = times * shear
    # Q is the imposed rotation; the polar rotation also includes shear's rotation.
    ax.plot(
        times,
        np.degrees(times * rotation),
        "k:",
        label="Imposed rotation (Q)",
    )
    ax.plot(
        times,
        np.degrees(times * rotation - np.arctan(gamma / 2.0)),
        label="Imposed + simple shear rotation (= Green-Naghdi)",
    )
    ax.plot(
        times,
        np.degrees(times * rotation - gamma / 2.0),
        "--",
        label="Integrated spin (= Jaumann)",
    )
    ax.set(
        xlabel="Load parameter t",
        ylabel="Rotation angle (degrees)",
        title="Polar rotation and material spin differ in shear",
    )
    ax.legend(fontsize="small")
    ax.grid(alpha=0.25)

    ax = axes[1, 1]
    for i, name in enumerate(RATES):
        ax.loglog(step_counts, errors[name], "o-", color=f"C{i}", label=name)
        if name in wrong_errors:
            ax.loglog(
                step_counts,
                wrong_errors[name],
                "--",
                color=f"C{i}",
                label=f"{name}, unchanged B:D",
            )
    ax.set(
        xlabel="Number of stress-update steps",
        ylabel="Maximum relative stress error (%)",
        title="Correct rate conversion is essential",
    )
    ax.legend(fontsize="small")
    ax.grid(which="both", alpha=0.25)
    fig.suptitle(
        "Shear, dilation and rotation: one hyperelastic material, three objective rates"
    )
    return fig


if __name__ == "__main__":
    # Steps in the rotation and stress-response plots; convergence uses N to 8N.
    n_steps = 8

    mineral = make_mineral()
    temperature = 1000.0
    pressure = 5.0e9
    mineral.hydrostatic.set_state(pressure, temperature)
    cell_vectors = mineral.hydrostatic.unrotated_cell_vectors.copy()
    mineral.set_output_frame("current")
    print(
        f"Objective Cauchy-stress rates at {temperature:.1f} K; "
        f"initial pressure {pressure/1.0e9:.1f} GPa."
    )
    rotation_experiment(mineral, cell_vectors, temperature, n_steps)
    deformation_experiment(mineral, cell_vectors, temperature, n_steps)
    plt.show()
