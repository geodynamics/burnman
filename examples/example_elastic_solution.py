import burnman
import numpy as np
import matplotlib.pyplot as plt
from burnman import minerals

"""

example_elastic_solution
--------------------------

This example shows how to create an elastic solution model
and output thermodynamic and thermoelastic quantities of
that model.

There are four main types of elastic solution currently implemented in
BurnMan:

1. Ideal solutions
2. Symmmetric solutions
3. Asymmetric solutions
4. Subregular solutions

These mirror the types of ordinary solution models
(those defined as a function of pressure and temperature,
rather than volume and temperature)

These solutions can potentially deal with:

* Disordered endmembers (more than one element on a crystallographic site)
* Site vacancies
* More than one valence/spin state of the same element on a site

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.ElasticSolution`
* :class:`burnman.ElasticSolutionModel`


*Demonstrates:*

* Different ways to define an elastic solution
* How to set composition and state
* How to output thermodynamic and thermoelastic properties

"""

if __name__ == "__main__":

    """
    First, let's create an instance of an elastic solution.
    This is done in a very similar way to the definition of
    normal solutions, except that we use the ElasticSolution
    class, and rather than volume interactions, we define pressure
    interactions.

    The following line creates an instance of a pyrope-grossular garnet.
    The model is chosen to have symmetric interaction parameters,
    with a constant interaction energy of -60 kJ/mol. This
    interaction energy corresponds to partial relaxation of the energy
    of interaction between pyrope-like and grossular-like structural
    groups in the solution structure.
    """
    gt_F = burnman.ElasticSolution(
        name="Garnet",
        solution_type="symmetric",
        endmembers=[
            [minerals.SLB_2011.pyrope(), "[Mg]3[Al]2Si3O12"],
            [minerals.SLB_2011.grossular(), "[Ca]3[Al]2Si3O12"],
        ],
        energy_interaction=[[-60.0e3]],
        entropy_interaction=[[0.0]],
        pressure_interaction=[[0.0e9]],
    )

    """
    For comparison, let's define a normal (P-T) subregular solution
    for pyrope grossular, matching that proposed in the study of
    Ganguly et al. (1996).
    """

    def mult(x, n):
        return [[[v * n for v in i] for i in j] for j in x]

    gt_G = burnman.Solution(
        name="Subregular pyrope-almandine-grossular " "garnet (Ganguly et al., 1996)",
        solution_type="subregular",
        endmembers=[
            [minerals.SLB_2011.pyrope(), "[Mg]3[Al]2Si3O12"],
            [minerals.SLB_2011.grossular(), "[Ca]3[Al]2Si3O12"],
        ],
        energy_interaction=mult([[[9834.0, 21627.0]]], 3.0),
        volume_interaction=mult([[[0.058e-5, 0.012e-5]]], 3.0),
        entropy_interaction=mult([[[5.78, 5.78]]], 3.0),
    )

    """
    Let's look at the properties of the two solutions at 1 GPa and 500 C
    (773.15 K).
    """
    P = 1.0e9
    T = 500.0 + 273.15

    xs = np.linspace(0.0, 1.0, 101)
    dGs_F = np.empty_like(xs)
    dGs_G = np.empty_like(xs)
    dVs_F = np.empty_like(xs)
    dVs_G = np.empty_like(xs)
    dSs_F = np.empty_like(xs)
    dSs_G = np.empty_like(xs)

    gt_G.set_state(P, T)
    gt_F.set_state(P, T)

    for i, x in enumerate(xs):
        gt_G.set_composition([x, 1.0 - x])
        gt_F.set_composition([x, 1.0 - x])

        dGs_F[i] = gt_F.molar_gibbs
        dGs_G[i] = gt_G.molar_gibbs
        dVs_F[i] = gt_F.molar_volume
        dVs_G[i] = gt_G.molar_volume
        dSs_F[i] = gt_F.molar_entropy
        dSs_G[i] = gt_G.molar_entropy

    x_py = 0.4
    proportions = np.array([x_py, 1.0 - x_py])
    gt_F.set_composition(proportions)
    gt_G.set_composition(proportions)

    """
    Print the activity coefficients at p(py) = 0.4
    """
    print(f"py, gr activities at {P/1.e9:.1f} GPa, {T:.0f} K and x_py={x_py:.2f}:")
    print(f"[{gt_F.activities[0]:.2f}, {gt_F.activities[1]:.2f}]")

    dGs_mech = (1.0 - xs) * dGs_F[0] + xs * dGs_F[-1]
    dVs_mech = (1.0 - xs) * dVs_F[0] + xs * dVs_F[-1]
    dSs_mech = (1.0 - xs) * dSs_F[0] + xs * dSs_F[-1]

    """
    Plot the properties
    """
    fig = plt.figure(figsize=(12, 4))
    ax = [fig.add_subplot(1, 3, i) for i in range(1, 4)]

    ax[0].plot(xs, (dGs_F - dGs_mech) / 1.0e3)
    ax[0].plot(xs, (dGs_G - dGs_mech) / 1.0e3)
    ax[1].plot(xs, (dVs_F - dVs_mech) * 1.0e6)
    ax[1].plot(xs, (dVs_G - dVs_mech) * 1.0e6)
    ax[2].plot(xs, dSs_F - dSs_mech, label="Elastic solution")
    ax[2].plot(xs, dSs_G - dSs_mech, label="Ganguly et al. (1996)")

    x_py = 0.4
    gt_F.set_composition([x_py, 1.0 - x_py])
    ax[0].plot([1.0, 0.0], (gt_F.partial_gibbs - [dGs_F[-1], dGs_F[0]]) / 1.0e3)
    ax[1].plot([1.0, 0.0], (gt_F.partial_volumes - [dVs_F[-1], dVs_F[0]]) * 1.0e6)
    ax[2].plot(
        [1.0, 0.0],
        gt_F.partial_entropies - [dSs_F[-1], dSs_F[0]],
        label="Partials at x(py)=0.4",
    )

    ax[2].legend()

    for i in range(3):
        ax[i].set_xlim(0.0, 1.0)
        ax[i].set_xlabel("x(py)")

    ax[0].set_ylabel("Gibbs excess (kJ/mol)")
    ax[1].set_ylabel("Volume excess (cm$^3$/mol)")
    ax[2].set_ylabel("Entropy excess (J/K/mol)")

    plt.show()

    fig = plt.figure(figsize=(12, 8))
    ax = [fig.add_subplot(2, 3, i) for i in range(1, 7)]

    ax[0].plot(xs, dGs_F / 1.0e3)
    ax[0].plot([1.0, 0.0], gt_F.partial_gibbs / 1.0e3)
    ax[1].plot(xs, dVs_F * 1.0e6)
    ax[1].plot([1.0, 0.0], gt_F.partial_volumes * 1.0e6)
    ax[2].plot(xs, dSs_F)
    ax[2].plot([1.0, 0.0], gt_F.partial_entropies)

    ax[3].plot(
        xs,
        (dGs_F - xs * gt_F.partial_gibbs[0] - (1.0 - xs) * gt_F.partial_gibbs[1])
        / 1.0e3,
    )
    ax[4].plot(
        xs,
        (dVs_F - xs * gt_F.partial_volumes[0] - (1.0 - xs) * gt_F.partial_volumes[1])
        * 1.0e6,
    )
    ax[5].plot(
        xs,
        (
            dSs_F
            - xs * gt_F.partial_entropies[0]
            - (1.0 - xs) * gt_F.partial_entropies[1]
        ),
    )

    for i in range(3, 6):
        ax[i].plot([0.0, 1.0], [0.0, 0.0])

    for i in range(6):
        ax[i].set_xlim(0.0, 1.0)
        ax[i].set_xlabel("x(py)")

    ax[0].set_ylabel("Gibbs (kJ/mol)")
    ax[1].set_ylabel("Volume (cm$^3$/mol)")
    ax[2].set_ylabel("Entropy (J/K/mol)")
    ax[3].set_ylabel("Gibbs - $\mu(0.4)$ (kJ/mol)")
    ax[4].set_ylabel("Volume - $d\mu(0.4)/dP$ (cm$^3$/mol)")
    ax[5].set_ylabel("Entropy + $d\mu(0.4)/dT$ (J/K/mol)")

    fig.set_tight_layout(True)

    plt.show()
