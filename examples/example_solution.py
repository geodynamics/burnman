# This file is part of BurnMan - a thermoelastic and thermodynamic toolkit
# for the Earth and Planetary Sciences
# Copyright (C) 2012 - 2015 by the BurnMan team, released under the GNU
# GPL v2 or later.


"""

example_solution
----------------

This example shows how to create different solution models and output
thermodynamic and thermoelastic quantities.

There are four main types of solution currently implemented in
BurnMan:

1. Ideal solutions
2. Symmmetric solutions
3. Asymmetric solutions
4. Subregular solutions

These solutions can potentially deal with:

* Disordered endmembers (more than one element on a crystallographic site)
* Site vacancies
* More than one valence/spin state of the same element on a site

*Uses:*

* :doc:`mineral_database`
* :class:`burnman.Solution`
* :class:`burnman.SolutionModel`


*Demonstrates:*

* Different ways to define a solution
* How to set composition and state
* How to output thermodynamic and thermoelastic properties

"""
import numpy as np
import matplotlib.pyplot as plt


import burnman
from burnman import minerals
from burnman import Solution
from burnman.classes.solutionmodel import (
    IdealSolution,
    SymmetricRegularSolution,
    AsymmetricRegularSolution,
    SubregularSolution,
)


if __name__ == "__main__":
    """
    First, let's pick a starting pressure and temperature
    """
    P = 1.0e5
    T = 1000.0

    """
    We can create an instance of a solution in two ways.
    We can create a single instance using the Solution constructor
    (here's an example of an ideal pyrope-almandine garnet) ...
    """
    g1 = Solution(
        name="Ideal pyrope-almandine garnet",
        solution_model=IdealSolution(
            endmembers=[
                [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
            ]
        ),
        molar_fractions=[0.5, 0.5],
    )

    g1.set_state(P, T)
    gibbs = g1.gibbs

    """
    ... or we can create a class that derives from solution
    (so that we can create multiple instances of the same solution)
    """

    class mg_fe_garnet(burnman.Solution):
        def __init__(self, molar_fractions=None):
            self.name = "Ideal pyrope-almandine garnet"
            self.solution_model = IdealSolution(
                endmembers=[
                    [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                    [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
                ]
            )

            burnman.Solution.__init__(self, molar_fractions=molar_fractions)

    """
    Initialisation can optionally include setting the composition
    of the solution, i.e.
    """
    g1 = mg_fe_garnet([0.0, 1.0])
    g1.set_state(P, T)
    alm_gibbs = g1.gibbs

    """
    Alternatively, the composition can be set after initialisation
    """
    g1.set_composition([1.0, 0.0])
    g1.set_state(P, T)
    py_gibbs = g1.gibbs

    """
    Let's plot some interesting stuff...
    """
    comp = np.linspace(0.001, 0.999, 100)
    g1_gibbs = np.empty_like(comp)
    g1_excess_gibbs = np.empty_like(comp)
    g1_pyrope_activity = np.empty_like(comp)
    g1_almandine_activity = np.empty_like(comp)
    g1_pyrope_gamma = np.empty_like(comp)
    g1_almandine_gamma = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, c]
        g1.set_composition(molar_fractions)
        g1.set_state(P, T)
        g1_gibbs[i] = g1.gibbs
        g1_excess_gibbs[i] = g1.excess_gibbs
        g1_pyrope_activity[i] = g1.activities[0]
        g1_almandine_activity[i] = g1.activities[1]
        g1_pyrope_gamma[i] = g1.activity_coefficients[0]
        g1_almandine_gamma[i] = g1.activity_coefficients[1]

    plt.plot(
        [0.0, 1.0],
        [py_gibbs / 1000.0, alm_gibbs / 1000.0],
        "b-",
        linewidth=1.0,
        label="Mechanical mixing",
    )
    plt.plot(comp, g1_gibbs / 1000.0, "r-", linewidth=1.0, label="Ideal solution")
    plt.title("Ideal pyrope-almandine join")
    plt.ylabel("Gibbs free energy of solution (kJ/mol)")
    plt.xlabel("Almandine fraction")
    plt.legend(loc="lower right")
    plt.show()

    plt.plot(comp, g1_pyrope_activity, "r-", linewidth=1.0, label="Pyrope activity")
    plt.plot(
        comp, g1_almandine_activity, "b-", linewidth=1.0, label="Almandine activity"
    )

    plt.plot(
        comp, g1_pyrope_gamma, "r-.", linewidth=1.0, label="Pyrope activity coefficient"
    )
    plt.plot(
        comp,
        g1_almandine_gamma,
        "b--",
        linewidth=1.0,
        label="Almandine activity coefficient",
    )
    plt.title("Ideal pyrope-almandine join")
    plt.ylim(0.0, 1.01)
    plt.ylabel("Activities")
    plt.xlabel("Almandine fraction")
    plt.legend(loc="lower right")
    plt.show()

    """
    Not included in this example document are ways to create solutions
    with spin transitions, vacancies and mixed valence states. However, the
    formula parsing in BurnMan is comprehensive, so any model should be
    reproducible.
    Some formatted formulae examples follow:
    - [Fe]O: simple wuestite
    - [Fef2/3Vac1/3]O: a theoretical ferric endmember of wuestite.
      Note the distinct element name for Fe3+.
    - [Fe2/3Vac1/3]O: an alternative ferric endmember.
      Mixing between this and [Fe]O will produce a
      different configurational entropy to the previous model
      because there is no distinction between Fe2+ and Fe3+.
    - [Fels]O: Low spin wuestite. Another example illustrating the free-form
      approach to element definition.
    """

    """
    The solution corresponding to the pyrope-almandine join is not quite ideal.
    It can be well-approximated with a symmetric regular solution model
    """
    g2 = burnman.Solution(
        name="Symmetric pyrope-almandine garnet",
        solution_model=SymmetricRegularSolution(
            endmembers=[
                [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
            ],
            energy_interaction=[[2.5e3]],
        ),
    )

    g2_excess_gibbs = np.empty_like(comp)
    g2_pyrope_activity = np.empty_like(comp)
    g2_almandine_activity = np.empty_like(comp)
    g2_pyrope_gamma = np.empty_like(comp)
    g2_almandine_gamma = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, c]
        g2.set_composition(molar_fractions)
        g2.set_state(P, T)
        g2_excess_gibbs[i] = g2.excess_gibbs
        g2_pyrope_activity[i] = g2.activities[0]
        g2_almandine_activity[i] = g2.activities[1]
        g2_pyrope_gamma[i] = g2.activity_coefficients[0]
        g2_almandine_gamma[i] = g2.activity_coefficients[1]

    plt.plot(comp, g1_excess_gibbs, "r-", linewidth=1.0, label="Ideal solution")
    plt.plot(
        comp,
        g2_excess_gibbs,
        "g-",
        linewidth=1.0,
        label="Symmetric solution, 2.5 kJ/mol",
    )
    plt.title("Pyrope-almandine join (model comparison)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Almandine fraction")
    plt.legend(loc="upper left")
    plt.show()

    plt.plot(comp, g2_pyrope_activity, "r-", linewidth=1.0, label="Pyrope activity")
    plt.plot(
        comp, g2_almandine_activity, "b-", linewidth=1.0, label="Almandine activity"
    )

    plt.plot(
        comp, g2_pyrope_gamma, "r-.", linewidth=1.0, label="Pyrope activity coefficient"
    )
    plt.plot(
        comp,
        g2_almandine_gamma,
        "b--",
        linewidth=1.0,
        label="Almandine activity coefficient",
    )
    plt.title("Non-ideal pyrope-almandine join")
    plt.ylabel("Activities")
    plt.xlabel("Almandine fraction")
    plt.legend(loc="lower right")
    plt.show()

    """
    Adding more endmembers is very straightforward.
    Interaction terms must be added in the order
    [[12,13,14,...],[23,24,...],[34,...],...]
    Here, the new endmember is majorite, illustrating the addition of
    endmembers with multiple occupancy on a single site
    (here Mg and Si to replace Al)

    Here we also illustrate the addition of excess entropy and volume terms.
    V and S excess lists are optional parameters to the solution
    initialisation.
    In the initialisation, the order of interaction terms is H, V, S
    (S excesses are least commonly constrained).

    An important note: In Holland and Powell datasets, the "DQF" terms
    are equivalent to the excesses in burnman with the *exception* of the
    temperature correction, which is the negative of the entropy correction.
    i.e.
    DQF[0] = H_excess
    DQF[1] = -S_excess
    DQF[2] = V_excess
    """

    g3 = Solution(
        name="Symmetric pyrope-almandine-majorite garnet",
        solution_model=SymmetricRegularSolution(
            endmembers=[
                [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.maj(), "[Mg]3[Mg1/2Si1/2]2Si3O12"],
            ],
            energy_interaction=[[2.5e3, 0.0e3], [10.0e3]],
            entropy_interaction=[[0.0e3, 0.0e3], [0.0e3]],
            volume_interaction=[[0.0e3, 0.0e3], [0.0e3]],
        ),
    )

    g3_configurational_entropy = np.empty_like(comp)
    g3_excess_entropy = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, 0.0, c]
        g3.set_composition(molar_fractions)
        g3.set_state(P, T)
        Sconf = g3.solution_model.configurational_entropy(molar_fractions)
        Gex = g3.solution_model._ideal_excess_partial_gibbs(T, molar_fractions)
        g3_configurational_entropy[i] = Sconf
        g3_excess_entropy[i] = -np.dot(molar_fractions, Gex) / T

    plt.plot(
        comp,
        g3_configurational_entropy,
        "g-",
        linewidth=1.0,
        label="Configurational entropy",
    )
    plt.plot(comp, g3_excess_entropy, "r-", linewidth=1.0, label="Excess entropy")
    plt.title("Pyrope-majorite join")
    plt.ylabel("Entropy (J/K/mol)")
    plt.xlabel("Majorite fraction")
    plt.legend(loc="upper left")
    plt.show()

    """
    The addition of grossular garnet (Ca-Al garnet) illustrates
    the use of asymmetric solution models.
    This model is also found in the HP_2011_ds62 database
    """

    g4 = Solution(
        name="garnet",
        solution_model=AsymmetricRegularSolution(
            endmembers=[
                [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.gr(), "[Ca]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.andr(), "[Ca]3[Fe]2Si3O12"],
            ],
            alphas=[1.0, 1.0, 2.7, 2.7],
            energy_interaction=[[2.5e3, 31.0e3, 53.2e3], [5.0e3, 37.24e3], [2.0e3]],
        ),
    )

    g4_excess_gibbs_400 = np.empty_like(comp)
    g4_excess_gibbs_800 = np.empty_like(comp)
    g4_excess_gibbs_1200 = np.empty_like(comp)

    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, 0.0, c, 0.0]
        g4.set_composition(molar_fractions)
        g4.set_state(P, 400.0)
        g4_excess_gibbs_400[i] = g4.excess_gibbs
        g4.set_state(P, 800.0)
        g4_excess_gibbs_800[i] = g4.excess_gibbs
        g4.set_state(P, 1200.0)
        g4_excess_gibbs_1200[i] = g4.excess_gibbs

    plt.plot(comp, g4_excess_gibbs_400, "r-", linewidth=1.0, label="400 K")
    plt.plot(comp, g4_excess_gibbs_800, "g-", linewidth=1.0, label="800 K")
    plt.plot(comp, g4_excess_gibbs_1200, "b-", linewidth=1.0, label="1200 K")
    plt.title("Pyrope-grossular join (asymmetric model)")
    plt.ylabel("Excess gibbs free energy of solution (J/mol)")
    plt.xlabel("Grossular fraction")
    plt.legend(loc="lower left")
    plt.show()

    """
    The subregular solution model (Helffrich and Wood, 1989)
    provides a more flexible way of constructing an asymmetric
    model. Here's a garnet model from Ganguly et al., 1996:

    (N.B. Published excesses are on a 4-oxygen (1-cation) basis,
    so we need to multiply each by 3.)
    """

    def mult(x, n):
        return [[[v * n for v in i] for i in j] for j in x]

    g5 = Solution(
        name="Subregular pyrope-almandine-grossular " "garnet (Ganguly et al., 1996)",
        solution_model=SubregularSolution(
            endmembers=[
                [minerals.HP_2011_ds62.py(), "[Mg]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.alm(), "[Fe]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.gr(), "[Ca]3[Al]2Si3O12"],
                [minerals.HP_2011_ds62.spss(), "[Mn]3[Al]2Si3O12"],
            ],
            energy_interaction=mult(
                [
                    [[2117.0, 695.0], [9834.0, 21627.0], [12083.0, 12083.0]],
                    [[6773.0, 873.0], [539.0, 539.0]],
                    [[0.0, 0.0]],
                ],
                3.0,
            ),
            volume_interaction=mult(
                [
                    [[0.07e-5, 0.0], [0.058e-5, 0.012e-5], [0.04e-5, 0.03e-5]],
                    [[0.03e-5, 0.0], [0.04e-5, 0.01e-5]],
                    [[0.0, 0.0]],
                ],
                3.0,
            ),
            entropy_interaction=mult(
                [
                    [[0.0, 0.0], [5.78, 5.78], [7.67, 7.67]],
                    [[1.69, 1.69], [0.0, 0.0]],
                    [[0.0, 0.0]],
                ],
                3.0,
            ),
        ),
    )

    g5_excess_enthalpy = np.empty_like(comp)
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, 0.0, c, 0.0]
        g5.set_composition(molar_fractions)
        g5.set_state(1.0e5, 298.15)
        g5_excess_enthalpy[i] = g5.excess_gibbs + 298.15 * g5.excess_entropy

    plt.plot(
        comp,
        g5_excess_enthalpy / 3.0,
        "r-",
        linewidth=1.0,
        label="Py-Gr excess enthalpy (J/cation-mole)",
    )

    plt.title("Asymmetric py-gr join (Ganguly et al., 1996; Figure 5)")
    plt.ylabel("Excess enthalpy of solution (J/cation-mol)")
    plt.xlabel("XCa")
    plt.legend(loc="lower left")
    plt.show()
