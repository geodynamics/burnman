# Benchmarks for the solid solution class
import burnman
from burnman import minerals
from burnman.classes.solutionmodel import (
    SymmetricRegularSolution,
    AsymmetricRegularSolution,
    SubregularSolution,
)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg


"""
Solvus shapes (a proxy for Gibbs free energy checking
"""

# van Laar parameter
# Figure 2a of Holland and Powell, 2003

# Temperature dependence
# Figure 2b of Holland and Powell, 2003

# A specific solvus example: sanidine-high albite
# Includes asymmetry and pressure, temperature dependence
# Figure 3 of Holland and Powell, 2003


"""
Excess properties
"""
# Configurational entropy
# Navrotsky and Kleppa, 1967


class o_d_spinel(burnman.SolidSolution):
    def __init__(self):
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [minerals.HP_2011_ds62.sp(), "[Mg][Al]2O4"],
                [minerals.HP_2011_ds62.sp(), "[Al][Mg1/2Al1/2]2O4"],
            ],
            energy_interaction=[[0.0]],
        )

        burnman.SolidSolution.__init__(self)


comp = np.linspace(0.001, 0.999, 100)
sp = o_d_spinel()
sp_entropies = np.empty_like(comp)
sp_entropies_NK1967 = np.empty_like(comp)
for i, c in enumerate(comp):
    molar_fractions = [1.0 - c, c]
    sp.set_composition(np.array(molar_fractions))
    sp.set_state(1e5, 298.15)
    sp_entropies[i] = sp.solution_model.configurational_entropy(molar_fractions)
    sp_entropies_NK1967[i] = -8.3145 * (
        c * np.log(c)
        + (1.0 - c) * np.log(1.0 - c)
        + c * np.log(c / 2.0)
        + (2.0 - c) * np.log(1.0 - c / 2.0)
    )  # eq. 7 in Navrotsky and Kleppa, 1967.

# fig1 = mpimg.imread('configurational_entropy.png')  # Uncomment these two lines if you want to overlay the plot on a screengrab from SLB2011
# plt.imshow(fig1, extent=[0.0, 1.0,0.,17.0], aspect='auto')
plt.plot(comp, sp_entropies_NK1967, "b-", linewidth=3.0)
plt.plot(comp, sp_entropies, "r--", linewidth=3.0)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 17.0)
plt.ylabel("Configurational entropy of solution (J/K/mol)")
plt.xlabel("fraction inverse spinel")
plt.show()

# Configurational entropy
# Figure 3b of Stixrude and Lithgow-Bertelloni, 2011


class orthopyroxene_red(burnman.SolidSolution):
    def __init__(self):
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [minerals.SLB_2011.enstatite(), "Mg[Mg][Si]SiO6"],
                [minerals.SLB_2011.mg_tschermaks(), "Mg[Al][Al]SiO6"],
            ],
            energy_interaction=[[0.0]],
        )

        burnman.SolidSolution.__init__(self)


class orthopyroxene_blue(burnman.SolidSolution):
    def __init__(self):
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [minerals.SLB_2011.enstatite(), "Mg[Mg]Si2O6"],
                [minerals.SLB_2011.mg_tschermaks(), "Mg[Al]AlSiO6"],
            ],
            energy_interaction=[[0.0]],
        )

        burnman.SolidSolution.__init__(self)


class orthopyroxene_long_dashed(burnman.SolidSolution):
    def __init__(self):
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [minerals.SLB_2011.enstatite(), "Mg[Mg]Si2O6"],
                [minerals.SLB_2011.mg_tschermaks(), "[Mg1/2Al1/2]2AlSiO6"],
            ],
            energy_interaction=[[10.0e3]],
        )

        burnman.SolidSolution.__init__(self)


class orthopyroxene_short_dashed(burnman.SolidSolution):
    def __init__(self):
        self.name = "orthopyroxene"
        self.solution_model = SymmetricRegularSolution(
            endmembers=[
                [minerals.SLB_2011.enstatite(), "Mg[Mg][Si]2O6"],
                [minerals.SLB_2011.mg_tschermaks(), "Mg[Al][Al1/2Si1/2]2O6"],
            ],
            energy_interaction=[[0.0]],
        )

        burnman.SolidSolution.__init__(self)


comp = np.linspace(0, 1.0, 100)
opx_models = [
    orthopyroxene_red(),
    orthopyroxene_blue(),
    orthopyroxene_long_dashed(),
    orthopyroxene_short_dashed(),
]
opx_entropies = [np.empty_like(comp) for model in opx_models]
for idx, model in enumerate(opx_models):
    for i, c in enumerate(comp):
        molar_fractions = [1.0 - c, c]
        model.set_composition(np.array(molar_fractions))
        model.set_state(0.0, 0.0)
        opx_entropies[idx][i] = model.solution_model.configurational_entropy(
            molar_fractions
        )

fig1 = mpimg.imread("configurational_entropy.png")
# Uncomment these two lines if you want to overlay the plot
# on a screengrab from SLB2011
plt.imshow(fig1, extent=[0.0, 1.0, 0.0, 17.0], aspect="auto")
plt.plot(comp, opx_entropies[0], "r--", linewidth=3.0)
plt.plot(comp, opx_entropies[1], "b--", linewidth=3.0)
plt.plot(comp, opx_entropies[2], "g--", linewidth=3.0)
plt.plot(comp, opx_entropies[3], "g-.", linewidth=3.0)
plt.xlim(0.0, 1.0)
plt.ylim(0.0, 17.0)
plt.ylabel("Configurational entropy of solution (J/K/mol)")
plt.xlabel("cats fraction")
plt.show()

# Excess volume of solution

# Excess energy of solution
# Figure 5 of Stixrude and Lithgow-Bertelloni, 2011


class clinopyroxene(burnman.SolidSolution):
    def __init__(self):
        self.name = "clinopyroxene"
        self.solution_model = AsymmetricRegularSolution(
            endmembers=[
                [minerals.SLB_2011.diopside(), "[Ca][Mg][Si]2O6"],
                [minerals.SLB_2011.ca_tschermaks(), "[Ca][Al][Si1/2Al1/2]2O6"],
            ],
            energy_interaction=[[26.0e3]],
            alphas=[1.0, 3.5],
        )

        burnman.SolidSolution.__init__(self)


cpx = clinopyroxene()

comp = np.linspace(0, 1.0, 100)
gibbs = np.empty_like(comp)

for i, c in enumerate(comp):
    cpx.set_composition(np.array([1.0 - c, c]))
    cpx.set_state(0.0, 0.0)
    gibbs[i] = cpx.excess_gibbs


fig1 = mpimg.imread("dicats.png")
# Uncomment these two lines if you want to overlay the plot
# on a screengrab from SLB2011
plt.imshow(fig1, extent=[0.0, 1.0, -2.0, 8.0], aspect="auto")

plt.plot(comp, gibbs / 1000.0, "b--", linewidth=3.0)
plt.xlim(0.0, 1.0)
plt.ylim(-2.0, 8.0)
plt.ylabel("Excess energy of solution (kJ/mol)")
plt.xlabel("cats fraction")
plt.show()


# Subregular model from Ganguly et al., 1996

mult = lambda x, n: [[[v * n for v in i] for i in j] for j in x]

g5 = burnman.SolidSolution(
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

fig1 = mpimg.imread("figures/Ganguly_1996_Fig_5.png")
plt.imshow(fig1, extent=[0.0, 1.0, -1000.0, 5000.0], aspect="auto")

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
